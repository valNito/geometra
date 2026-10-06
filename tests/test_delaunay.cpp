#include "geometra/delaunay.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <list>
#include <numbers>
#include <random>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "geometra/convex_hull.hpp"
#include "geometra/predicates.hpp"
#include "geotest.hpp"

namespace {

using geo::delaunay_triangulation;
using geo::incircle;
using geo::no_neighbor;
using geo::Orientation;
using geo::orientation;
using geo::Point2D;
using geo::Triangle;
using geo::Triangulation;
using Points = std::vector<Point2D>;
using Triple = std::array<std::size_t, 3>;

bool on_closed_segment(Point2D p, Point2D q, Point2D r) {
    if (orientation(q, r, p) != Orientation::collinear) return false;
    const auto [lo, hi] = std::minmax(q, r);
    return lo <= p && p <= hi;  // sobre la recta, el orden lexicográfico es monótono
}

/// Para cada punto distinto de la entrada, el menor índice en que aparece.
std::vector<std::size_t> representatives(const Points& pts) {
    std::vector<std::size_t> ids(pts.size());
    for (std::size_t i = 0; i < ids.size(); ++i) ids[i] = i;
    std::ranges::sort(ids, [&](std::size_t a, std::size_t b) {
        return pts[a] != pts[b] ? pts[a] < pts[b] : a < b;
    });
    const auto dup = std::ranges::unique(ids, [&](std::size_t a, std::size_t b) {
        return pts[a] == pts[b];
    });
    ids.erase(dup.begin(), dup.end());
    return ids;
}

/// El triángulo vecino de `t` por su arista i debe tener esa misma arista en sentido
/// contrario y apuntar de vuelta a `t`. Devuelve el vértice del vecino opuesto a ella.
std::size_t check_neighbor(const Triangulation& tri, std::size_t t, std::size_t i) {
    const Triangle& tr = tri.triangles[t];
    const std::size_t a = tr.v[(i + 1) % 3];
    const std::size_t b = tr.v[(i + 2) % 3];
    const Triangle& other = tri.triangles[tr.neighbor[i]];
    for (std::size_t j = 0; j < 3; ++j) {
        if (other.v[(j + 1) % 3] == b && other.v[(j + 2) % 3] == a) {
            CHECK_EQ(other.neighbor[j], t);
            return other.v[j];
        }
    }
    CHECK_MSG(false, "el vecino no comparte la arista");
    return a;
}

/// Comprueba todo el contrato de delaunay_triangulation sobre `input`. Con `global`
/// además verifica la propiedad de Delaunay contra todos los puntos (O(n·T)); sin él,
/// solo la local entre triángulos vecinos, que en una triangulación válida es equivalente.
void check_delaunay(const Points& input, const Triangulation& tri, bool global = true) {
    CHECK(tri.points == input);
    const std::size_t n = input.size();
    const std::vector<std::size_t> reps = representatives(input);
    const Points hull = geo::convex_hull(input);
    const std::size_t count = tri.triangles.size();

    if (hull.size() < 3) {  // menos de 3 puntos distintos o todos colineales
        CHECK_EQ(count, std::size_t{0});
        return;
    }

    // Euler: T = 2n - h - 2, con h los puntos distintos sobre el borde de la envolvente.
    std::size_t on_boundary = 0;
    for (const std::size_t id : reps) {
        for (std::size_t i = 0; i < hull.size(); ++i) {
            if (on_closed_segment(input[id], hull[i], hull[(i + 1) % hull.size()])) {
                ++on_boundary;
                break;
            }
        }
    }
    CHECK_EQ(count, 2 * reps.size() - on_boundary - 2);

    std::vector<char> used(n, 0);
    std::set<std::pair<std::size_t, std::size_t>> directed;
    std::size_t boundary_edges = 0;
    std::size_t bad_shape = 0;
    std::size_t outside_boundary = 0;
    std::size_t local_violations = 0;

    for (std::size_t t = 0; t < count; ++t) {
        const Triangle& tr = tri.triangles[t];
        REQUIRE(tr.v[0] < n && tr.v[1] < n && tr.v[2] < n);
        const Point2D a = input[tr.v[0]];
        const Point2D b = input[tr.v[1]];
        const Point2D c = input[tr.v[2]];
        if (!(tr.v[0] < tr.v[1] && tr.v[0] < tr.v[2])) ++bad_shape;
        if (orientation(a, b, c) != Orientation::counterclockwise) ++bad_shape;
        if (t > 0 && !(tri.triangles[t - 1].v < tr.v)) ++bad_shape;

        for (std::size_t i = 0; i < 3; ++i) {
            used[tr.v[i]] = 1;
            const std::size_t from = tr.v[(i + 1) % 3];
            const std::size_t to = tr.v[(i + 2) % 3];
            if (!directed.insert({from, to}).second) ++bad_shape;  // arista dirigida repetida

            if (tr.neighbor[i] == no_neighbor) {
                // Arista de borde: ningún punto puede quedar a su derecha.
                ++boundary_edges;
                for (const std::size_t id : reps) {
                    if (orientation(input[from], input[to], input[id]) == Orientation::clockwise) {
                        ++outside_boundary;
                    }
                }
            } else {
                REQUIRE(tr.neighbor[i] < count);
                const std::size_t opposite = check_neighbor(tri, t, i);
                if (incircle(a, b, c, input[opposite]) > 0.0) ++local_violations;
            }
        }
    }
    CHECK_EQ(bad_shape, std::size_t{0});
    CHECK_EQ(boundary_edges, on_boundary);
    CHECK_EQ(outside_boundary, std::size_t{0});
    CHECK_EQ(local_violations, std::size_t{0});

    // Exactamente los representantes son vértices: ni repetidos ni puntos olvidados.
    std::vector<char> is_rep(n, 0);
    for (const std::size_t id : reps) is_rep[id] = 1;
    CHECK(used == is_rep);

    if (global) {
        std::size_t violations = 0;
        for (const Triangle& tr : tri.triangles) {
            for (const std::size_t id : reps) {
                if (incircle(input[tr.v[0]], input[tr.v[1]], input[tr.v[2]], input[id]) > 0.0) {
                    ++violations;
                }
            }
        }
        CHECK_EQ(violations, std::size_t{0});
    }
}

/// Triángulos como ternas de puntos, rotadas para empezar por el menor y ordenadas:
/// permite comparar triangulaciones de entradas con distinto orden o con repetidos.
std::vector<std::array<Point2D, 3>> geometric(const Triangulation& tri) {
    std::vector<std::array<Point2D, 3>> out;
    for (const Triangle& tr : tri.triangles) {
        std::array<Point2D, 3> p{tri.points[tr.v[0]], tri.points[tr.v[1]], tri.points[tr.v[2]]};
        std::ranges::rotate(p, std::ranges::min_element(p));
        out.push_back(p);
    }
    std::ranges::sort(out);
    return out;
}

std::vector<Triple> vertex_triples(const Triangulation& tri) {
    std::vector<Triple> out;
    for (const Triangle& tr : tri.triangles) out.push_back(tr.v);
    return out;
}

/// Triangulación de Delaunay por definición, en O(n^4): todas las ternas cuya
/// circunferencia deja estrictamente fuera al resto, empezando por el menor índice. Solo
/// vale en posición general (sin tres puntos colineales ni cuatro concíclicos); si no la
/// hay, devuelve false.
bool brute_force(const Points& pts, std::vector<Triple>& out) {
    const std::size_t n = pts.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            for (std::size_t k = j + 1; k < n; ++k) {
                const Orientation o = orientation(pts[i], pts[j], pts[k]);
                if (o == Orientation::collinear) return false;
                const Triple t = o == Orientation::counterclockwise ? Triple{i, j, k}
                                                                    : Triple{i, k, j};
                bool empty = true;
                for (std::size_t l = 0; l < n; ++l) {
                    if (l == i || l == j || l == k) continue;
                    const double s = incircle(pts[t[0]], pts[t[1]], pts[t[2]], pts[l]);
                    if (s == 0.0) return false;
                    if (s > 0.0) empty = false;
                }
                if (empty) out.push_back(t);
            }
        }
    }
    std::ranges::sort(out);
    return true;
}

Points grid(int k, Point2D origin = {0, 0}) {
    Points pts;
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            pts.push_back(origin + Point2D{static_cast<double>(i), static_cast<double>(j)});
        }
    }
    return pts;
}

Points random_points(std::mt19937_64& rng, std::size_t n, double lo, double hi) {
    std::uniform_real_distribution<double> u(lo, hi);
    Points pts(n);
    for (Point2D& p : pts) p = {u(rng), u(rng)};
    return pts;
}

Points on_circle(int n, Point2D center, double radius) {
    Points pts;
    for (int i = 0; i < n; ++i) {
        const double a = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n);
        pts.push_back(center + Point2D{radius * std::cos(a), radius * std::sin(a)});
    }
    return pts;
}

}  // namespace

TEST_CASE("delaunay: menos de tres puntos distintos no dan triángulos") {
    for (const Points& pts : {Points{}, Points{{1, 2}}, Points{{0, 0}, {1, 1}},
                              Points{{3, 3}, {3, 3}, {3, 3}, {3, 3}},
                              Points{{0, 0}, {1, 0}, {0, 0}, {1, 0}}}) {
        const Triangulation tri = delaunay_triangulation(pts);
        CHECK(tri.points == pts);
        CHECK(tri.triangles.empty());
    }
}

TEST_CASE("delaunay: un triángulo, en cualquier orden de entrada") {
    const Triangulation ccw = delaunay_triangulation({{0, 0}, {1, 0}, {0, 1}});
    REQUIRE(ccw.triangles.size() == 1);
    CHECK(ccw.triangles[0].v == (Triple{0, 1, 2}));
    CHECK(ccw.triangles[0].neighbor == (Triple{no_neighbor, no_neighbor, no_neighbor}));

    // Entrada en sentido horario: se reordena a antihorario empezando por el menor índice.
    const Triangulation cw = delaunay_triangulation({{0, 0}, {0, 1}, {1, 0}});
    REQUIRE(cw.triangles.size() == 1);
    CHECK(cw.triangles[0].v == (Triple{0, 2, 1}));

    Points pts{{5, -1}, {2, 7}, {-3, 0}};
    std::ranges::sort(pts);
    do {
        check_delaunay(pts, delaunay_triangulation(pts));
    } while (std::ranges::next_permutation(pts).found);
}

TEST_CASE("delaunay: colineales, con o sin repetidos, no dan triángulos") {
    const Points cases[] = {
        {{0, 0}, {1, 0}, {2, 0}},
        {{3, 5}, {3, -1}, {3, 2}, {3, 5}, {3, 0}},
        {{0, 0}, {4, 4}, {1, 1}, {3, 3}, {2, 2}, {1, 1}, {0, 0}},
        {{1e15, 1e15 / 4}, {-2e15, -2e15 / 4}, {0.5, 0.125}},
    };
    for (const Points& pts : cases) {
        CHECK(delaunay_triangulation(pts).triangles.empty());
    }
}

TEST_CASE("delaunay: cuadrado (cuatro concíclicos) y cuadrado con centro") {
    const Points square{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    const Triangulation sq = delaunay_triangulation(square);
    REQUIRE(sq.triangles.size() == 2);
    check_delaunay(square, sq);
    for (const Triangle& t : sq.triangles) {
        CHECK_EQ(std::ranges::count(t.neighbor, no_neighbor), 2);
    }

    Points centered = square;
    centered.push_back({0.5, 0.5});
    const Triangulation c = delaunay_triangulation(centered);
    REQUIRE(c.triangles.size() == 4);
    check_delaunay(centered, c);
    for (const Triangle& t : c.triangles) {
        CHECK(std::ranges::find(t.v, std::size_t{4}) != t.v.end());  // todos tocan el centro
        CHECK_EQ(std::ranges::count(t.neighbor, no_neighbor), 1);
    }
}

TEST_CASE("delaunay: los puntos repetidos usan el de menor índice") {
    const Points pts{{0, 0}, {4, 0}, {0, 4}, {4, 0}, {0, 0}, {1, 1}, {4, 4}, {1, 1}};
    const Triangulation tri = delaunay_triangulation(pts);
    check_delaunay(pts, tri);
    for (const Triangle& t : tri.triangles) {
        for (const std::size_t v : t.v) CHECK(v != 3 && v != 4 && v != 7);
    }
}

TEST_CASE("delaunay: los colineales sobre el borde de la envolvente son vértices") {
    Points pts;
    for (int i = 0; i <= 6; ++i) {
        const auto t = static_cast<double>(i);
        pts.push_back({t, 0});
        pts.push_back({6, t});
        pts.push_back({6 - t, 6});
        pts.push_back({0, 6 - t});
    }
    pts.push_back({2.5, 3.5});
    pts.push_back({4, 1.5});
    check_delaunay(pts, delaunay_triangulation(pts));
}

TEST_CASE("delaunay: rejillas enteras (todo concíclico), ordenadas y desordenadas") {
    // En una rejilla unitaria de k x k, toda triangulación de Delaunay parte cada celda en
    // dos: 2(k-1)² triángulos de área exactamente 1/2, que con coordenadas enteras se
    // comprueba sin error de redondeo.
    std::mt19937_64 rng(17);
    for (int k = 2; k <= 12; ++k) {
        Points pts = grid(k);
        for (int round = 0; round < 3; ++round) {
            const Triangulation tri = delaunay_triangulation(pts);
            const auto cells = static_cast<std::size_t>((k - 1) * (k - 1));
            CHECK_EQ(tri.triangles.size(), 2 * cells);
            check_delaunay(pts, tri);
            std::size_t wrong_area = 0;
            for (const Triangle& t : tri.triangles) {
                if (geo::orient2d(pts[t.v[0]], pts[t.v[1]], pts[t.v[2]]) != 1.0) ++wrong_area;
            }
            CHECK_EQ(wrong_area, std::size_t{0});
            std::ranges::shuffle(pts, rng);
        }
    }
}

TEST_CASE("delaunay: circunferencia de 1000 puntos (todos en el borde)") {
    const Points pts = on_circle(1000, {0, 0}, 1.0);
    const Triangulation tri = delaunay_triangulation(pts);
    CHECK_EQ(tri.triangles.size(), std::size_t{998});
    check_delaunay(pts, tri);
}

TEST_CASE("delaunay: coincide con la fuerza bruta en conjuntos pequeños") {
    std::mt19937_64 rng(8);
    std::uniform_int_distribution<std::size_t> size(3, 9);
    int compared = 0;
    for (int trial = 0; trial < 1500; ++trial) {
        const Points pts = random_points(rng, size(rng), 0.0, 1.0);
        std::vector<Triple> expected;
        if (!brute_force(pts, expected)) continue;
        ++compared;
        const Triangulation tri = delaunay_triangulation(pts);
        CHECK(vertex_triples(tri) == expected);
        check_delaunay(pts, tri);
    }
    CHECK(compared > 1000);
}

TEST_CASE("delaunay: propiedades en conjuntos aleatorios grandes") {
    std::mt19937_64 rng(2026);
    const Points uniform = random_points(rng, 2000, -1.0, 1.0);
    check_delaunay(uniform, delaunay_triangulation(uniform));

    std::normal_distribution<double> g(0.0, 1.0);
    Points gaussian(30000);
    for (Point2D& p : gaussian) p = {g(rng), g(rng)};
    check_delaunay(gaussian, delaunay_triangulation(gaussian), false);

    // Muchos repetidos y concíclicos: enteros en una rejilla chica.
    std::uniform_int_distribution<int> small(0, 40);
    Points lattice(5000);
    for (Point2D& p : lattice) {
        p = {static_cast<double>(small(rng)), static_cast<double>(small(rng))};
    }
    check_delaunay(lattice, delaunay_triangulation(lattice), false);
}

TEST_CASE("delaunay: invariante ante permutaciones y duplicación de la entrada") {
    std::mt19937_64 rng(77);
    const Points inputs[] = {grid(9), random_points(rng, 500, 0.0, 1.0)};
    for (const Points& pts : inputs) {
        const auto expected = geometric(delaunay_triangulation(pts));

        Points shuffled = pts;
        std::ranges::shuffle(shuffled, rng);
        CHECK(geometric(delaunay_triangulation(shuffled)) == expected);

        const Points reversed(pts.rbegin(), pts.rend());
        CHECK(geometric(delaunay_triangulation(reversed)) == expected);

        Points doubled = pts;
        doubled.insert(doubled.end(), shuffled.begin(), shuffled.end());
        const Triangulation tri = delaunay_triangulation(doubled);
        CHECK(geometric(tri) == expected);
        check_delaunay(doubled, tri);
    }
}

TEST_CASE("delaunay: robustez con puntos casi concíclicos y casi colineales") {
    // Circunferencia de radio 1 centrada en (1e9, 1e9): el redondeo de las coordenadas
    // deja los puntos a distancias minúsculas de la circunferencia, y muchas decisiones
    // caen en la evaluación exacta.
    const Points circle = on_circle(400, {1e9, 1e9}, 1.0);
    check_delaunay(circle, delaunay_triangulation(circle));

    // Puntos sobre y = x/3 con x entre 1e15 y 2e15: la división redondea y en el último
    // bit, así que quedan a distancias minúsculas de la recta.
    Points line;
    for (int i = 0; i < 500; ++i) {
        const double x = 1e15 + 1e12 * static_cast<double>(i);
        line.push_back({x, x / 3.0});
    }
    check_delaunay(line, delaunay_triangulation(line));
    line.push_back({1.5e15, 1e15});
    line.push_back({1.5e15, -1e15});
    check_delaunay(line, delaunay_triangulation(line));

    // Rejilla entera trasladada a 2^40: concíclicos exactos lejos del origen.
    const Points far = grid(10, {std::ldexp(1.0, 40), -std::ldexp(1.0, 40)});
    const Triangulation tri = delaunay_triangulation(far);
    CHECK_EQ(tri.triangles.size(), std::size_t{162});
    check_delaunay(far, tri);

    // Puntos casi alineados entre dos extremos aleatorios.
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> coord(-1e9, 1e9);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int trial = 0; trial < 100; ++trial) {
        const Point2D a{coord(rng), coord(rng)};
        const Point2D b{coord(rng), coord(rng)};
        Points pts;
        for (int i = 0; i < 40; ++i) {
            const double t = unit(rng);
            pts.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
        }
        check_delaunay(pts, delaunay_triangulation(pts));
    }
}

TEST_CASE("delaunay: lanza std::invalid_argument con coordenadas no finitas") {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(delaunay_triangulation({{0, 0}, {nan, 1}, {1, 0}}), std::invalid_argument);
    CHECK_THROWS_AS(delaunay_triangulation({{0, 0}, {1, inf}}), std::invalid_argument);
    CHECK_THROWS_AS(delaunay_triangulation({{-inf, 0}}), std::invalid_argument);
}

TEST_CASE("delaunay: acepta vector, array, initializer_list, list, span y vistas de ranges") {
    const Points vec{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}};
    const auto expected = geometric(delaunay_triangulation(vec));
    REQUIRE(expected.size() == 4);

    const std::array<Point2D, 5> arr{{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}}};
    CHECK(geometric(delaunay_triangulation(arr)) == expected);
    CHECK(geometric(delaunay_triangulation(std::span<const Point2D>{vec})) == expected);
    CHECK(geometric(delaunay_triangulation({{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}})) ==
          expected);

    const std::list<Point2D> lst(vec.begin(), vec.end());
    CHECK(geometric(delaunay_triangulation(lst)) == expected);

    const auto view = vec | std::views::transform([](Point2D p) { return p; });
    CHECK(geometric(delaunay_triangulation(view)) == expected);
}

GEOTEST_MAIN

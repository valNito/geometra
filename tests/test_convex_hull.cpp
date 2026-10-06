#include "geometra/convex_hull.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <list>
#include <numbers>
#include <random>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "geometra/polygon.hpp"
#include "geometra/predicates.hpp"
#include "geotest.hpp"

namespace {

using geo::convex_hull;
using geo::Orientation;
using geo::orientation;
using geo::Point2D;
using Points = std::vector<Point2D>;

Points sorted_unique(Points pts) {
    std::ranges::sort(pts);
    pts.erase(std::ranges::unique(pts).begin(), pts.end());
    return pts;
}

bool on_closed_segment(Point2D p, Point2D q, Point2D r) {
    if (orientation(q, r, p) != Orientation::collinear) return false;
    const auto [lo, hi] = std::minmax(q, r);
    return lo <= p && p <= hi;  // sobre la recta, el orden lexicográfico es monótono
}

bool in_closed_triangle(Point2D p, Point2D q, Point2D r, Point2D s) {
    const Orientation o = orientation(q, r, s);
    if (o == Orientation::collinear) {
        return on_closed_segment(p, q, r) || on_closed_segment(p, r, s) ||
               on_closed_segment(p, q, s);
    }
    const auto same_side = [&](Point2D u, Point2D v) {
        const Orientation w = orientation(u, v, p);
        return w == o || w == Orientation::collinear;
    };
    return same_side(q, r) && same_side(r, s) && same_side(s, q);
}

// Vértices extremos por definición (Carathéodory en 2D): p es vértice de la envolvente
// si y solo si no pertenece a ningún segmento ni triángulo formado por los demás puntos.
// O(n^4): solo para conjuntos pequeños. Devuelve los vértices en orden lexicográfico.
Points brute_force_hull_vertices(const Points& input) {
    const Points pts = sorted_unique(input);
    Points result;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        Points others;
        for (std::size_t j = 0; j < pts.size(); ++j) {
            if (j != i) others.push_back(pts[j]);
        }
        bool extreme = true;
        for (std::size_t j = 0; j < others.size() && extreme; ++j) {
            for (std::size_t k = j + 1; k < others.size() && extreme; ++k) {
                if (on_closed_segment(pts[i], others[j], others[k])) extreme = false;
                for (std::size_t l = k + 1; l < others.size() && extreme; ++l) {
                    if (in_closed_triangle(pts[i], others[j], others[k], others[l])) {
                        extreme = false;
                    }
                }
            }
        }
        if (extreme) result.push_back(pts[i]);
    }
    return result;
}

bool strictly_convex_ccw(const Points& h) {
    const std::size_t n = h.size();
    if (n < 3) return true;
    for (std::size_t i = 0; i < n; ++i) {
        if (orientation(h[i], h[(i + 1) % n], h[(i + 2) % n]) != Orientation::counterclockwise) {
            return false;
        }
    }
    return true;
}

bool contains_all(const Points& h, const Points& pts) {
    if (h.empty()) return pts.empty();
    if (h.size() == 1) return std::ranges::all_of(pts, [&](Point2D p) { return p == h[0]; });
    if (h.size() == 2) {
        return std::ranges::all_of(pts, [&](Point2D p) { return on_closed_segment(p, h[0], h[1]); });
    }
    for (const Point2D p : pts) {
        for (std::size_t i = 0; i < h.size(); ++i) {
            if (orientation(h[i], h[(i + 1) % h.size()], p) == Orientation::clockwise) return false;
        }
    }
    return true;
}

bool subset_of_input(const Points& h, const Points& pts) {
    const Points s = sorted_unique(pts);
    return std::ranges::all_of(h, [&](Point2D v) { return std::ranges::binary_search(s, v); });
}

bool has_duplicates(const Points& h) { return sorted_unique(h).size() != h.size(); }

// Propiedades que debe cumplir toda envolvente devuelta por la librería.
void check_invariants(const Points& pts, const Points& h) {
    CHECK(contains_all(h, pts));
    CHECK(subset_of_input(h, pts));
    CHECK(strictly_convex_ccw(h));
    CHECK(!has_duplicates(h));
    if (!pts.empty()) CHECK_EQ(h.front(), *std::ranges::min_element(pts));
    // La orientación CCW ya la verifica strictly_convex_ccw con el predicado exacto.
    // signed_area es una magnitud en double: en envolventes casi degeneradas puede
    // redondear a <= 0 aunque la envolvente sea correcta, así que no sirve de invariante.
}

bool same_set(const Points& a, const Points& b) { return sorted_unique(a) == sorted_unique(b); }

Points random_points(std::mt19937_64& rng, std::size_t n, double lo, double hi) {
    std::uniform_real_distribution<double> u(lo, hi);
    Points pts(n);
    for (Point2D& p : pts) p = {u(rng), u(rng)};
    return pts;
}

}  // namespace

TEST_CASE("convex_hull: conjunto vacío") {
    const Points empty;
    CHECK(convex_hull(empty).empty());
    CHECK(convex_hull(std::span<const Point2D>{}).empty());
    CHECK(convex_hull(std::initializer_list<Point2D>{}).empty());
}

TEST_CASE("convex_hull: un punto") {
    const Points one{{3.0, -1.0}};
    CHECK_EQ(convex_hull(one), one);
}

TEST_CASE("convex_hull: dos puntos distintos") {
    const Points two{{2.0, 1.0}, {0.0, 0.0}};
    const Points h = convex_hull(two);
    CHECK_EQ(h, (Points{{0.0, 0.0}, {2.0, 1.0}}));
    check_invariants(two, h);
}

TEST_CASE("convex_hull: dos puntos iguales") {
    const Points two{{1.0, 1.0}, {1.0, 1.0}};
    CHECK_EQ(convex_hull(two), (Points{{1.0, 1.0}}));
}

TEST_CASE("convex_hull: tres puntos colineales") {
    CHECK_EQ(convex_hull({{0, 0}, {1, 1}, {2, 2}}), (Points{{0, 0}, {2, 2}}));
    CHECK_EQ(convex_hull({{2, 2}, {0, 0}, {1, 1}}), (Points{{0, 0}, {2, 2}}));
    CHECK_EQ(convex_hull({{0, 0}, {1, 0}, {2, 0}}), (Points{{0, 0}, {2, 0}}));      // horizontal
    CHECK_EQ(convex_hull({{5, 3}, {5, 1}, {5, 2}}), (Points{{5, 1}, {5, 3}}));      // vertical
    CHECK_EQ(convex_hull({{-1, 2}, {3, -6}, {1, -2}}), (Points{{-1, 2}, {3, -6}}));  // decreciente
}

TEST_CASE("convex_hull: muchos colineales con duplicados") {
    Points pts;
    for (int i = 0; i < 100; ++i) {
        pts.push_back({static_cast<double>(i), 2.0 * i});
        pts.push_back({static_cast<double>(i), 2.0 * i});
    }
    std::mt19937_64 rng(1);
    std::ranges::shuffle(pts, rng);
    const Points h = convex_hull(pts);
    CHECK_EQ(h, (Points{{0, 0}, {99, 198}}));
    check_invariants(pts, h);
}

TEST_CASE("convex_hull: todos los puntos iguales") {
    const Points pts(50, Point2D{1.5, 2.5});
    CHECK_EQ(convex_hull(pts), (Points{{1.5, 2.5}}));
}

TEST_CASE("convex_hull: triángulo en orden antihorario desde el menor lexicográfico") {
    CHECK_EQ(convex_hull({{2, 0}, {1, 1}, {0, 0}}), (Points{{0, 0}, {2, 0}, {1, 1}}));
    CHECK_EQ(convex_hull({{0, 0}, {1, 1}, {2, 0}}), (Points{{0, 0}, {2, 0}, {1, 1}}));
    CHECK_EQ(convex_hull({{1, 1}, {0, 0}, {2, 0}, {1, 1}}), (Points{{0, 0}, {2, 0}, {1, 1}}));
}

TEST_CASE("convex_hull: cuadrado con puntos interiores, sobre las aristas y duplicados") {
    const Points pts{{0, 0},     {1, 0},      {1, 1},   {0, 1},    // esquinas
                     {0.5, 0.5}, {0.25, 0.75}, {0.9, 0.1},        // interiores
                     {0.5, 0},   {1, 0.5},    {0, 0.5}, {0.5, 1},  // sobre las aristas
                     {0, 0},     {1, 1},      {0, 1}};             // duplicados
    const Points h = convex_hull(pts);
    CHECK_EQ(h, (Points{{0, 0}, {1, 0}, {1, 1}, {0, 1}}));
    check_invariants(pts, h);
}

TEST_CASE("convex_hull: puntos en posición convexa (parábola) salen todos y ordenados") {
    Points pts;
    for (int i = 0; i < 100; ++i) {
        const double x = static_cast<double>(i);
        pts.push_back({x, x * x});
    }
    const Points expected = pts;  // ya está en orden antihorario desde (0, 0)
    std::mt19937_64 rng(5);
    std::ranges::shuffle(pts, rng);
    const Points h = convex_hull(pts);
    CHECK_EQ(h, expected);
    check_invariants(pts, h);
}

TEST_CASE("convex_hull: circunferencia de 1000 puntos (todos son vértices)") {
    constexpr std::size_t n = 1000;
    Points pts;
    for (std::size_t i = 0; i < n; ++i) {
        const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n);
        pts.push_back({100.0 * std::cos(angle), 100.0 * std::sin(angle)});
    }
    std::mt19937_64 rng(11);
    std::ranges::shuffle(pts, rng);
    const Points h = convex_hull(pts);
    CHECK_EQ(h.size(), n);
    check_invariants(pts, h);
}

TEST_CASE("convex_hull: invariante ante permutaciones, inversión y duplicación de la entrada") {
    std::mt19937_64 rng(2025);
    std::uniform_int_distribution<std::size_t> size(3, 200);
    for (int trial = 0; trial < 50; ++trial) {
        Points pts = random_points(rng, size(rng), -10.0, 10.0);
        const Points h = convex_hull(pts);
        check_invariants(pts, h);

        std::ranges::shuffle(pts, rng);
        CHECK_EQ(convex_hull(pts), h);

        std::ranges::reverse(pts);
        CHECK_EQ(convex_hull(pts), h);

        Points doubled = pts;
        doubled.insert(doubled.end(), pts.begin(), pts.end());
        CHECK_EQ(convex_hull(doubled), h);
    }
}

TEST_CASE("convex_hull: coincide con la fuerza bruta en rejillas pequeñas") {
    // Rejilla 5x5 con hasta 10 puntos: abundan los duplicados y las colinealidades.
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int> size(0, 10);
    std::uniform_int_distribution<int> coord(0, 4);
    for (int trial = 0; trial < 3000; ++trial) {
        Points pts;
        const int n = size(rng);
        for (int i = 0; i < n; ++i) {
            pts.push_back({static_cast<double>(coord(rng)), static_cast<double>(coord(rng))});
        }
        const Points h = convex_hull(pts);
        const Points expected = brute_force_hull_vertices(pts);
        CHECK_MSG(same_set(h, expected), "entrada: " + geotest::to_string(pts) +
                                             "\n    obtenido: " + geotest::to_string(h) +
                                             "\n    esperado: " + geotest::to_string(expected));
        check_invariants(pts, h);
    }
}

TEST_CASE("convex_hull: propiedades en conjuntos aleatorios grandes") {
    std::mt19937_64 rng(777);

    const Points uniform = random_points(rng, 20000, -1000.0, 1000.0);
    check_invariants(uniform, convex_hull(uniform));

    std::normal_distribution<double> normal(0.0, 1.0);
    Points gaussian(20000);
    for (Point2D& p : gaussian) p = {normal(rng), normal(rng)};
    check_invariants(gaussian, convex_hull(gaussian));

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    Points disk(20000);
    for (Point2D& p : disk) {
        const double r = std::sqrt(unit(rng));
        const double angle = 2.0 * std::numbers::pi * unit(rng);
        p = {r * std::cos(angle), r * std::sin(angle)};
    }
    check_invariants(disk, convex_hull(disk));
}

TEST_CASE("convex_hull: robustez con puntos casi colineales de gran magnitud") {
    // Puntos sobre y = x/3 con x entre 1e15 y 2e15. La división redondea y en el último
    // bit (~0.06), así que los puntos quedan a distancias minúsculas de la recta mientras
    // los productos cruzados rondan 1e29, muy por encima de lo que double puede resolver.
    // Con un producto cruz ingenuo la "envolvente" sale con giros a la derecha; con el
    // predicado exacto es estrictamente convexa.
    Points pts;
    for (int i = 0; i < 1000; ++i) {
        const double x = 1e15 + 1e12 * static_cast<double>(i);
        pts.push_back({x, x / 3.0});
    }
    std::mt19937_64 rng(3);
    std::ranges::shuffle(pts, rng);
    const Points h = convex_hull(pts);
    CHECK(h.size() >= 2);
    check_invariants(pts, h);

    // Lo mismo con dos puntos lejanos por encima y por debajo de la recta.
    pts.push_back({1.5e15, 1e15});
    pts.push_back({1.5e15, -1e15});
    check_invariants(pts, convex_hull(pts));

    // Y con ternas aleatorias de doubles casi colineales (ruido de redondeo puro).
    std::uniform_real_distribution<double> coord(-1e9, 1e9);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int trial = 0; trial < 200; ++trial) {
        const Point2D a{coord(rng), coord(rng)};
        const Point2D b{coord(rng), coord(rng)};
        Points line;
        for (int i = 0; i < 50; ++i) {
            const double t = unit(rng);
            line.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
        }
        check_invariants(line, convex_hull(line));
    }
}

TEST_CASE("convex_hull: lanza std::invalid_argument con coordenadas no finitas") {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(convex_hull({{0, 0}, {nan, 1}, {1, 0}}), std::invalid_argument);
    CHECK_THROWS_AS(convex_hull({{0, 0}, {1, inf}}), std::invalid_argument);
    CHECK_THROWS_AS(convex_hull({{-inf, 0}}), std::invalid_argument);
}

TEST_CASE("convex_hull: acepta vector, array, initializer_list, list, span y vistas de ranges") {
    const Points expected{{0, 0}, {1, 0}, {1, 1}, {0, 1}};

    const Points vec{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}};
    CHECK_EQ(convex_hull(vec), expected);
    CHECK_EQ(convex_hull(Points{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}}), expected);  // rvalue

    const std::array<Point2D, 5> arr{{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}}};
    CHECK_EQ(convex_hull(arr), expected);

    CHECK_EQ(convex_hull({{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}}), expected);

    const std::list<Point2D> lst{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}};
    CHECK_EQ(convex_hull(lst), expected);

    CHECK_EQ(convex_hull(std::span<const Point2D>{vec}), expected);
    CHECK_EQ(convex_hull(std::span<const Point2D>{vec}.first(4)), expected);

    const Points with_outliers{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}, {-5, 0}, {7, 7}};
    auto inside_unit_square = with_outliers | std::views::filter([](Point2D p) {
                                  return p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1;
                              });
    CHECK_EQ(convex_hull(inside_unit_square), expected);
}

GEOTEST_MAIN

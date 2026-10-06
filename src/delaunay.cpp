// Triangulación de Delaunay por divide y vencerás.
//
// Algoritmo: L. J. Guibas y J. Stolfi, "Primitives for the Manipulation of General
// Subdivisions and the Computation of Voronoi Diagrams", ACM Transactions on Graphics
// 4(2):74-123, 1985. Implementación propia, sobre la estructura quad-edge del mismo
// artículo.
//
// 1. Se ordenan los puntos lexicográficamente y se descartan los repetidos (queda el de
//    menor índice).
// 2. Se triangula recursivamente cada mitad; los casos base son 2 y 3 puntos.
// 3. Para fusionar ambas mitades se busca su tangente común inferior y desde ella se
//    "teje" hacia arriba: en cada paso se elige entre el candidato de la izquierda y el
//    de la derecha con el test del círculo, borrando antes las aristas de cada mitad que
//    dejan de ser de Delaunay.
// 4. Cada cara acotada de la subdivisión resultante es un triángulo; recorriéndolas se
//    arman los triángulos y sus adyacencias.
//
// Todas las decisiones pasan por orientation() e incircle(), que son exactos, de modo que
// el resultado es combinatoriamente correcto para los puntos tal cual se representan en
// double, incluidos los colineales y los concíclicos.

#include "geometra/delaunay.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "geometra/predicates.hpp"

namespace geo {
namespace {

/// Arista dirigida de la estructura quad-edge: 4 * (número de arista) + rotación. Las
/// rotaciones 0 y 2 son los dos sentidos de la arista primal; 1 y 3, los de su dual.
using EdgeRef = std::size_t;

constexpr EdgeRef rot(EdgeRef e) noexcept { return (e & ~EdgeRef{3}) | ((e + 1) & 3); }
constexpr EdgeRef sym(EdgeRef e) noexcept { return (e & ~EdgeRef{3}) | ((e + 2) & 3); }
constexpr EdgeRef rot_inv(EdgeRef e) noexcept { return (e & ~EdgeRef{3}) | ((e + 3) & 3); }

bool ccw(Point2D a, Point2D b, Point2D c) noexcept {
    return orientation(a, b, c) == Orientation::counterclockwise;
}

/// `d` está estrictamente dentro de la circunferencia por a, b, c (en sentido antihorario).
bool in_circle(Point2D a, Point2D b, Point2D c, Point2D d) noexcept {
    return incircle(a, b, c, d) > 0.0;
}

/// Subdivisión del plano en quad-edges cuyos vértices son índices a un conjunto de
/// puntos fijo. Las aristas borradas se reutilizan.
class Subdivision {
public:
    explicit Subdivision(std::span<const Point2D> points) : points_(points) {}

    void reserve(std::size_t edges) {
        next_.reserve(4 * edges);
        org_.reserve(4 * edges);
        alive_.reserve(edges);
    }

    [[nodiscard]] EdgeRef onext(EdgeRef e) const noexcept { return next_[e]; }
    [[nodiscard]] EdgeRef oprev(EdgeRef e) const noexcept { return rot(onext(rot(e))); }
    [[nodiscard]] EdgeRef lnext(EdgeRef e) const noexcept { return rot(onext(rot_inv(e))); }
    [[nodiscard]] EdgeRef rprev(EdgeRef e) const noexcept { return onext(sym(e)); }

    [[nodiscard]] std::size_t org(EdgeRef e) const noexcept { return org_[e]; }
    [[nodiscard]] std::size_t dest(EdgeRef e) const noexcept { return org_[sym(e)]; }
    [[nodiscard]] Point2D point(std::size_t v) const noexcept { return points_[v]; }
    [[nodiscard]] Point2D org_point(EdgeRef e) const noexcept { return points_[org(e)]; }
    [[nodiscard]] Point2D dest_point(EdgeRef e) const noexcept { return points_[dest(e)]; }

    /// `p` queda estrictamente a la derecha / izquierda de la recta dirigida de `e`.
    [[nodiscard]] bool right_of(Point2D p, EdgeRef e) const noexcept {
        return ccw(p, dest_point(e), org_point(e));
    }
    [[nodiscard]] bool left_of(Point2D p, EdgeRef e) const noexcept {
        return ccw(p, org_point(e), dest_point(e));
    }

    /// Cantidad de aristas dirigidas reservadas (vivas o no): cota de los EdgeRef válidos.
    [[nodiscard]] std::size_t slots() const noexcept { return next_.size(); }
    [[nodiscard]] bool alive(EdgeRef e) const noexcept { return alive_[e / 4] != 0; }

    /// Arista aislada de `from` a `to`.
    EdgeRef make_edge(std::size_t from, std::size_t to) {
        EdgeRef e = 0;
        if (free_.empty()) {
            e = next_.size();
            next_.resize(e + 4);
            org_.resize(e + 4);
            alive_.push_back(0);
        } else {
            e = free_.back();
            free_.pop_back();
        }
        next_[e] = e;
        next_[e + 1] = e + 3;
        next_[e + 2] = e + 2;
        next_[e + 3] = e + 1;
        org_[e] = from;
        org_[e + 2] = to;
        alive_[e / 4] = 1;
        return e;
    }

    /// Operador Splice de Guibas y Stolfi: une o separa los anillos de `a` y `b`.
    void splice(EdgeRef a, EdgeRef b) noexcept {
        const EdgeRef alpha = rot(onext(a));
        const EdgeRef beta = rot(onext(b));
        std::swap(next_[a], next_[b]);
        std::swap(next_[alpha], next_[beta]);
    }

    /// Nueva arista del destino de `a` al origen de `b`, con las mismas caras a la
    /// izquierda que `a` y `b`.
    EdgeRef connect(EdgeRef a, EdgeRef b) {
        const EdgeRef e = make_edge(dest(a), org(b));
        splice(e, lnext(a));
        splice(sym(e), b);
        return e;
    }

    void remove(EdgeRef e) {
        splice(e, oprev(e));
        splice(sym(e), oprev(sym(e)));
        alive_[e / 4] = 0;
        free_.push_back(e & ~EdgeRef{3});
    }

private:
    std::span<const Point2D> points_;
    std::vector<EdgeRef> next_;      // Onext de cada arista dirigida
    std::vector<std::size_t> org_;   // vértice de origen (solo rotaciones 0 y 2)
    std::vector<unsigned char> alive_;
    std::vector<EdgeRef> free_;
};

/// Aristas de la envolvente convexa de una triangulación parcial: `left` sale del vértice
/// de más a la izquierda en sentido antihorario y `right` sale del de más a la derecha en
/// sentido horario.
struct HullEdges {
    EdgeRef left;
    EdgeRef right;
};

class DelaunayBuilder {
public:
    DelaunayBuilder(Subdivision& mesh, std::span<const std::size_t> sorted_ids)
        : mesh_(mesh), ids_(sorted_ids) {}

    /// Triangula ids_[lo, hi), con al menos 2 puntos.
    HullEdges build(std::size_t lo, std::size_t hi) {
        Subdivision& m = mesh_;
        const std::size_t n = hi - lo;

        if (n == 2) {
            const EdgeRef a = m.make_edge(ids_[lo], ids_[lo + 1]);
            return {a, sym(a)};
        }
        if (n == 3) {
            const EdgeRef a = m.make_edge(ids_[lo], ids_[lo + 1]);
            const EdgeRef b = m.make_edge(ids_[lo + 1], ids_[lo + 2]);
            m.splice(sym(a), b);
            switch (orientation(m.point(ids_[lo]), m.point(ids_[lo + 1]), m.point(ids_[lo + 2]))) {
                case Orientation::counterclockwise:
                    m.connect(b, a);
                    return {a, sym(b)};
                case Orientation::clockwise: {
                    const EdgeRef c = m.connect(b, a);
                    return {sym(c), c};
                }
                case Orientation::collinear: break;
            }
            return {a, sym(b)};
        }

        const std::size_t mid = lo + n / 2;
        auto [ldo, ldi] = build(lo, mid);
        auto [rdi, rdo] = build(mid, hi);

        // Tangente común inferior de ambas mitades.
        while (true) {
            if (m.left_of(m.org_point(rdi), ldi)) {
                ldi = m.lnext(ldi);
            } else if (m.right_of(m.org_point(ldi), rdi)) {
                rdi = m.rprev(rdi);
            } else {
                break;
            }
        }

        EdgeRef basel = m.connect(sym(rdi), ldi);
        if (m.org(ldi) == m.org(ldo)) ldo = sym(basel);
        if (m.org(rdi) == m.org(rdo)) rdo = basel;

        // Un candidato es válido si queda por encima de la base actual.
        const auto valid = [&](EdgeRef e) { return m.right_of(m.dest_point(e), basel); };

        while (true) {
            EdgeRef lcand = m.onext(sym(basel));
            if (valid(lcand)) {
                while (in_circle(m.dest_point(basel), m.org_point(basel), m.dest_point(lcand),
                                 m.dest_point(m.onext(lcand)))) {
                    const EdgeRef t = m.onext(lcand);
                    m.remove(lcand);
                    lcand = t;
                }
            }

            EdgeRef rcand = m.oprev(basel);
            if (valid(rcand)) {
                while (in_circle(m.dest_point(basel), m.org_point(basel), m.dest_point(rcand),
                                 m.dest_point(m.oprev(rcand)))) {
                    const EdgeRef t = m.oprev(rcand);
                    m.remove(rcand);
                    rcand = t;
                }
            }

            const bool lvalid = valid(lcand);
            const bool rvalid = valid(rcand);
            if (!lvalid && !rvalid) break;  // la base es la tangente común superior

            if (!lvalid || (rvalid && in_circle(m.dest_point(lcand), m.org_point(lcand),
                                                m.org_point(rcand), m.dest_point(rcand)))) {
                basel = m.connect(rcand, sym(basel));
            } else {
                basel = m.connect(sym(basel), sym(lcand));
            }
        }
        return {ldo, rdo};
    }

private:
    Subdivision& mesh_;
    std::span<const std::size_t> ids_;
};

/// Recorre las caras acotadas de la subdivisión (todas triángulos) y arma la salida en
/// forma canónica: cada triángulo empieza por su menor vértice y la lista se ordena.
std::vector<Triangle> extract_triangles(const Subdivision& m) {
    constexpr std::size_t none = no_neighbor;
    std::vector<std::size_t> face(m.slots(), none);  // triángulo a la izquierda de cada arista
    std::vector<std::array<EdgeRef, 3>> edges;

    for (EdgeRef q = 0; q < m.slots(); q += 4) {
        if (!m.alive(q)) continue;
        for (const EdgeRef e0 : {q, sym(q)}) {
            if (face[e0] != none) continue;
            const EdgeRef e1 = m.lnext(e0);
            const EdgeRef e2 = m.lnext(e1);
            // La cara exterior también puede tener tres aristas, pero en sentido horario.
            if (m.lnext(e2) != e0 || !ccw(m.org_point(e0), m.org_point(e1), m.org_point(e2))) {
                continue;
            }
            face[e0] = face[e1] = face[e2] = edges.size();
            edges.push_back({e0, e1, e2});
        }
    }

    // La arista ei va de v[i] a v[i+1], así que es la opuesta a v[i+2].
    std::vector<Triangle> raw(edges.size());
    for (std::size_t t = 0; t < edges.size(); ++t) {
        const auto [e0, e1, e2] = edges[t];
        Triangle tri{{m.org(e0), m.org(e1), m.org(e2)},
                     {face[sym(e1)], face[sym(e2)], face[sym(e0)]}};
        const auto first = std::ranges::min_element(tri.v) - tri.v.begin();
        std::ranges::rotate(tri.v, tri.v.begin() + first);
        std::ranges::rotate(tri.neighbor, tri.neighbor.begin() + first);
        raw[t] = tri;
    }

    std::vector<std::size_t> order(raw.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return raw[a].v < raw[b].v; });

    std::vector<std::size_t> new_id(raw.size());
    for (std::size_t pos = 0; pos < order.size(); ++pos) new_id[order[pos]] = pos;

    std::vector<Triangle> triangles;
    triangles.reserve(raw.size());
    for (const std::size_t old : order) {
        Triangle tri = raw[old];
        for (std::size_t& nb : tri.neighbor) {
            if (nb != none) nb = new_id[nb];
        }
        triangles.push_back(tri);
    }
    return triangles;
}

Triangulation divide_and_conquer(std::vector<Point2D> pts) {
    for (const Point2D& p : pts) {
        if (!is_finite(p)) {
            throw std::invalid_argument(
                "geo::delaunay_triangulation: las coordenadas deben ser finitas");
        }
    }

    Triangulation result;
    result.points = std::move(pts);
    const std::span<const Point2D> points{result.points};

    // Índices en orden lexicográfico de sus puntos; entre repetidos, el menor índice primero.
    std::vector<std::size_t> ids(points.size());
    std::iota(ids.begin(), ids.end(), std::size_t{0});
    std::ranges::sort(ids, [&](std::size_t a, std::size_t b) {
        if (points[a] != points[b]) return points[a] < points[b];
        return a < b;
    });
    const auto repeated = std::ranges::unique(
        ids, [&](std::size_t a, std::size_t b) { return points[a] == points[b]; });
    ids.erase(repeated.begin(), repeated.end());

    if (ids.size() < 3) return result;

    Subdivision mesh(points);
    mesh.reserve(3 * ids.size());
    DelaunayBuilder(mesh, ids).build(0, ids.size());
    result.triangles = extract_triangles(mesh);
    return result;
}

}  // namespace

Triangulation delaunay_triangulation(std::span<const Point2D> points) {
    return divide_and_conquer(std::vector<Point2D>(points.begin(), points.end()));
}

Triangulation delaunay_triangulation(std::initializer_list<Point2D> points) {
    return divide_and_conquer(std::vector<Point2D>(points));
}

}  // namespace geo

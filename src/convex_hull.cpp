// Envolvente convexa por cadena monótona.
//
// Algoritmo: A. M. Andrew, "Another efficient algorithm for convex hulls in two
// dimensions", Information Processing Letters 9(5):216-219, 1979. Implementación propia.
//
// 1. Se ordenan los puntos lexicográficamente (x, luego y) y se eliminan los duplicados.
// 2. Recorriendo de izquierda a derecha se construye la cadena inferior: cada nuevo punto
//    expulsa del final de la cadena a los que dejan de formar un giro a la izquierda.
// 3. Recorriendo de derecha a izquierda se construye la cadena superior del mismo modo.
// 4. Ambas cadenas comparten sus extremos; concatenadas sin repetirlos forman la
//    envolvente en sentido antihorario.
//
// Al expulsar también cuando el giro es nulo (colineal), solo sobreviven los vértices
// extremos. Todas las decisiones pasan por `orientation()`, que es exacto, de modo que el
// resultado es combinatoriamente correcto para los puntos tal cual se representan en
// double, por degenerada que sea la entrada.

#include "geometra/convex_hull.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

#include "geometra/predicates.hpp"

namespace geo {
namespace {

/// Añade `p` al final de una cadena convexa, retirando antes los vértices que con `p`
/// dejarían de girar a la izquierda. El primer vértice de la cadena nunca se retira.
void extend_chain(std::vector<Point2D>& chain, Point2D p) {
    while (chain.size() >= 2 &&
           orientation(chain[chain.size() - 2], chain.back(), p) != Orientation::counterclockwise) {
        chain.pop_back();
    }
    chain.push_back(p);
}

std::vector<Point2D> monotone_chain(std::vector<Point2D> pts) {
    for (const Point2D& p : pts) {
        if (!is_finite(p)) {
            throw std::invalid_argument("geo::convex_hull: las coordenadas deben ser finitas");
        }
    }

    std::ranges::sort(pts);
    pts.erase(std::ranges::unique(pts).begin(), pts.end());

    const std::size_t n = pts.size();
    if (n <= 2) return pts;  // 0, 1 o 2 puntos distintos ya son su propia envolvente

    std::vector<Point2D> lower;
    std::vector<Point2D> upper;
    lower.reserve(n);
    upper.reserve(n);

    for (const Point2D& p : pts) extend_chain(lower, p);                      // izquierda → derecha
    for (auto it = pts.rbegin(); it != pts.rend(); ++it) extend_chain(upper, *it);  // derecha → izquierda

    // `lower` termina en el punto donde empieza `upper`, y `upper` termina en el punto
    // donde empieza `lower`: se omite el último de cada una al concatenar.
    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());
    return lower;
}

}  // namespace

std::vector<Point2D> convex_hull(std::span<const Point2D> points) {
    return monotone_chain(std::vector<Point2D>(points.begin(), points.end()));
}

std::vector<Point2D> convex_hull(std::initializer_list<Point2D> points) {
    return monotone_chain(std::vector<Point2D>(points));
}

}  // namespace geo

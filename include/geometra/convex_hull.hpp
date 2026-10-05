#pragma once

#include <concepts>
#include <initializer_list>
#include <ranges>
#include <span>
#include <vector>

#include "geometra/point.hpp"

namespace geo {

/// Envolvente convexa (cierre convexo) de un conjunto de puntos del plano.
///
/// Contrato de salida:
///  - Los vértices se devuelven en sentido antihorario, empezando por el punto
///    lexicográficamente menor (menor `x`; a igual `x`, menor `y`).
///  - Solo se devuelven los vértices extremos: sin duplicados y sin puntos colineales
///    sobre las aristas.
///  - Casos degenerados: 0 puntos → vacío; 1 punto → ese punto; 2 puntos distintos →
///    ambos; todos los puntos iguales → un punto; todos colineales → los dos extremos.
///  - Las coordenadas deben ser finitas; si no, lanza `std::invalid_argument`.
///
/// Algoritmo: cadena monótona de Andrew con predicado de orientación exacto, O(n log n).
/// La entrada no se modifica.
[[nodiscard]] std::vector<Point2D> convex_hull(std::span<const Point2D> points);

/// Sobrecarga para listas literales: `geo::convex_hull({{0, 0}, {1, 0}, {0, 1}})`.
[[nodiscard]] std::vector<Point2D> convex_hull(std::initializer_list<Point2D> points);

/// Sobrecarga para cualquier rango de `Point2D` que no sea contiguo en memoria
/// (`std::list`, vistas de `std::ranges`, ...). Los contenedores contiguos (`std::vector`,
/// `std::array`, arreglos C) usan directamente la versión con `std::span`.
template <std::ranges::input_range R>
    requires(!std::convertible_to<R, std::span<const Point2D>>) &&
            std::convertible_to<std::ranges::range_reference_t<R>, Point2D>
[[nodiscard]] std::vector<Point2D> convex_hull(R&& points) {
    std::vector<Point2D> copy;
    if constexpr (std::ranges::sized_range<R>) {
        copy.reserve(static_cast<std::size_t>(std::ranges::size(points)));
    }
    for (auto&& p : points) {
        copy.push_back(static_cast<Point2D>(p));
    }
    return convex_hull(std::span<const Point2D>{copy});
}

}  // namespace geo

#pragma once

#include <concepts>
#include <initializer_list>
#include <ranges>
#include <span>
#include <vector>

#include "geometra/point.hpp"
#include "geometra/triangulation.hpp"

namespace geo {

/// Triangulación de Delaunay de un conjunto de puntos del plano.
///
/// Contrato de salida:
///  - `points` es una copia de la entrada, en el mismo orden: los índices de los
///    triángulos se refieren a las posiciones de la entrada original.
///  - Triángulos en sentido antihorario con sus adyacencias, en la forma que describe
///    `Triangulation`. Se ordenan lexicográficamente por `v`, así que la salida es
///    determinista.
///  - Propiedad de Delaunay: ningún punto queda estrictamente dentro de la circunferencia
///    circunscrita de un triángulo, según `incircle()`, que es exacto.
///  - La unión de los triángulos es la envolvente convexa, y todo punto distinto es
///    vértice de algún triángulo, incluidos los colineales sobre el borde.
///  - Si hay puntos repetidos, solo se usa el de menor índice; los demás no aparecen.
///  - Con cuatro o más puntos concíclicos la triangulación de Delaunay no es única: se
///    devuelve una de ellas, que depende solo del conjunto de puntos y no de su orden en
///    la entrada.
///  - Casos degenerados: menos de 3 puntos distintos, o todos colineales → sin triángulos.
///  - Las coordenadas deben ser finitas; si no, lanza `std::invalid_argument`. Valen las
///    precondiciones de magnitud de `incircle()`.
///
/// Algoritmo: divide y vencerás de Guibas y Stolfi con predicados exactos, O(n log n).
/// La entrada no se modifica.
[[nodiscard]] Triangulation delaunay_triangulation(std::span<const Point2D> points);

/// Sobrecarga para listas literales:
/// `geo::delaunay_triangulation({{0, 0}, {1, 0}, {0, 1}})`.
[[nodiscard]] Triangulation delaunay_triangulation(std::initializer_list<Point2D> points);

/// Sobrecarga para cualquier rango de `Point2D` que no sea contiguo en memoria
/// (`std::list`, vistas de `std::ranges`, ...). Los contenedores contiguos usan
/// directamente la versión con `std::span`.
template <std::ranges::input_range R>
    requires(!std::convertible_to<R, std::span<const Point2D>>) &&
            std::convertible_to<std::ranges::range_reference_t<R>, Point2D>
[[nodiscard]] Triangulation delaunay_triangulation(R&& points) {
    std::vector<Point2D> copy;
    if constexpr (std::ranges::sized_range<R>) {
        copy.reserve(static_cast<std::size_t>(std::ranges::size(points)));
    }
    for (auto&& p : points) {
        copy.push_back(static_cast<Point2D>(p));
    }
    return delaunay_triangulation(std::span<const Point2D>{copy});
}

}  // namespace geo

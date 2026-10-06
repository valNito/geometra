#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

#include "geometra/point.hpp"

namespace geo {

/// Valor de `Triangle::neighbor` para una arista sin triángulo vecino (en el borde).
inline constexpr std::size_t no_neighbor = std::numeric_limits<std::size_t>::max();

/// Triángulo de una `Triangulation`, dado por índices.
struct Triangle {
    /// Vértices: índices a `Triangulation::points`, en sentido antihorario. `v[0]` es el
    /// menor de los tres.
    std::array<std::size_t, 3> v{};

    /// `neighbor[i]` es el índice del triángulo que comparte la arista opuesta a `v[i]`
    /// (la que une `v[(i + 1) % 3]` con `v[(i + 2) % 3]`), o `no_neighbor` si esa arista
    /// está en el borde de la triangulación.
    std::array<std::size_t, 3> neighbor{};

    friend constexpr bool operator==(const Triangle&, const Triangle&) = default;
};

/// Malla de triángulos con adyacencias.
///
/// Es un tipo de valor autocontenido: guarda los puntos junto con los triángulos, de modo
/// que los índices de `triangles` siempre se refieren a `points`.
struct Triangulation {
    std::vector<Point2D> points;
    std::vector<Triangle> triangles;
};

}  // namespace geo

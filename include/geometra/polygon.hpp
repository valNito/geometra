#pragma once

#include <cmath>
#include <span>

#include "geometra/point.hpp"

namespace geo {

/// Área con signo de un polígono simple dado por sus vértices en orden (fórmula del
/// cordón de zapato). Positiva si los vértices van en sentido antihorario, negativa en
/// sentido horario. Con menos de 3 vértices devuelve 0.
[[nodiscard]] double signed_area(std::span<const Point2D> polygon) noexcept;

/// Área (sin signo) de un polígono simple.
[[nodiscard]] inline double area(std::span<const Point2D> polygon) noexcept {
    return std::fabs(signed_area(polygon));
}

}  // namespace geo

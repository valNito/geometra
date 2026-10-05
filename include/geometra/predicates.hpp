#pragma once

#include <ostream>

#include "geometra/point.hpp"

namespace geo {

/// Giro de la terna ordenada (a, b, c).
enum class Orientation : int {
    clockwise = -1,        ///< c queda a la derecha de la recta dirigida a→b
    collinear = 0,         ///< a, b y c están alineados (incluye puntos repetidos)
    counterclockwise = 1,  ///< c queda a la izquierda de la recta dirigida a→b
};

[[nodiscard]] constexpr Orientation opposite(Orientation o) noexcept {
    return static_cast<Orientation>(-static_cast<int>(o));
}

[[nodiscard]] constexpr const char* to_string(Orientation o) noexcept {
    switch (o) {
        case Orientation::clockwise: return "clockwise";
        case Orientation::collinear: return "collinear";
        case Orientation::counterclockwise: return "counterclockwise";
    }
    return "?";
}

inline std::ostream& operator<<(std::ostream& os, Orientation o) { return os << to_string(o); }

/// Determinante de orientación
///
///     | a.x - c.x   a.y - c.y |
///     | b.x - c.x   b.y - c.y |
///
/// es decir, el doble del área con signo del triángulo (a, b, c).
///
/// El signo del valor devuelto es siempre el correcto: positivo, negativo o exactamente
/// cero según la posición real de los puntos, con independencia del redondeo de coma
/// flotante (aritmética adaptativa de Shewchuk). La magnitud es aproximada salvo en los
/// casos casi degenerados, en los que se calcula exactamente.
///
/// Precondición: coordenadas finitas y sin desbordamiento en los productos.
[[nodiscard]] double orient2d(Point2D a, Point2D b, Point2D c) noexcept;

/// Orientación exacta de la terna (a, b, c). Véase `orient2d()`.
[[nodiscard]] inline Orientation orientation(Point2D a, Point2D b, Point2D c) noexcept {
    const double det = orient2d(a, b, c);
    if (det > 0.0) return Orientation::counterclockwise;
    if (det < 0.0) return Orientation::clockwise;
    return Orientation::collinear;
}

}  // namespace geo

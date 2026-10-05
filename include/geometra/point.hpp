#pragma once

#include <cmath>
#include <compare>
#include <ostream>

namespace geo {

/// Punto del plano con coordenadas `double`. También se usa como vector de
/// desplazamiento (la diferencia entre dos puntos).
///
/// Es un agregado: `Point2D{1.0, 2.0}` y `Point2D{.x = 1.0, .y = 2.0}` son válidos.
/// El orden que define `<=>` es lexicográfico (primero `x`; a igual `x`, `y`). Es el
/// orden canónico que usan los algoritmos de la librería.
struct Point2D {
    double x = 0.0;
    double y = 0.0;

    friend constexpr bool operator==(const Point2D&, const Point2D&) = default;
    friend constexpr auto operator<=>(const Point2D&, const Point2D&) = default;

    constexpr Point2D& operator+=(Point2D o) noexcept {
        x += o.x;
        y += o.y;
        return *this;
    }
    constexpr Point2D& operator-=(Point2D o) noexcept {
        x -= o.x;
        y -= o.y;
        return *this;
    }
    constexpr Point2D& operator*=(double s) noexcept {
        x *= s;
        y *= s;
        return *this;
    }
    constexpr Point2D& operator/=(double s) noexcept {
        x /= s;
        y /= s;
        return *this;
    }
};

[[nodiscard]] constexpr Point2D operator+(Point2D a, Point2D b) noexcept {
    return {a.x + b.x, a.y + b.y};
}
[[nodiscard]] constexpr Point2D operator-(Point2D a, Point2D b) noexcept {
    return {a.x - b.x, a.y - b.y};
}
[[nodiscard]] constexpr Point2D operator-(Point2D a) noexcept { return {-a.x, -a.y}; }
[[nodiscard]] constexpr Point2D operator*(Point2D a, double s) noexcept {
    return {a.x * s, a.y * s};
}
[[nodiscard]] constexpr Point2D operator*(double s, Point2D a) noexcept {
    return {s * a.x, s * a.y};
}
[[nodiscard]] constexpr Point2D operator/(Point2D a, double s) noexcept {
    return {a.x / s, a.y / s};
}

/// Producto punto.
[[nodiscard]] constexpr double dot(Point2D a, Point2D b) noexcept {
    return a.x * b.x + a.y * b.y;
}

/// Producto cruz 2D (componente z de a × b).
///
/// Su signo indica el giro de `a` hacia `b`, pero evaluado en coma flotante puede ser
/// incorrecto en configuraciones casi degeneradas. Para tomar decisiones combinatorias
/// (¿gira a la izquierda o a la derecha?) use `orientation()` de predicates.hpp, que es
/// exacto; `cross()` sirve para magnitudes (áreas, momentos).
[[nodiscard]] constexpr double cross(Point2D a, Point2D b) noexcept {
    return a.x * b.y - a.y * b.x;
}

[[nodiscard]] constexpr double squared_norm(Point2D v) noexcept { return dot(v, v); }
[[nodiscard]] inline double norm(Point2D v) noexcept { return std::hypot(v.x, v.y); }
[[nodiscard]] constexpr double squared_distance(Point2D a, Point2D b) noexcept {
    return squared_norm(a - b);
}
[[nodiscard]] inline double distance(Point2D a, Point2D b) noexcept { return norm(a - b); }
[[nodiscard]] constexpr Point2D midpoint(Point2D a, Point2D b) noexcept {
    return {0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
}

/// `true` si ambas coordenadas son finitas (ni NaN ni infinito).
[[nodiscard]] inline bool is_finite(Point2D p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y);
}

inline std::ostream& operator<<(std::ostream& os, Point2D p) {
    return os << '(' << p.x << ", " << p.y << ')';
}

}  // namespace geo

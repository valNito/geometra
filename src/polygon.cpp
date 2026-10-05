#include "geometra/polygon.hpp"

#include <cstddef>

namespace geo {

double signed_area(std::span<const Point2D> polygon) noexcept {
    const std::size_t n = polygon.size();
    if (n < 3) return 0.0;

    // Se suma relativo al primer vértice: reduce la cancelación catastrófica cuando las
    // coordenadas son grandes respecto al tamaño del polígono.
    const Point2D origin = polygon[0];
    double twice_area = 0.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        twice_area += cross(polygon[i] - origin, polygon[i + 1] - origin);
    }
    return 0.5 * twice_area;
}

}  // namespace geo

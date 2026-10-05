// Ejemplo mínimo de uso de Geometra: envolvente convexa de un conjunto de puntos.
//
// Desde la raíz del proyecto:
//   cmake -S . -B build && cmake --build build && ./build/examples/basic_usage

#include <geometra/geometra.hpp>

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

int main() {
    const std::vector<geo::Point2D> points = {
        {0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0},  // esquinas de un cuadrado
        {2.0, 2.0}, {1.0, 3.0}, {3.0, 1.0},              // puntos interiores
        {2.0, 0.0}, {4.0, 2.0},                          // sobre las aristas (no son vértices)
        {0.0, 0.0}, {4.0, 4.0},                          // duplicados
    };

    // Toda la complejidad queda detrás de esta llamada.
    const std::vector<geo::Point2D> hull = geo::convex_hull(points);

    std::cout << "Geometra " << geo::version_string << "\n\n";
    std::cout << "Entrada: " << points.size() << " puntos\n";
    std::cout << "Envolvente convexa: " << hull.size()
              << " vértices, en sentido antihorario desde el menor lexicográfico:\n";
    for (const geo::Point2D& p : hull) {
        std::cout << "  " << p << '\n';
    }
    std::cout << "Área: " << geo::area(hull) << "\n\n";

    // Con el predicado de orientación (exacto) se clasifica cualquier punto respecto a la
    // envolvente: está dentro si queda a la izquierda de todas las aristas.
    const auto classify = [&hull](geo::Point2D p) -> std::string {
        bool on_boundary = false;
        for (std::size_t i = 0; i < hull.size(); ++i) {
            const geo::Point2D a = hull[i];
            const geo::Point2D b = hull[(i + 1) % hull.size()];
            switch (geo::orientation(a, b, p)) {
                case geo::Orientation::clockwise: return "fuera";
                case geo::Orientation::collinear: on_boundary = true; break;
                case geo::Orientation::counterclockwise: break;
            }
        }
        return on_boundary ? "en el borde" : "dentro";
    };

    for (const geo::Point2D p : {geo::Point2D{1.0, 1.0}, geo::Point2D{4.0, 2.0},
                                 geo::Point2D{5.0, 1.0}, geo::Point2D{4.0, 4.0}}) {
        std::cout << p << " está " << classify(p) << '\n';
    }
    return 0;
}

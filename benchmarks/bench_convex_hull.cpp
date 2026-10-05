// Benchmark básico de geo::convex_hull.
//
// Mide la mediana de varias repeticiones para distintos tamaños y distribuciones de
// puntos, y la compara con el costo de solo ordenar y desduplicar los mismos puntos
// (la cota inferior natural del algoritmo de cadena monótona).
//
// Uso: bench_convex_hull [n_maximo]      (por defecto 1 000 000; compilar en Release)

#include <geometra/geometra.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <random>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using Points = std::vector<geo::Point2D>;
using Generator = Points (*)(std::size_t, std::mt19937_64&);

Points uniform_square(std::size_t n, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Points pts(n);
    for (auto& p : pts) p = {u(rng), u(rng)};
    return pts;
}

Points uniform_disk(std::size_t n, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Points pts(n);
    for (auto& p : pts) {
        const double r = std::sqrt(u(rng));
        const double a = 2.0 * std::numbers::pi * u(rng);
        p = {r * std::cos(a), r * std::sin(a)};
    }
    return pts;
}

// Peor caso para el tamaño de la salida: todos los puntos son vértices.
Points on_circle(std::size_t n, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> u(0.0, 2.0 * std::numbers::pi);
    Points pts(n);
    for (auto& p : pts) {
        const double a = u(rng);
        p = {std::cos(a), std::sin(a)};
    }
    return pts;
}

Points gaussian(std::size_t n, std::mt19937_64& rng) {
    std::normal_distribution<double> g(0.0, 1.0);
    Points pts(n);
    for (auto& p : pts) p = {g(rng), g(rng)};
    return pts;
}

// Coordenadas enteras en una rejilla 1000x1000: muchos duplicados y colinealidades.
Points integer_grid(std::size_t n, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> u(0, 999);
    Points pts(n);
    for (auto& p : pts) p = {static_cast<double>(u(rng)), static_cast<double>(u(rng))};
    return pts;
}

std::size_t sink = 0;  // evita que el optimizador elimine el trabajo

template <class F>
double median_ms(F&& f, int reps) {
    std::vector<double> times;
    times.reserve(static_cast<std::size_t>(reps));
    for (int i = 0; i < reps; ++i) {
        const auto t0 = Clock::now();
        f();
        const auto t1 = Clock::now();
        times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    std::sort(times.begin(), times.end());
    return times[times.size() / 2];
}

int repetitions(std::size_t n) {
    if (n <= 1'000) return 50;
    if (n <= 10'000) return 20;
    if (n <= 100'000) return 7;
    return 3;
}

}  // namespace

int main(int argc, char** argv) {
    const std::size_t max_n = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 1'000'000;
    const std::size_t sizes[] = {1'000, 10'000, 100'000, 1'000'000};

    struct Distribution {
        const char* name;
        Generator generate;
    };
    const Distribution distributions[] = {
        {"uniforme en cuadrado", &uniform_square},
        {"uniforme en disco", &uniform_disk},
        {"circunferencia", &on_circle},
        {"gaussiana", &gaussian},
        {"rejilla entera 1000x1000", &integer_grid},
    };

    std::printf("Geometra %s - benchmark de convex_hull (mediana de R repeticiones)\n\n",
                geo::version_string);
    std::printf("%-26s %10s %9s %4s %12s %10s %12s %8s\n", "distribucion", "n", "vertices",
                "R", "hull (ms)", "ns/punto", "sort+uniq", "hull/sort");
    std::printf("%-26s %10s %9s %4s %12s %10s %12s %8s\n", "", "", "", "", "", "", "(ms)", "");

    for (const Distribution& dist : distributions) {
        for (const std::size_t n : sizes) {
            if (n > max_n) continue;
            std::mt19937_64 rng(12345);
            const Points pts = dist.generate(n, rng);
            const int reps = repetitions(n);

            const double hull_ms = median_ms(
                [&] {
                    const Points h = geo::convex_hull(pts);
                    sink += h.size();
                },
                reps);

            const double sort_ms = median_ms(
                [&] {
                    Points copy(pts);
                    std::sort(copy.begin(), copy.end());
                    copy.erase(std::unique(copy.begin(), copy.end()), copy.end());
                    sink += copy.size();
                },
                reps);

            const std::size_t vertices = geo::convex_hull(pts).size();
            std::printf("%-26s %10zu %9zu %4d %12.3f %10.1f %12.3f %8.2f\n", dist.name, n,
                        vertices, reps, hull_ms, hull_ms * 1e6 / static_cast<double>(n), sort_ms,
                        hull_ms / sort_ms);
        }
        std::printf("\n");
    }

    std::printf("(control: %zu)\n", sink);
    return 0;
}

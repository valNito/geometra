#include "geometra/predicates.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "exact_reference.hpp"
#include "geotest.hpp"

namespace {

using geo::Orientation;
using geo::orient2d;
using geo::orientation;
using geo::Point2D;

constexpr Orientation kCcw = Orientation::counterclockwise;
constexpr Orientation kCw = Orientation::clockwise;
constexpr Orientation kCol = Orientation::collinear;

// Los tests con referencia entera usan coordenadas con |valor| < 2^29: las diferencias
// caben en 30 bits y los productos del determinante en 60, así que int64 es exacto.
constexpr std::int64_t kMaxCoord = std::int64_t{1} << 29;

std::int64_t as_int(double v) { return static_cast<std::int64_t>(v); }
double as_double(std::int64_t v) { return static_cast<double>(v); }
Point2D int_point(std::int64_t x, std::int64_t y) { return {as_double(x), as_double(y)}; }

Orientation sign_of(std::int64_t d) { return d > 0 ? kCcw : (d < 0 ? kCw : kCol); }
Orientation sign_of(double d) { return d > 0.0 ? kCcw : (d < 0.0 ? kCw : kCol); }

// Referencia exacta en enteros.
std::int64_t exact_det(Point2D a, Point2D b, Point2D c) {
    return (as_int(a.x) - as_int(c.x)) * (as_int(b.y) - as_int(c.y)) -
           (as_int(a.y) - as_int(c.y)) * (as_int(b.x) - as_int(c.x));
}

// La fórmula ingenua en double. Sirve para comprobar que los casos de prueba son de
// verdad difíciles: debe equivocarse con frecuencia donde el predicado exacto no.
Orientation naive_orientation(Point2D a, Point2D b, Point2D c) {
    return sign_of((a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x));
}

struct Bezout {
    std::int64_t g;
    std::int64_t x;
    std::int64_t y;  // a*x + b*y = g
};

Bezout extended_gcd(std::int64_t a, std::int64_t b) {
    if (b == 0) return {a, 1, 0};
    const Bezout r = extended_gcd(b, a % b);
    return {r.g, r.y, r.x - (a / b) * r.y};
}

std::int64_t floor_div(std::int64_t a, std::int64_t b) {  // b > 0
    const std::int64_t q = a / b;
    return (a % b != 0 && a < 0) ? q - 1 : q;
}

double nudge(double v, int ulps) {
    for (; ulps > 0; --ulps) v = std::nextafter(v, std::numeric_limits<double>::infinity());
    for (; ulps < 0; ++ulps) v = std::nextafter(v, -std::numeric_limits<double>::infinity());
    return v;
}

// Un predicado exacto respeta las simetrías del determinante: las permutaciones cíclicas
// conservan el signo y las transposiciones lo invierten. Una implementación ingenua las
// viola en casos casi degenerados.
void check_all_symmetries(Point2D a, Point2D b, Point2D c, Orientation expected) {
    CHECK_EQ(orientation(a, b, c), expected);
    CHECK_EQ(orientation(b, c, a), expected);
    CHECK_EQ(orientation(c, a, b), expected);
    CHECK_EQ(orientation(b, a, c), geo::opposite(expected));
    CHECK_EQ(orientation(a, c, b), geo::opposite(expected));
    CHECK_EQ(orientation(c, b, a), geo::opposite(expected));
}

}  // namespace

TEST_CASE("orientation: casos básicos") {
    check_all_symmetries({0, 0}, {1, 0}, {0, 1}, kCcw);
    check_all_symmetries({0, 0}, {0, 1}, {1, 0}, kCw);
    check_all_symmetries({0, 0}, {1, 1}, {2, 2}, kCol);
    check_all_symmetries({0, 0}, {1, 0}, {2, 0}, kCol);   // horizontal
    check_all_symmetries({3, 0}, {3, 1}, {3, -5}, kCol);  // vertical
    check_all_symmetries({-1, -1}, {2, 0}, {0, 3}, kCcw);

    CHECK_EQ(orient2d({0, 0}, {1, 0}, {0, 1}), 1.0);
    CHECK_EQ(orient2d({0, 0}, {2, 0}, {0, 3}), 6.0);
    CHECK_EQ(orient2d({0, 0}, {0, 3}, {2, 0}), -6.0);
}

TEST_CASE("orientation: puntos repetidos son colineales") {
    const Point2D a{1.5, -2.0};
    const Point2D b{3.25, 4.0};
    CHECK_EQ(orientation(a, a, b), kCol);
    CHECK_EQ(orientation(a, b, a), kCol);
    CHECK_EQ(orientation(b, a, a), kCol);
    CHECK_EQ(orientation(a, a, a), kCol);
    CHECK_EQ(orient2d(a, a, b), 0.0);
}

TEST_CASE("orientation: opposite y to_string") {
    CHECK_EQ(geo::opposite(kCcw), kCw);
    CHECK_EQ(geo::opposite(kCw), kCcw);
    CHECK_EQ(geo::opposite(kCol), kCol);
    CHECK_EQ(std::string(geo::to_string(kCcw)), std::string("counterclockwise"));
    CHECK_EQ(std::string(geo::to_string(kCw)), std::string("clockwise"));
    CHECK_EQ(std::string(geo::to_string(kCol)), std::string("collinear"));
}

TEST_CASE("orientation: colineales exactos de gran magnitud dan exactamente cero") {
    std::mt19937_64 rng(7);
    std::uniform_int_distribution<std::int64_t> base(0, (std::int64_t{1} << 28) - 1);
    std::uniform_int_distribution<std::int64_t> small(-(std::int64_t{1} << 13),
                                                      (std::int64_t{1} << 13) - 1);
    for (int i = 0; i < 5000; ++i) {
        const std::int64_t ox = base(rng);
        const std::int64_t oy = base(rng);
        const std::int64_t dx = small(rng);
        const std::int64_t dy = small(rng);
        const std::int64_t m = small(rng);
        const std::int64_t k = small(rng);
        const Point2D a = int_point(ox, oy);
        const Point2D b = int_point(ox + m * dx, oy + m * dy);
        const Point2D c = int_point(ox + k * dx, oy + k * dy);
        REQUIRE(exact_det(a, b, c) == 0);
        CHECK_EQ(orient2d(a, b, c), 0.0);
        CHECK_EQ(orient2d(c, a, b), 0.0);
        CHECK_EQ(orient2d(b, a, c), 0.0);
    }
    // Misma abscisa: colineales sea cual sea el valor de x.
    CHECK_EQ(orient2d({0.1, 1e15}, {0.1, -3.7}, {0.1, 2.5e-15}), 0.0);
}

TEST_CASE("orientation: determinantes exactos ±1 con productos de ~2^56 (casos Bezout)") {
    // Para p, q coprimos existe (r, s) con p*s - q*r = 1 (identidad de Bezout). Con
    // a = o, b = o + (p, q), c = o + (r, s) el determinante exacto vale 1: el triángulo
    // tiene área 1/2 aunque sus lados midan ~2^28. Como los productos intermedios rondan
    // 2^56, el redondeo del double (~2^3) supera al propio determinante y la fórmula
    // ingenua falla en una fracción significativa de los casos (~18 % en la práctica).
    std::mt19937_64 rng(2024);
    std::uniform_int_distribution<std::int64_t> big(std::int64_t{1} << 27,
                                                    (std::int64_t{1} << 28) - 1);
    std::uniform_int_distribution<std::int64_t> base(0, (std::int64_t{1} << 28) - 1);

    int trials = 0;
    int naive_wrong = 0;
    while (trials < 20000) {
        const std::int64_t p = big(rng);
        const std::int64_t q = big(rng);
        const Bezout bz = extended_gcd(p, q);
        if (bz.g != 1) continue;
        ++trials;

        // p*x + q*y = 1  =>  con (r, s) = (-y, x) vale p*s - q*r = 1. Desplazar (r, s) a lo
        // largo de (p, q) no altera el determinante; se elige el múltiplo con r en [0, p).
        const std::int64_t k = -floor_div(-bz.y, p);
        const std::int64_t r = -bz.y + k * p;
        const std::int64_t s = bz.x + k * q;

        const std::int64_t ox = base(rng);
        const std::int64_t oy = base(rng);
        const Point2D a = int_point(ox, oy);
        const Point2D b = int_point(ox + p, oy + q);
        const Point2D c = int_point(ox + r, oy + s);
        REQUIRE(exact_det(a, b, c) == 1);

        check_all_symmetries(a, b, c, kCcw);

        // Trasladar los tres puntos por un vector entero no cambia la orientación.
        const Point2D v = int_point(base(rng), base(rng));
        CHECK_EQ(orientation(a + v, b + v, c + v), kCcw);
        CHECK_EQ(orientation(a - v, b - v, c - v), kCcw);

        if (naive_orientation(a, b, c) != kCcw) ++naive_wrong;
    }
    CHECK_MSG(naive_wrong > trials / 20,
              "la fórmula ingenua acertó demasiado: falló solo " + std::to_string(naive_wrong) +
                  " de " + std::to_string(trials) + " casos");
}

TEST_CASE("orientation: casi colineales aleatorios contra la referencia entera") {
    std::mt19937_64 rng(99);
    std::uniform_int_distribution<std::int64_t> coord(0, kMaxCoord - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_int_distribution<std::int64_t> jitter(-2, 2);
    const auto clamp = [](std::int64_t v) { return std::clamp(v, std::int64_t{0}, kMaxCoord - 1); };

    for (int i = 0; i < 50000; ++i) {
        const std::int64_t ax = coord(rng);
        const std::int64_t ay = coord(rng);
        const std::int64_t bx = coord(rng);
        const std::int64_t by = coord(rng);
        const double t = unit(rng);
        const std::int64_t cx =
            clamp(std::llround(as_double(ax) + t * as_double(bx - ax)) + jitter(rng));
        const std::int64_t cy =
            clamp(std::llround(as_double(ay) + t * as_double(by - ay)) + jitter(rng));

        const Point2D a = int_point(ax, ay);
        const Point2D b = int_point(bx, by);
        const Point2D c = int_point(cx, cy);
        const Orientation expected = sign_of(exact_det(a, b, c));
        CHECK_EQ(orientation(a, b, c), expected);
        CHECK_EQ(orientation(c, a, b), expected);
        CHECK_EQ(orientation(b, a, c), geo::opposite(expected));
    }
}

TEST_CASE("orientation: ejemplo de Kettner et al. perturbado por ulps") {
    // Valores tomados de L. Kettner, K. Mehlhorn, S. Pion, S. Schirra y C. Yap, "Classroom
    // examples of robustness problems in geometric computations", Computational Geometry:
    // Theory and Applications 40(1):61-78, 2008 (solo los datos numéricos del ejemplo).
    // p = (0.5, 0.5), q = (12, 12), r = (24, 24). El determinante exacto es 12*(p.y - p.x),
    // así que al mover p.x e p.y por i y j ulps el signo correcto es el de (j - i). En
    // double, p.x - 24 absorbe la perturbación y la fórmula ingenua responde "colineal".
    const Point2D q{12.0, 12.0};
    const Point2D r{24.0, 24.0};
    int naive_wrong = 0;
    for (int i = -6; i <= 6; ++i) {
        for (int j = -6; j <= 6; ++j) {
            const Point2D p{nudge(0.5, i), nudge(0.5, j)};
            const Orientation expected = j > i ? kCcw : (j < i ? kCw : kCol);
            check_all_symmetries(p, q, r, expected);
            if (naive_orientation(p, q, r) != expected) ++naive_wrong;
        }
    }
    CHECK_MSG(naive_wrong > 0, "la fórmula ingenua acertó en todo el ejemplo de Kettner");
}

TEST_CASE("orientation: antisimetría y simetría cíclica en ternas casi colineales de doubles") {
    std::mt19937_64 rng(31337);
    std::uniform_real_distribution<double> coord(-1e6, 1e6);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_int_distribution<int> ulps(-3, 3);

    for (int i = 0; i < 20000; ++i) {
        const Point2D a{coord(rng), coord(rng)};
        const Point2D b{coord(rng), coord(rng)};
        const double t = unit(rng);
        const Point2D c{nudge(a.x + t * (b.x - a.x), ulps(rng)),
                        nudge(a.y + t * (b.y - a.y), ulps(rng))};
        check_all_symmetries(a, b, c, orientation(a, b, c));
    }
}

// ---- incircle ----------------------------------------------------------------------

namespace {

using geo::incircle;

int sign(double v) { return v > 0.0 ? 1 : (v < 0.0 ? -1 : 0); }

// El filtro de la etapa A sin la cota de error: lo que respondería una implementación
// ingenua. Sirve para comprobar que los casos de prueba son de verdad difíciles.
int naive_incircle_sign(Point2D a, Point2D b, Point2D c, Point2D d) {
    const double adx = a.x - d.x;
    const double ady = a.y - d.y;
    const double bdx = b.x - d.x;
    const double bdy = b.y - d.y;
    const double cdx = c.x - d.x;
    const double cdy = c.y - d.y;
    return sign((adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) +
                (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy) +
                (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady));
}

// Simetrías del determinante 4x4 de las filas (x, y, x² + y², 1): las permutaciones
// cíclicas de (a, b, c) conservan el signo y cualquier transposición de dos de los cuatro
// puntos lo invierte.
void check_incircle_symmetries(Point2D a, Point2D b, Point2D c, Point2D d, int expected) {
    CHECK_EQ(sign(incircle(a, b, c, d)), expected);
    CHECK_EQ(sign(incircle(b, c, a, d)), expected);
    CHECK_EQ(sign(incircle(c, a, b, d)), expected);
    CHECK_EQ(sign(incircle(b, a, c, d)), -expected);
    CHECK_EQ(sign(incircle(a, c, b, d)), -expected);
    CHECK_EQ(sign(incircle(d, b, c, a)), -expected);
    CHECK_EQ(sign(incircle(a, d, c, b)), -expected);
}

}  // namespace

TEST_CASE("incircle: casos básicos") {
    // Circunferencia unidad por a, b, c en sentido antihorario.
    const Point2D a{1, 0};
    const Point2D b{0, 1};
    const Point2D c{-1, 0};
    CHECK_EQ(incircle(a, b, c, {0, 0}), 2.0);  // el centro: dentro
    CHECK(incircle(a, b, c, {0.5, -0.5}) > 0.0);
    CHECK(incircle(a, b, c, {2, 0}) < 0.0);      // fuera
    CHECK(incircle(a, b, c, {1, 1}) < 0.0);
    CHECK_EQ(incircle(a, b, c, {0, -1}), 0.0);   // sobre la circunferencia
    CHECK_EQ(incircle(a, b, c, a), 0.0);         // un vértice está sobre ella

    // En sentido horario el signo se invierte.
    CHECK(incircle(a, c, b, {0, 0}) < 0.0);
    CHECK(incircle(a, c, b, {2, 0}) > 0.0);

    check_incircle_symmetries(a, b, c, {0.25, 0.25}, 1);
    check_incircle_symmetries(a, b, c, {3, -1}, -1);
    check_incircle_symmetries(a, b, c, {0, -1}, 0);
}

TEST_CASE("incircle: rectángulos de doubles arbitrarios son exactamente concíclicos") {
    // Los cuatro vértices de un rectángulo de lados paralelos a los ejes siempre son
    // concíclicos, sean cuales sean sus coordenadas. Con magnitudes muy distintas las
    // diferencias respecto de d no son exactas, de modo que se ejercita la evaluación
    // exacta sin trasladar.
    std::mt19937_64 rng(1234);
    std::uniform_real_distribution<double> mantissa(1.0, 2.0);
    std::uniform_int_distribution<int> exponent(-40, 40);
    std::bernoulli_distribution negative(0.5);
    const auto random_double = [&] {
        const double v = std::ldexp(mantissa(rng), exponent(rng));
        return negative(rng) ? -v : v;
    };

    int naive_wrong = 0;
    for (int i = 0; i < 3000; ++i) {
        const double x1 = random_double();
        const double x2 = random_double();
        const double y1 = random_double();
        const double y2 = random_double();
        if (x1 == x2 || y1 == y2) continue;
        const Point2D a{x1, y1};
        const Point2D b{x2, y1};
        const Point2D c{x2, y2};
        const Point2D d{x1, y2};
        REQUIRE(exact::incircle_sign(a, b, c, d) == 0);
        CHECK_EQ(incircle(a, b, c, d), 0.0);
        CHECK_EQ(incircle(d, b, a, c), 0.0);
        CHECK_EQ(incircle(c, d, b, a), 0.0);
        if (naive_incircle_sign(a, b, c, d) != 0) ++naive_wrong;
    }
    CHECK_MSG(naive_wrong > 0, "la fórmula ingenua acertó en todos los rectángulos");
}

TEST_CASE("incircle: casi concíclicos aleatorios contra la referencia de enteros grandes") {
    // Tres puntos sobre una circunferencia y un cuarto también sobre ella, movido unos
    // pocos ulps. El redondeo de las coordenadas ya rompe la conciclicidad exacta, así que
    // el signo correcto solo se puede saber con aritmética exacta: aquí, la de
    // exact_reference.hpp, que no comparte código con la librería.
    struct Regime {
        double center;
        double radius;
    };
    // Centro lejano con radio pequeño (diferencias exactas: evaluación trasladada) y
    // centro cercano al origen (diferencias inexactas: evaluación sin trasladar).
    const Regime regimes[] = {{1e6, 1.0}, {1e9, 1e3}, {0.5, 1.0}, {1e-3, 1e3}, {3.0, 1e-6}};

    std::mt19937_64 rng(4242);
    std::uniform_real_distribution<double> angle(0.0, 2.0 * std::numbers::pi);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    std::uniform_int_distribution<int> ulps(-3, 3);

    int trials = 0;
    int naive_wrong = 0;
    for (const Regime& reg : regimes) {
        for (int i = 0; i < 4000; ++i) {
            const Point2D center{reg.center * unit(rng), reg.center * unit(rng)};
            const auto on_circle = [&] {
                const double t = angle(rng);
                return center + Point2D{reg.radius * std::cos(t), reg.radius * std::sin(t)};
            };
            const Point2D a = on_circle();
            Point2D b = on_circle();
            Point2D c = on_circle();
            if (orientation(a, b, c) == kCw) std::swap(b, c);
            const Point2D q = on_circle();
            const Point2D d{nudge(q.x, ulps(rng)), nudge(q.y, ulps(rng))};

            const int expected = exact::incircle_sign(a, b, c, d);
            check_incircle_symmetries(a, b, c, d, expected);
            ++trials;
            if (naive_incircle_sign(a, b, c, d) != expected) ++naive_wrong;
        }
    }
    CHECK_MSG(naive_wrong > trials / 100,
              "la fórmula ingenua acertó demasiado: falló solo " + std::to_string(naive_wrong) +
                  " de " + std::to_string(trials) + " casos");
}

TEST_CASE("incircle: puntos enteros de una circunferencia, concíclicos lejos del origen") {
    // Los puntos enteros de x² + y² = 1105² (1105 = 5·13·17 tiene muchas representaciones)
    // trasladados a ~2^40: todas las cuaternas son concíclicas y el determinante debe ser
    // exactamente cero. Como el filtro inicial no puede certificar un cero, cada una pasa
    // por la evaluación exacta sobre las diferencias.
    std::vector<Point2D> circle;
    constexpr std::int64_t r = 1105;
    for (std::int64_t x = -r; x <= r; ++x) {
        const std::int64_t y2 = r * r - x * x;
        const auto y = static_cast<std::int64_t>(std::llround(std::sqrt(as_double(y2))));
        if (y * y != y2) continue;
        circle.push_back(int_point(x, y));
        if (y != 0) circle.push_back(int_point(x, -y));
    }
    REQUIRE(circle.size() >= 32);

    const Point2D center{std::ldexp(1.0, 40) + 3.0, -std::ldexp(1.0, 40) + 7.0};
    for (Point2D& p : circle) p += center;

    std::mt19937_64 rng(5);
    std::uniform_int_distribution<std::size_t> pick(0, circle.size() - 1);
    int checked = 0;
    while (checked < 5000) {
        const Point2D a = circle[pick(rng)];
        const Point2D b = circle[pick(rng)];
        const Point2D c = circle[pick(rng)];
        const Point2D d = circle[pick(rng)];
        if (orientation(a, b, c) == kCol) continue;  // hay puntos repetidos
        ++checked;
        CHECK_EQ(incircle(a, b, c, d), 0.0);

        // Mover d un ulp en x hacia el centro lo deja estrictamente dentro.
        if (d.x == center.x) continue;
        const Point2D inward{nudge(d.x, d.x > center.x ? -1 : 1), d.y};
        const int expected = orientation(a, b, c) == kCcw ? 1 : -1;
        CHECK_EQ(sign(incircle(a, b, c, inward)), expected);
        CHECK_EQ(exact::incircle_sign(a, b, c, inward), expected);
    }
}

GEOTEST_MAIN

#include "geometra/point.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include "geotest.hpp"

using geo::Point2D;

TEST_CASE("Point2D: construcción y valores por defecto") {
    const Point2D origin;
    CHECK_EQ(origin.x, 0.0);
    CHECK_EQ(origin.y, 0.0);

    const Point2D p{1.5, -2.0};
    CHECK_EQ(p.x, 1.5);
    CHECK_EQ(p.y, -2.0);

    const Point2D q{.x = 3.0, .y = 4.0};
    CHECK_EQ(q.x, 3.0);
    CHECK_EQ(q.y, 4.0);
}

TEST_CASE("Point2D: aritmética de vectores") {
    const Point2D a{1.0, 2.0};
    const Point2D b{3.0, -4.0};

    CHECK_EQ(a + b, (Point2D{4.0, -2.0}));
    CHECK_EQ(a - b, (Point2D{-2.0, 6.0}));
    CHECK_EQ(-a, (Point2D{-1.0, -2.0}));
    CHECK_EQ(a * 2.0, (Point2D{2.0, 4.0}));
    CHECK_EQ(2.0 * a, (Point2D{2.0, 4.0}));
    CHECK_EQ(b / 2.0, (Point2D{1.5, -2.0}));

    Point2D c = a;
    c += b;
    CHECK_EQ(c, (Point2D{4.0, -2.0}));
    c -= b;
    CHECK_EQ(c, a);
    c *= 3.0;
    CHECK_EQ(c, (Point2D{3.0, 6.0}));
    c /= 3.0;
    CHECK_EQ(c, a);
}

TEST_CASE("Point2D: igualdad y orden lexicográfico") {
    // Los paréntesis extra protegen las comas de las llaves frente al preprocesador.
    CHECK((Point2D{1.0, 2.0} == Point2D{1.0, 2.0}));
    CHECK((Point2D{1.0, 2.0} != Point2D{2.0, 1.0}));
    CHECK((Point2D{0.0, 0.0} == Point2D{-0.0, 0.0}));  // -0.0 == 0.0 en IEEE 754

    CHECK((Point2D{1.0, 9.0} < Point2D{2.0, 0.0}));  // decide x
    CHECK((Point2D{1.0, 2.0} < Point2D{1.0, 3.0}));  // a igual x, decide y
    CHECK(!(Point2D{1.0, 2.0} < Point2D{1.0, 2.0}));
    CHECK((Point2D{1.0, 2.0} <= Point2D{1.0, 2.0}));
    CHECK((Point2D{2.0, 0.0} > Point2D{1.0, 9.0}));

    std::vector<Point2D> pts{{2, 1}, {1, 2}, {1, 1}, {0, 5}, {2, 0}};
    std::sort(pts.begin(), pts.end());
    const std::vector<Point2D> expected{{0, 5}, {1, 1}, {1, 2}, {2, 0}, {2, 1}};
    CHECK_EQ(pts, expected);
}

TEST_CASE("Point2D: dot, cross, normas y distancias") {
    CHECK_EQ(geo::dot({1, 2}, {3, 4}), 11.0);
    CHECK_EQ(geo::cross({1, 0}, {0, 1}), 1.0);
    CHECK_EQ(geo::cross({0, 1}, {1, 0}), -1.0);
    CHECK_EQ(geo::cross({2, 3}, {4, 6}), 0.0);  // paralelos

    CHECK_EQ(geo::squared_norm({3, 4}), 25.0);
    CHECK_EQ(geo::norm({3, 4}), 5.0);
    CHECK_EQ(geo::squared_distance({1, 1}, {4, 5}), 25.0);
    CHECK_EQ(geo::distance({1, 1}, {4, 5}), 5.0);
    CHECK_EQ(geo::midpoint({0, 0}, {2, 4}), (Point2D{1, 2}));

    // norm usa hypot: no desborda con componentes grandes.
    CHECK(std::isfinite(geo::norm({1e200, 1e200})));
}

TEST_CASE("Point2D: is_finite") {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    CHECK(geo::is_finite({0.0, 0.0}));
    CHECK(geo::is_finite({-1e308, 1e308}));
    CHECK(!geo::is_finite({nan, 0.0}));
    CHECK(!geo::is_finite({0.0, nan}));
    CHECK(!geo::is_finite({inf, 0.0}));
    CHECK(!geo::is_finite({0.0, -inf}));
}

TEST_CASE("Point2D: las operaciones básicas son constexpr") {
    static_assert(Point2D{1, 2} + Point2D{3, 4} == Point2D{4, 6});
    static_assert(Point2D{1, 2} - Point2D{3, 4} == Point2D{-2, -2});
    static_assert(geo::cross({1, 0}, {0, 1}) == 1.0);
    static_assert(geo::dot({1, 2}, {3, 4}) == 11.0);
    static_assert(Point2D{1, 2} < Point2D{1, 3});
    static_assert(geo::squared_distance({0, 0}, {3, 4}) == 25.0);
    CHECK(true);
}

TEST_CASE("Point2D: salida por flujo") {
    std::ostringstream os;
    os << Point2D{1.5, -2.0};
    CHECK_EQ(os.str(), std::string("(1.5, -2)"));
}

GEOTEST_MAIN

#pragma once

// Mini framework de tests de Geometra, sin dependencias externas.
//
//   TEST_CASE("nombre") { CHECK(cond); CHECK_EQ(a, b); REQUIRE(cond); }
//   GEOTEST_MAIN   // al final del archivo
//
// El ejecutable acepta un filtro opcional por subcadena: ./test_x "colineal".
// Devuelve 0 si todos los casos pasan y 1 en caso contrario.

#include <cstdio>
#include <exception>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace geotest {

struct TestCase {
    const char* name;
    void (*fn)();
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> cases;
    return cases;
}

struct Registrar {
    Registrar(const char* name, void (*fn)()) { registry().push_back({name, fn}); }
};

struct State {
    int checks = 0;
    int failures = 0;
    const char* current = "";
};

inline State& state() {
    static State s;
    return s;
}

struct RequireFailure final : std::exception {
    const char* what() const noexcept override { return "REQUIRE falló"; }
};

template <class T>
std::string to_string(const T& value) {
    std::ostringstream os;
    os << std::boolalpha << std::setprecision(17) << value;
    return os.str();
}

template <class T>
std::string to_string(const std::vector<T>& values) {
    std::string out = "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) out += ", ";
        out += to_string(values[i]);
    }
    return out + "]";
}

inline void record(bool ok, const char* file, int line, const char* expr,
                   const std::string& detail) {
    ++state().checks;
    if (ok) return;
    ++state().failures;
    std::fprintf(stderr, "  FALLO %s:%d en \"%s\"\n    %s\n", file, line, state().current, expr);
    if (!detail.empty()) std::fprintf(stderr, "    %s\n", detail.c_str());
}

inline int run(int argc, char** argv) {
    const std::string_view filter = argc > 1 ? argv[1] : "";
    int ran = 0;
    int failed_cases = 0;

    for (const TestCase& tc : registry()) {
        if (!filter.empty() && std::string_view(tc.name).find(filter) == std::string_view::npos) {
            continue;
        }
        state().current = tc.name;
        const int before = state().failures;
        try {
            tc.fn();
        } catch (const RequireFailure&) {
            // ya registrado por REQUIRE
        } catch (const std::exception& e) {
            record(false, "?", 0, "excepción inesperada", e.what());
        } catch (...) {
            record(false, "?", 0, "excepción inesperada", "(desconocida)");
        }
        ++ran;
        const bool ok = state().failures == before;
        failed_cases += ok ? 0 : 1;
        std::printf("[%s] %s\n", ok ? " OK " : "FAIL", tc.name);
    }

    std::printf("\n%d casos, %d aserciones, %d fallos\n", ran, state().checks, state().failures);
    return failed_cases == 0 ? 0 : 1;
}

}  // namespace geotest

#define GEOTEST_CAT_(a, b) a##b
#define GEOTEST_CAT(a, b) GEOTEST_CAT_(a, b)

#define TEST_CASE(name)                                                             \
    static void GEOTEST_CAT(geotest_fn_, __LINE__)();                               \
    static const ::geotest::Registrar GEOTEST_CAT(geotest_reg_, __LINE__)(          \
        name, &GEOTEST_CAT(geotest_fn_, __LINE__));                                 \
    static void GEOTEST_CAT(geotest_fn_, __LINE__)()

#define CHECK(expr) ::geotest::record(static_cast<bool>(expr), __FILE__, __LINE__, #expr, "")

#define CHECK_MSG(expr, msg) \
    ::geotest::record(static_cast<bool>(expr), __FILE__, __LINE__, #expr, (msg))

#define CHECK_EQ(a, b)                                                                  \
    do {                                                                                \
        const auto& geotest_a_ = (a);                                                   \
        const auto& geotest_b_ = (b);                                                   \
        ::geotest::record(geotest_a_ == geotest_b_, __FILE__, __LINE__, #a " == " #b,  \
                          "izquierda: " + ::geotest::to_string(geotest_a_) +            \
                              "\n    derecha:   " + ::geotest::to_string(geotest_b_));  \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                          \
    do {                                                                               \
        const double geotest_a_ = (a);                                                 \
        const double geotest_b_ = (b);                                                 \
        const double geotest_d_ =                                                      \
            geotest_a_ > geotest_b_ ? geotest_a_ - geotest_b_ : geotest_b_ - geotest_a_; \
        ::geotest::record(geotest_d_ <= (tol), __FILE__, __LINE__, #a " ≈ " #b,        \
                          "izquierda: " + ::geotest::to_string(geotest_a_) +           \
                              "\n    derecha:   " + ::geotest::to_string(geotest_b_)); \
    } while (0)

#define REQUIRE(expr)                                                                    \
    do {                                                                                 \
        const bool geotest_ok_ = static_cast<bool>(expr);                                \
        ::geotest::record(geotest_ok_, __FILE__, __LINE__, "REQUIRE(" #expr ")", "");    \
        if (!geotest_ok_) throw ::geotest::RequireFailure{};                             \
    } while (0)

#define CHECK_THROWS_AS(expr, Exception)                                                 \
    do {                                                                                 \
        bool geotest_caught_ = false;                                                    \
        try {                                                                            \
            (void)(expr);                                                                \
        } catch (const Exception&) {                                                     \
            geotest_caught_ = true;                                                      \
        } catch (...) {                                                                  \
        }                                                                                \
        ::geotest::record(geotest_caught_, __FILE__, __LINE__,                           \
                          "CHECK_THROWS_AS(" #expr ", " #Exception ")", "");             \
    } while (0)

#define GEOTEST_MAIN \
    int main(int argc, char** argv) { return ::geotest::run(argc, argv); }

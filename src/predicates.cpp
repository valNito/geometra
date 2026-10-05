// Predicado de orientación con signo exacto.
//
// PROCEDENCIA Y LICENCIA
// ----------------------
// Este archivo es una adaptación a C++ de las rutinas orient2d(), orient2dadapt() y
// fast_expansion_sum_zeroelim() del archivo predicates.c de Jonathan Richard Shewchuk
// (School of Computer Science, Carnegie Mellon University, 18 de mayo de 1996), que su
// autor colocó en el dominio público ("Placed in the public domain by Jonathan Richard
// Shewchuk"; su página web añade "This code is in the public domain"). Las constantes
// kResultErrBound y kCcwErrBound{A,B,C} son las que ese archivo calcula en exactinit() y
// que el artículo citado demuestra. Se conserva deliberadamente la estructura de las
// etapas A-D y la secuencia exacta de operaciones, porque de ellas depende la prueba de
// corrección; la expresión en C++ (tipos, lambdas, std::span, comentarios) es propia.
// Diferencias respecto del original: two_product() usa std::fma en lugar de la partición
// de Dekker (Split), y fast_expansion_sum_zeroelim() no lee más allá de los arreglos.
// Nada de este archivo proviene de Triangle, del mismo autor, cuya licencia sí es
// restrictiva. Los detalles y referencias completas están en NOTICE, en la raíz.
//
// Referencias:
//  - J. R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and Fast Robust
//    Geometric Predicates", Discrete & Computational Geometry 18(3):305-363, 1997.
//    Versión preliminar: informe técnico CMU-CS-96-140, 1996.
//    https://www.cs.cmu.edu/~quake/robust.html
//  - D. E. Knuth, The Art of Computer Programming, vol. 2, §4.2.2, 1969 (TwoSum).
//  - T. J. Dekker, "A floating-point technique for extending the available precision",
//    Numerische Mathematik 18:224-242, 1971 (FastTwoSum).
//  - T. Ogita, S. M. Rump y S. Oishi, "Accurate Sum and Dot Product", SIAM Journal on
//    Scientific Computing 26(6):1955-1988, 2005 (TwoProduct mediante FMA).
//
// FUNCIONAMIENTO
// --------------
// La evaluación es adaptativa: primero se calcula el determinante en double junto con una
// cota de su error de redondeo. Si el valor supera la cota, su signo es correcto y se
// devuelve de inmediato (es el caso habitual). Si no, se recalcula con precisión creciente
// —etapas B y C— y, solo en el peor caso, se obtiene el valor exacto como expansión (D).
//
// REQUISITOS DEL ENTORNO
// ----------------------
// Cada operación en double debe redondearse individualmente según
// IEEE 754. La librería se compila con -ffp-contract=off (sin contracción implícita a FMA)
// y no debe compilarse con -ffast-math ni con precisión extendida x87.

#include "geometra/predicates.hpp"

#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

static_assert(std::numeric_limits<double>::is_iec559, "Geometra requiere doubles IEEE 754.");
static_assert(FLT_EVAL_METHOD == 0 || FLT_EVAL_METHOD == 1,
              "Geometra requiere que las operaciones en double se evalúen en precisión double "
              "(FLT_EVAL_METHOD 0 o 1): la precisión extendida x87 rompe los predicados exactos.");

namespace geo {
namespace {

// Épsilon de Shewchuk: 2^-53, la mitad del épsilon de la máquina.
constexpr double kEpsilon = std::numeric_limits<double>::epsilon() / 2.0;
constexpr double kResultErrBound = (3.0 + 8.0 * kEpsilon) * kEpsilon;
constexpr double kCcwErrBoundA = (3.0 + 16.0 * kEpsilon) * kEpsilon;
constexpr double kCcwErrBoundB = (2.0 + 12.0 * kEpsilon) * kEpsilon;
constexpr double kCcwErrBoundC = (9.0 + 64.0 * kEpsilon) * kEpsilon * kEpsilon;

/// Par (hi, lo) tal que hi + lo es, exactamente, el resultado de una operación.
struct TwoDouble {
    double hi;
    double lo;
};

// ---- Transformaciones libres de error (Knuth, Dekker) --------------------------------

inline TwoDouble two_sum(double a, double b) noexcept {
    const double x = a + b;
    const double bv = x - a;
    const double av = x - bv;
    const double br = b - bv;
    const double ar = a - av;
    return {x, ar + br};
}

/// Versión rápida de two_sum; requiere |a| >= |b|.
inline TwoDouble fast_two_sum(double a, double b) noexcept {
    const double x = a + b;
    const double bv = x - a;
    return {x, b - bv};
}

inline TwoDouble two_diff(double a, double b) noexcept {
    const double x = a - b;
    const double bv = a - x;
    const double av = x + bv;
    const double br = bv - b;
    const double ar = a - av;
    return {x, ar + br};
}

/// Error de la resta x = fl(a - b) ya calculada: a - b == x + two_diff_tail(a, b, x).
inline double two_diff_tail(double a, double b, double x) noexcept {
    const double bv = a - x;
    const double av = x + bv;
    const double br = bv - b;
    const double ar = a - av;
    return ar + br;
}

inline TwoDouble two_product(double a, double b) noexcept {
    const double x = a * b;
    return {x, std::fma(a, b, -x)};
}

/// (a1 + a0) - (b1 + b0) como expansión de 4 componentes, de menor a mayor magnitud.
inline std::array<double, 4> two_two_diff(double a1, double a0, double b1, double b0) noexcept {
    const TwoDouble i0 = two_diff(a0, b0);
    const TwoDouble j0 = two_sum(a1, i0.hi);
    const TwoDouble i1 = two_diff(j0.lo, b1);
    const TwoDouble j1 = two_sum(j0.hi, i1.hi);
    return {i0.lo, i1.lo, j1.lo, j1.hi};
}

/// Suma de dos expansiones no solapadas (componentes de menor a mayor magnitud),
/// eliminando los ceros del resultado. Devuelve cuántos componentes escribió en `h`,
/// siempre al menos uno. `h` debe tener capacidad para e.size() + f.size().
std::size_t fast_expansion_sum_zeroelim(std::span<const double> e, std::span<const double> f,
                                        std::span<double> h) noexcept {
    std::size_t ei = 0;
    std::size_t fi = 0;
    std::size_t hi = 0;

    // Extrae el siguiente componente de menor magnitud entre las cabezas de e y f.
    const auto next = [&]() noexcept {
        const bool take_e =
            fi >= f.size() || (ei < e.size() && std::fabs(e[ei]) <= std::fabs(f[fi]));
        return take_e ? e[ei++] : f[fi++];
    };
    const auto emit = [&](double v) noexcept {
        if (v != 0.0) h[hi++] = v;
    };

    double q = next();
    if (ei < e.size() && fi < f.size()) {
        TwoDouble s = fast_two_sum(next(), q);
        q = s.hi;
        emit(s.lo);
        while (ei < e.size() && fi < f.size()) {
            s = two_sum(q, next());
            q = s.hi;
            emit(s.lo);
        }
    }
    while (ei < e.size()) {
        const TwoDouble s = two_sum(q, e[ei++]);
        q = s.hi;
        emit(s.lo);
    }
    while (fi < f.size()) {
        const TwoDouble s = two_sum(q, f[fi++]);
        q = s.hi;
        emit(s.lo);
    }
    if (q != 0.0 || hi == 0) h[hi++] = q;
    return hi;
}

/// Etapas B, C y D de orient2d. Solo se llega aquí si el filtro inicial (A) no pudo
/// garantizar el signo. `detsum` es la suma de los valores absolutos de ambos productos.
double orient2d_adaptive(Point2D a, Point2D b, Point2D c, double detsum) noexcept {
    const double acx = a.x - c.x;
    const double bcx = b.x - c.x;
    const double acy = a.y - c.y;
    const double bcy = b.y - c.y;

    // Etapa B: productos exactos de las diferencias redondeadas.
    const TwoDouble left = two_product(acx, bcy);
    const TwoDouble right = two_product(acy, bcx);
    const std::array<double, 4> B = two_two_diff(left.hi, left.lo, right.hi, right.lo);

    double det = (B[0] + B[1]) + B[2] + B[3];
    double errbound = kCcwErrBoundB * detsum;
    if (det >= errbound || -det >= errbound) return det;

    // Errores de redondeo de las cuatro diferencias.
    const double acxtail = two_diff_tail(a.x, c.x, acx);
    const double bcxtail = two_diff_tail(b.x, c.x, bcx);
    const double acytail = two_diff_tail(a.y, c.y, acy);
    const double bcytail = two_diff_tail(b.y, c.y, bcy);

    // Si las diferencias fueron exactas, det ya es exacto.
    if (acxtail == 0.0 && acytail == 0.0 && bcxtail == 0.0 && bcytail == 0.0) return det;

    // Etapa C: corrección de primer orden con los errores de las diferencias.
    errbound = kCcwErrBoundC * detsum + kResultErrBound * std::fabs(det);
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if (det >= errbound || -det >= errbound) return det;

    // Etapa D: valor exacto como expansión de hasta 16 componentes.
    std::array<double, 8> C1{};
    std::array<double, 12> C2{};
    std::array<double, 16> D{};

    TwoDouble s = two_product(acxtail, bcy);
    TwoDouble t = two_product(acytail, bcx);
    std::array<double, 4> u = two_two_diff(s.hi, s.lo, t.hi, t.lo);
    const std::size_t c1len = fast_expansion_sum_zeroelim(B, u, C1);

    s = two_product(acx, bcytail);
    t = two_product(acy, bcxtail);
    u = two_two_diff(s.hi, s.lo, t.hi, t.lo);
    const std::size_t c2len = fast_expansion_sum_zeroelim(std::span{C1}.first(c1len), u, C2);

    s = two_product(acxtail, bcytail);
    t = two_product(acytail, bcxtail);
    u = two_two_diff(s.hi, s.lo, t.hi, t.lo);
    const std::size_t dlen = fast_expansion_sum_zeroelim(std::span{C2}.first(c2len), u, D);

    // El componente más significativo de una expansión no solapada tiene su mismo signo.
    return D[dlen - 1];
}

}  // namespace

double orient2d(Point2D a, Point2D b, Point2D c) noexcept {
    const double detleft = (a.x - c.x) * (b.y - c.y);
    const double detright = (a.y - c.y) * (b.x - c.x);
    const double det = detleft - detright;

    // Etapa A: si los productos tienen signos distintos (o uno es cero) la resta no puede
    // cancelar y el signo de `det` es correcto. Si no, se aplica la cota de error.
    double detsum = 0.0;
    if (detleft > 0.0) {
        if (detright <= 0.0) return det;
        detsum = detleft + detright;
    } else if (detleft < 0.0) {
        if (detright >= 0.0) return det;
        detsum = -detleft - detright;
    } else {
        return det;
    }

    const double errbound = kCcwErrBoundA * detsum;
    if (det >= errbound || -det >= errbound) return det;

    return orient2d_adaptive(a, b, c, detsum);
}

}  // namespace geo

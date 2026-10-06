#pragma once

// Referencia exacta para los tests, independiente de la aritmética de expansiones de la
// librería: enteros con signo de precisión arbitraria y conversión exacta de doubles.
// Solo implementa lo que necesitan los tests (suma, resta, producto y signo) y no busca
// eficiencia.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "geometra/point.hpp"

namespace exact {

class BigInt {
public:
    BigInt() = default;

    /// m · 2^shift, con shift >= 0.
    static BigInt from_scaled(std::int64_t m, int shift) {
        BigInt r;
        if (m == 0) return r;
        r.negative_ = m < 0;
        const std::uint64_t mag =
            m < 0 ? std::uint64_t{0} - static_cast<std::uint64_t>(m) : static_cast<std::uint64_t>(m);
        const auto s = static_cast<unsigned>(shift % 32);
        r.limbs_.assign(static_cast<std::size_t>(shift / 32), 0U);
        constexpr std::uint64_t kMask = 0xffffffffULL;
        const std::uint64_t lo = mag & kMask;
        const std::uint64_t hi = mag >> 32;
        if (s == 0) {
            r.limbs_.push_back(static_cast<std::uint32_t>(lo));
            r.limbs_.push_back(static_cast<std::uint32_t>(hi));
        } else {
            r.limbs_.push_back(static_cast<std::uint32_t>((lo << s) & kMask));
            r.limbs_.push_back(static_cast<std::uint32_t>(((lo >> (32 - s)) | (hi << s)) & kMask));
            r.limbs_.push_back(static_cast<std::uint32_t>(hi >> (32 - s)));
        }
        trim(r.limbs_);
        return r;
    }

    [[nodiscard]] int sign() const { return limbs_.empty() ? 0 : (negative_ ? -1 : 1); }

    friend BigInt operator-(BigInt a) {
        if (!a.limbs_.empty()) a.negative_ = !a.negative_;
        return a;
    }

    friend BigInt operator+(const BigInt& a, const BigInt& b) {
        BigInt r;
        if (a.negative_ == b.negative_) {
            r.limbs_ = add(a.limbs_, b.limbs_);
            r.negative_ = a.negative_;
        } else if (compare(a.limbs_, b.limbs_) >= 0) {
            r.limbs_ = sub(a.limbs_, b.limbs_);
            r.negative_ = a.negative_;
        } else {
            r.limbs_ = sub(b.limbs_, a.limbs_);
            r.negative_ = b.negative_;
        }
        if (r.limbs_.empty()) r.negative_ = false;
        return r;
    }

    friend BigInt operator-(const BigInt& a, const BigInt& b) { return a + (-b); }

    friend BigInt operator*(const BigInt& a, const BigInt& b) {
        BigInt r;
        if (a.limbs_.empty() || b.limbs_.empty()) return r;
        r.limbs_ = mul(a.limbs_, b.limbs_);
        r.negative_ = a.negative_ != b.negative_;
        return r;
    }

private:
    using Limbs = std::vector<std::uint32_t>;  // little-endian, sin ceros a la izquierda

    static void trim(Limbs& a) {
        while (!a.empty() && a.back() == 0) a.pop_back();
    }

    static int compare(const Limbs& a, const Limbs& b) {
        if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
        for (std::size_t i = a.size(); i-- > 0;) {
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        }
        return 0;
    }

    static Limbs add(const Limbs& a, const Limbs& b) {
        Limbs r(std::max(a.size(), b.size()) + 1, 0U);
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < r.size(); ++i) {
            const std::uint64_t sum = carry + (i < a.size() ? a[i] : 0U) + (i < b.size() ? b[i] : 0U);
            r[i] = static_cast<std::uint32_t>(sum);
            carry = sum >> 32;
        }
        trim(r);
        return r;
    }

    /// a - b con a >= b.
    static Limbs sub(const Limbs& a, const Limbs& b) {
        Limbs r(a.size(), 0U);
        std::int64_t borrow = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            std::int64_t diff = static_cast<std::int64_t>(a[i]) - borrow -
                                static_cast<std::int64_t>(i < b.size() ? b[i] : 0U);
            borrow = diff < 0 ? 1 : 0;
            if (diff < 0) diff += std::int64_t{1} << 32;
            r[i] = static_cast<std::uint32_t>(diff);
        }
        trim(r);
        return r;
    }

    static Limbs mul(const Limbs& a, const Limbs& b) {
        Limbs r(a.size() + b.size(), 0U);
        for (std::size_t i = 0; i < a.size(); ++i) {
            std::uint64_t carry = 0;
            for (std::size_t j = 0; j < b.size(); ++j) {
                const std::uint64_t cur =
                    r[i + j] + static_cast<std::uint64_t>(a[i]) * b[j] + carry;
                r[i + j] = static_cast<std::uint32_t>(cur);
                carry = cur >> 32;
            }
            r[i + b.size()] = static_cast<std::uint32_t>(carry);
        }
        trim(r);
        return r;
    }

    bool negative_ = false;
    Limbs limbs_;
};

/// Exponente e tal que v = m · 2^e con m entero de a lo sumo 53 bits (v finito, no nulo).
inline int exponent_of(double v) {
    int e = 0;
    (void)std::frexp(v, &e);
    return e - 53;
}

/// v · 2^-base_exp como entero exacto. Requiere base_exp <= exponent_of(v).
inline BigInt scaled(double v, int base_exp) {
    if (v == 0.0) return {};
    int e = 0;
    const double f = std::frexp(v, &e);
    const auto m = static_cast<std::int64_t>(std::ldexp(f, 53));  // exacto: |m| < 2^53
    return BigInt::from_scaled(m, (e - 53) - base_exp);
}

/// Signo exacto del determinante de geo::incircle(a, b, c, d), evaluado en su forma
/// trasladada 3x3 con enteros grandes. Todas las coordenadas se escalan por la misma
/// potencia de 2, lo que multiplica el determinante por un factor positivo.
inline int incircle_sign(geo::Point2D a, geo::Point2D b, geo::Point2D c, geo::Point2D d) {
    const std::array<double, 8> coords{a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y};
    int base = std::numeric_limits<int>::max();
    for (const double v : coords) {
        if (v != 0.0) base = std::min(base, exponent_of(v));
    }
    if (base == std::numeric_limits<int>::max()) return 0;  // los cuatro puntos en el origen

    const auto s = [base](double v) { return scaled(v, base); };
    const BigInt adx = s(a.x) - s(d.x);
    const BigInt ady = s(a.y) - s(d.y);
    const BigInt bdx = s(b.x) - s(d.x);
    const BigInt bdy = s(b.y) - s(d.y);
    const BigInt cdx = s(c.x) - s(d.x);
    const BigInt cdy = s(c.y) - s(d.y);

    const BigInt alift = adx * adx + ady * ady;
    const BigInt blift = bdx * bdx + bdy * bdy;
    const BigInt clift = cdx * cdx + cdy * cdy;
    const BigInt det = alift * (bdx * cdy - cdx * bdy) + blift * (cdx * ady - adx * cdy) +
                       clift * (adx * bdy - bdx * ady);
    return det.sign();
}

}  // namespace exact

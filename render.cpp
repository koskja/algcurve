#include "render.hpp"

using OffsetPolynomial = Polynomial<Polynomial<double, 2>, 2>;
template<typename T>
using Point = std::pair<T, T>;

static usize get_max_individual_degree(OffsetPolynomial& poly) {
    usize max_deg = 0;
    for (const auto& [_, poly_coeff] : poly.coefficients) {
        for (const auto& [mono, _] : poly_coeff.coefficients) {
            if (mono.exponents[0] > max_deg) max_deg = mono.exponents[0];
            if (mono.exponents[1] > max_deg) max_deg = mono.exponents[1];
        }
    }
    return max_deg;
}

usize align_to_simd(usize x) {
    return (x + (SIMD_ALIGN - 1)) & ~(SIMD_ALIGN - 1);
}

void get_powers(std::span<double> powers, double x, usize n, bool aligned_to_simd=false) {
    if (aligned_to_simd) {
        assert((((usize)powers.data()) & (SIMD_ALIGN - 1) == 0));
        assert(powers.size() % SIMD_NUM_VALUES == 0);
    }
    if (n == 0) return;
    powers[0] = 1.0; 
    if (n == 1) return;

    usize remainder = n % SIMD_NUM_VALUES;
    usize to_simd_num = (remainder == 0) ? 0 : SIMD_NUM_VALUES - remainder;
    usize i = 1;

    for (; i < SIMD_NUM_VALUES && i < n; ++i) {
        powers[i] = powers[i - 1] * x;
    }
    if (i == n) return;
    double x_mul_step = powers[SIMD_NUM_VALUES - 1] * x;
    for (; i < n; i += SIMD_NUM_VALUES) {
        for (usize j = 0; j < SIMD_NUM_VALUES; ++j) {
             if (i + j >= n) break;
             powers[i + j] = powers[i + j - SIMD_NUM_VALUES] * x_mul_step;
        }
    }
    for (; i < n; ++i) {
        powers[i] = powers[i - 1] * x;
    }
}

// TODO: use Texture
struct PreparedPowers {
    usize stored_powers;
    double *powers;
    static PreparedPowers from_values(std::span<double> values, usize max_power) {
        return PreparedPowers(values, max_power);
    }
    static PreparedPowers from_range(double min, double max, usize n_values, usize max_power) {
        return PreparedPowers(min, max, n_values, max_power);
    }
    PreparedPowers(double min, double max, usize n_values, usize max_power) {
        stored_powers = align_to_simd(max_power + 1);
        powers = (double*) aligned_alloc(SIMD_ALIGN, stored_powers * n_values * sizeof(double));
        parallel_for(0, n_values, [&](usize i) {
            double value = min + i * (max - min) / (double)n_values;
            std::span<double> current_powers(powers + i * stored_powers, stored_powers);
            get_powers(current_powers, value, max_power + 1, true);
        });
    }
    PreparedPowers(std::span<double> values, usize max_power) {
        stored_powers = align_to_simd(max_power + 1);
        powers = (double*) aligned_alloc(SIMD_ALIGN, stored_powers * values.size() * sizeof(double));
        parallel_for(0, values.size(), [&](usize i) {
            std::span<double> current_powers(powers + i * stored_powers, stored_powers);
            get_powers(current_powers, values[i], max_power + 1, true);
        });
    }
    ~PreparedPowers() {
        free(powers);
    }
    std::span<const double> get(usize index) const {
        return std::span<const double>(powers + index * stored_powers, stored_powers);
    }
};

struct PreparedLattices {
    std::vector<PreparedPowers> powers;
    OffsetPolynomial polynomial;
    double width;
    PreparedLattices(double min, double max, OffsetPolynomial poly, usize max_granularity) 
     : polynomial(poly) {
        width = max - min;
        usize max_deg = get_max_individual_degree(poly);
        powers = parallel_nat_map(max_granularity, [&](usize granularity) {
            usize size = 1 << granularity;
            double step = width / size;
            double start = step / 2;
            double end = width - step / 2;
            return PreparedPowers::from_range(start, end, size, max_deg);
        })
    }
    Polynomial<double, 2> eval(usize granularity, Point<usize> at) {
        auto& p = powers[granularity];
        auto xpowers = p.get(at.first);
        auto ypowers = p.get(at.second);
        polynomial.map([&](Polynomial<double, 2> pcoeff) {
            double sum = 0;
            for (const auto& [mon, cof] : pcoeff.coefficients) {
                sum += cof * xpowers[mon.exponents[0]] * ypowers[mon.exponents[1]];
            }
            return sum;
        })
    }
}

bool may_have_root(const Polynomial<double, 2>& poly, std::span<double> delta_powers) {
    double sum = 0;
    for (const auto& [mon, cof] : poly.coefficients) {
        usize degree = mon.degree();
        if (degree == 0) {
            sum += abs(cof);
            continue;
        }
        sum -= abs(cof) * delta_powers[degree]; // TODO: maybe sum first and multiply afterwards?
    }
    return sum <= 0;
}

std::vector<bool> are_points_viable(std::span<Point<usize>> points, usize granularity, PreparedLattices& lattices) {

}
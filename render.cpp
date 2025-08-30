#include "render.hpp"
#include "thread.hpp"

/// A polynomial in two variables, [x, y], where each coefficient is a polynomial in two variables, [dx, dy].
/// By evaluating the coefficients at a "point" [dx, dy], we get a polynomial in two variables, [x, y].
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

/// Get the first `n` powers of `x` and store them in `powers`.
void get_powers(std::span<double> powers, double x, usize n, bool aligned_to_simd=false) {
    if (aligned_to_simd) {
        assert(((((usize)powers.data()) & (SIMD_ALIGN - 1)) == 0));
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

/// A container for precalculated powers of a list of values. 
/// Essentially a thin wrapper around a 2D array of doubles (a list of rows),
/// where each row is aligned to `SIMD_ALIGN`.
struct PreparedPowers {
    usize num_rows;
    usize row_len;
    double *powers;
    static PreparedPowers from_values(std::span<double> values, usize max_power) {
        return PreparedPowers(values, max_power);
    }
    static PreparedPowers from_range(double min, double max, usize n_values, usize max_power) {
        return PreparedPowers(min, max, n_values, max_power);
    }
    PreparedPowers(double min, double max, usize n_values, usize max_power) : num_rows(n_values) {
        row_len = align_to_simd(max_power + 1);
        powers = (double*) aligned_alloc(SIMD_ALIGN, row_len * n_values * sizeof(double));
        parallel_for( n_values, [&](usize i) {
            double value = min + i * (max - min) / (double)n_values;
            std::span<double> current_powers(powers + i * row_len, row_len);
            get_powers(current_powers, value, max_power + 1, true);
        }, 1024 / (max_power + 1));
    }
    PreparedPowers(std::span<double> values, usize max_power) : num_rows(values.size()) {
        row_len = align_to_simd(max_power + 1);
        powers = (double*) aligned_alloc(SIMD_ALIGN, row_len * values.size() * sizeof(double));
        parallel_for( values.size(), [&](usize i) {
            std::span<double> current_powers(powers + i * row_len, row_len);
            get_powers(current_powers, values[i], max_power + 1, true);
        }, 1024 / (max_power + 1));
    }
    ~PreparedPowers() {
        free(powers);
    }
    std::span<const double> get(usize index) const {
        assert(index < num_rows);
        return std::span<const double>(powers + index * row_len, row_len);
    }
};

/// A container for precalculated data for offsetting a polynomial (calculating the substitution p(x, y) -> p(x - dx, y - dy)).
/// The specific calculation required for checking for roots of a polynomial (`may_have_root`)
/// is done by checking points on a lattice (a square grid of points with a side length of 2^granularity).
/// This struct stores a list of PreparedPowers, one for each granularity.
struct PreparedLattices {
    std::vector<PreparedPowers> powers;
    OffsetPolynomial polynomial;
    /// The width of the bounding box of the lattice. 
    double width;
    PreparedLattices(double min, double max, OffsetPolynomial poly, usize max_granularity) 
     : polynomial(poly) {
        width = max - min;
        usize max_deg = get_max_individual_degree(poly);
        powers = parallel_nat_map<PreparedPowers>(max_granularity, [&](usize granularity) {
            usize size = 1 << granularity;
            double step = width / size;
            double start = step / 2;
            double end = width - step / 2;
            return PreparedPowers::from_range(start, end, size, max_deg);
        });
    }
    Polynomial<double, 2> eval(usize granularity, Point<usize> at) {
        auto& p = powers[granularity];
        auto xpowers = p.get(at.first);
        auto ypowers = p.get(at.second);
        return polynomial.map<double>([&](Polynomial<double, 2> pcoeff) -> double {
            double sum = 0;
            for (const auto& [mon, cof] : pcoeff.coefficients) {
                sum += cof * xpowers[mon.exponents[0]] * ypowers[mon.exponents[1]];
            }
            return sum;
        });
    }
};

/// Check if the polynomial may have a root in the box `[-delta, delta]^2`.
/// Returns false if it is CERTAIN that the polynomial DOES NOT have a root in the box.
/// Returns true if it is POSSIBLE that the polynomial DOES HAVE a root in the box.
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
    auto delta = lattices.width / (1 << granularity);
    auto delta_powers = std::vector<double>(get_max_individual_degree(lattices.polynomial) + 1);
    get_powers(delta_powers, delta, delta_powers.size());
    return parallel_map<bool>([&](Point<usize> point) -> bool {
        auto p = lattices.eval(granularity, point);
        return may_have_root(p, delta_powers);
    }, points);
}

std::vector<Point<usize>> subdivide_viable_points(std::span<Point<usize>> points, std::span<bool> viable, usize old_granularity) {
    assert(viable.size() == points.size());
    auto result = std::vector<Point<usize>>();
    for (usize i = 0; i < points.size(); ++i) {
        if (!viable[i]) continue;
        auto [x, y] = points[i];
        result.emplace_back(x * 2, y * 2);
        result.emplace_back(x * 2 + 1, y * 2);
        result.emplace_back(x * 2, y * 2 + 1);
        result.emplace_back(x * 2 + 1, y * 2 + 1);
    }
    return result;
}
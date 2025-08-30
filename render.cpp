#include "render.hpp"
#include "thread.hpp"

template <typename T> using Point = std::pair<T, T>;

OffsetPolynomial to_offset_polynomial(const Polynomial<double, 2>& poly) {
    // Promote the 2D polynomial p(x, y) to a 4D polynomial p(x, y, z, w)
    Polynomial<double, 4> p4;
    for (const auto& [mon2, cof] : poly.coefficients) {
        Monomial<4> mon4;
        mon4.exponents[0] = mon2.exponents[0];
        mon4.exponents[1] = mon2.exponents[1];
        mon4.exponents[2] = 0;
        mon4.exponents[3] = 0;
        p4.coefficients[mon4] += cof;
    }

    // Build substitutions: x -> x - z, y -> y - w (z = dx, w = dy)
    Polynomial<double, 4> xsub;
    xsub.set(1, 0, 0, 0) = 1.0;  // +x
    xsub.set(0, 0, 1, 0) = -1.0; // -z
    Polynomial<double, 4> ysub;
    ysub.set(0, 1, 0, 0) = 1.0;  // +y
    ysub.set(0, 0, 0, 1) = -1.0; // -w

    // Apply substitutions
    p4 = p4.substitute(0, xsub);
    p4 = p4.substitute(1, ysub);

    // Unnest inner variables (dx=z at index 2, dy=w at index 3 -> becomes index 2 after first unnest)
    auto nested_once = p4.unnest_inner(2);           // split out z
    auto nested_twice = nested_once.unnest_inner(2); // split out w (now at index 2)

    // Merge the two 1D inner polys (z and w) into a single 2D inner poly (dx, dy)
    return merge_coeffs<double, 2>(nested_twice);
}

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
void get_powers(std::span<double> powers, double x, usize n, bool aligned_to_simd) {
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

PreparedPowers::PreparedPowers() : num_rows(0), row_len(0), powers(nullptr) {}

PreparedPowers::PreparedPowers(PreparedPowers&& other) noexcept
    : num_rows(other.num_rows), row_len(other.row_len), powers(other.powers) {
    other.num_rows = 0;
    other.row_len = 0;
    other.powers = nullptr;
}

PreparedPowers& PreparedPowers::operator=(PreparedPowers&& other) noexcept {
    if (this == &other) return *this;
    if (powers) {
        free(powers);
    }
    num_rows = other.num_rows;
    row_len = other.row_len;
    powers = other.powers;
    other.num_rows = 0;
    other.row_len = 0;
    other.powers = nullptr;
    return *this;
}

PreparedPowers PreparedPowers::from_values(std::span<double> values, usize max_power) {
    return PreparedPowers(values, max_power);
}

PreparedPowers PreparedPowers::from_range(double min, double max, usize n_values, usize max_power) {
    return PreparedPowers(min, max, n_values, max_power);
}

PreparedPowers::PreparedPowers(double min, double max, usize n_values, usize max_power) : num_rows(n_values) {
    row_len = align_to_simd(max_power + 1);
    powers = (double *)aligned_alloc(SIMD_ALIGN, row_len * n_values * sizeof(double));
    parallel_for(
        n_values,
        [&](usize i) {
            double value = min + i * (max - min) / (double)n_values;
            std::span<double> current_powers(powers + i * row_len, row_len);
            get_powers(current_powers, value, max_power + 1, true);
        },
        1024 / (max_power + 1));
}

PreparedPowers::PreparedPowers(std::span<double> values, usize max_power) : num_rows(values.size()) {
    row_len = align_to_simd(max_power + 1);
    powers = (double *)aligned_alloc(SIMD_ALIGN, row_len * values.size() * sizeof(double));
    parallel_for(
        values.size(),
        [&](usize i) {
            std::span<double> current_powers(powers + i * row_len, row_len);
            get_powers(current_powers, values[i], max_power + 1, true);
        },
        1024 / (max_power + 1));
}

PreparedPowers::~PreparedPowers() {
    free(powers);
}

std::span<const double> PreparedPowers::get(usize index) const {
    assert(index < num_rows);
    return std::span<const double>(powers + index * row_len, row_len);
}

PreparedLattices::PreparedLattices(double min, double max, usize max_degree, usize max_granularity)
    : max_degree(max_degree) {
    width = max - min;
    powers.reserve(max_granularity + 1);
    for (usize granularity = 0; granularity <= max_granularity; ++granularity) {
        usize size = 1 << granularity;
        double step = width / size;
        double start = step / 2;
        double end = width - step / 2;
        powers.push_back(PreparedPowers::from_range(start + min, end + min, size, max_degree));
    }
}

Polynomial<double, 2> PreparedLattices::eval(usize granularity, Point<usize> at, OffsetPolynomial& poly) {
    auto& p = powers[granularity];
    auto xpowers = p.get(at.first);
    auto ypowers = p.get(at.second);
    assert(get_max_individual_degree(poly) <= max_degree);
    return poly.map<double>([&](Polynomial<double, 2> pcoeff) -> double {
        double sum = 0;
        for (const auto& [mon, cof] : pcoeff.coefficients) {
            sum += cof * xpowers[mon.exponents[0]] * ypowers[mon.exponents[1]];
        }
        return sum;
    });
}

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

std::vector<u8> are_points_viable(std::span<Point<usize>> points,
                                  usize granularity,
                                  PreparedLattices& lattices,
                                  OffsetPolynomial& poly) {
    auto delta = lattices.width / (1 << granularity);
    auto delta_powers = std::vector<double>(lattices.max_degree + 1);
    get_powers(delta_powers, delta, delta_powers.size());
    return parallel_map<u8>(
        [&](Point<usize> point) -> u8 {
            auto p = lattices.eval(granularity, point, poly);
            return may_have_root(p, delta_powers) ? 1 : 0;
        },
        points);
}

std::vector<Point<usize>> subdivide_viable_points(std::span<Point<usize>> points, std::span<u8> viable) {
    assert(viable.size() == points.size());
    auto result = std::vector<Point<usize>>();
    for (usize i = 0; i < points.size(); ++i) {
        if (viable[i] == 0) continue;
        auto [x, y] = points[i];
        result.emplace_back(x * 2, y * 2);
        result.emplace_back(x * 2 + 1, y * 2);
        result.emplace_back(x * 2, y * 2 + 1);
        result.emplace_back(x * 2 + 1, y * 2 + 1);
    }
    return result;
}

Texture2D<BlackWhite> render_image(OffsetPolynomial& poly, PreparedLattices& lattices, ImageParams& params) {
    assert(params.width == params.height);
    auto img = Texture2D<BlackWhite>(params.width, params.height);
    auto max_granularity = lattices.powers.size() - 1;
    assert(max_granularity == std::log2(params.width));
    auto granularity = 0;
    std::vector<Point<usize>> points = {std::pair<usize, usize>(0, 0)};
    while (granularity < max_granularity) {
        auto viable = are_points_viable(points, granularity, lattices, poly);
        points = subdivide_viable_points(points, viable);
        granularity++;
    }
    for (auto& point : points) {
        img(point.first, point.second).set_v(1);
    }
    return img;
}

std::vector<Texture2D<BlackWhite>> render_images(PreparedLattices& lattices, AnimationParams& params) {
    std::cout << "Rendering " << params.lambdas.size() << " images" << std::endl;
    return parallel_map<Texture2D<BlackWhite>>(
        [&](double lambda) {
            auto base = params.p1 * lambda + params.p2 * (1 - lambda);
            auto offset_poly = to_offset_polynomial(base);
            return render_image(offset_poly, lattices, params.image);
        },
        std::span(params.lambdas));
}
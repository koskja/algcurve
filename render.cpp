#include "render.hpp"
#include "paper.hpp"
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

enum OssifiedPolynomialType { DENSE, SPARSE };

OssifiedPolynomialType decide_ossified_polynomial_type(const OffsetPolynomial& poly) {
    if (poly.coefficients.empty()) return SPARSE;

    double total_fullness = 0.0;
    usize num_entries = 0;

    for (const auto& [_, inner] : poly.coefficients) {
        usize xdeg = inner.max_degree(0);
        usize ydeg = inner.max_degree(1);
        usize total_slots = (xdeg + 1) * (ydeg + 1);
        usize used_slots = inner.coefficients.size();
        double fullness = total_slots == 0 ? 0.0 : (double)used_slots / (double)total_slots;
        total_fullness += fullness;
        num_entries += 1;
    }

    double average_fullness = num_entries == 0 ? 0.0 : total_fullness / (double)num_entries;
    return average_fullness > 0.5 ? DENSE : SPARSE;
}

template <typename T> usize align_to_simd(usize x) {
    auto byte_len = x * sizeof(T);
    auto aligned_byte_len = (byte_len + (SIMD_ALIGN - 1)) & ~(SIMD_ALIGN - 1);
    return aligned_byte_len / sizeof(T);
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

PreparedPowers PreparedPowers::from_values(std::span<double> values, usize max_power) {
    return PreparedPowers(values, max_power);
}

PreparedPowers PreparedPowers::from_range(double min, double max, usize n_values, usize max_power) {
    return PreparedPowers(min, max, n_values, max_power);
}

PreparedPowers::PreparedPowers(double min, double max, usize n_values, usize max_power) : num_rows(n_values) {
    row_len = align_to_simd<double>(max_power + 1);
    powers = SimdHeapArray<double, SIMD_ALIGN>(row_len * n_values);
    parallel_for(
        n_values,
        [&](usize i) {
            double value = min + i * (max - min) / (double)n_values;
            auto current_powers = powers.slice_len(i * row_len, row_len);
            get_powers(current_powers, value, max_power + 1, true);
        },
        1024 / (max_power + 1));
}

PreparedPowers::PreparedPowers(std::span<double> values, usize max_power) : num_rows(values.size()) {
    row_len = align_to_simd<double>(max_power + 1);
    powers = SimdHeapArray<double, SIMD_ALIGN>(row_len * values.size());
    parallel_for(
        values.size(),
        [&](usize i) {
            auto current_powers = powers.slice_len(i * row_len, row_len);
            get_powers(current_powers, values[i], max_power + 1, true);
        },
        1024 / (max_power + 1));
}

std::span<const double> PreparedPowers::get(usize index) const {
    assert(index < num_rows);
    return powers.slice_len(index * row_len, row_len);
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

template <typename P>
    requires std::is_base_of_v<OssifiedPolynomial, P>
Polynomial<double, 2>
PreparedLattices::eval(usize granularity, Point<usize> at, const OssifiedOffsetPolynomial<P>& poly) {
    auto& p = powers[granularity];
    auto xpowers = p.get(at.first);
    auto ypowers = p.get(at.second);
    return poly.template map<double>([&](const P& pcoeff) -> double { return pcoeff.eval(xpowers, ypowers); });
}

template <typename P>
    requires std::is_base_of_v<OssifiedPolynomial, P>
std::vector<u8> are_points_viable(std::span<Point<usize>> points,
                                  usize granularity,
                                  PreparedLattices& lattices,
                                  OssifiedOffsetPolynomial<P>& poly) {
    auto delta = lattices.width / (2 << granularity);
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

template <typename P>
    requires std::is_base_of_v<OssifiedPolynomial, P>
Texture2D<BlackWhite> render_image(OssifiedOffsetPolynomial<P>& poly, PreparedLattices& lattices, ImageParams& params) {
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
            auto type = decide_ossified_polynomial_type(offset_poly);
            if (type == DENSE) {
                auto dense_poly = offset_poly.map<DenseOssifiedPolynomial>(
                    [](const Polynomial<double, 2>& p) -> DenseOssifiedPolynomial {
                        return DenseOssifiedPolynomial(p);
                    });
                return render_image(dense_poly, lattices, params.image);
            } else {
                auto sparse_poly = offset_poly.map<SparseOssifiedPolynomial>(
                    [](const Polynomial<double, 2>& p) -> SparseOssifiedPolynomial {
                        return SparseOssifiedPolynomial(p);
                    });
                return render_image(sparse_poly, lattices, params.image);
            }
        },
        std::span(params.lambdas));
}

DenseOssifiedPolynomial::DenseOssifiedPolynomial() {
    x_degree = 0;
    y_degree = 0;
    grid = SimdHeapArray<double, SIMD_ALIGN>();
}

DenseOssifiedPolynomial::DenseOssifiedPolynomial(Polynomial<double, 2> polynomial) {
    x_degree = polynomial.max_degree(0);
    y_degree = polynomial.max_degree(1);
    const usize xdim = x_degree + 1;
    const usize ydim = y_degree + 1;
    const usize total = xdim * ydim;
    grid = SimdHeapArray<double, SIMD_ALIGN>(total);
    for (usize i = 0; i < total; ++i) grid[i] = 0.0;

    for (const auto& [monomial, coefficient] : polynomial.coefficients) {
        const usize xi = static_cast<usize>(monomial.exponents[0]);
        const usize yi = static_cast<usize>(monomial.exponents[1]);
        grid[xi * ydim + yi] += coefficient;
    }
}

double DenseOssifiedPolynomial::eval(std::span<const double> xpowers, std::span<const double> ypowers) const {
    // Compute: sum_j xpowers[j] * dot(grid_row_j, ypowers)

    const double *__restrict xp = xpowers.data();
    const double *__restrict yp = ypowers.data();

    // assert(reinterpret_cast<uintptr_t>(xp) % SIMD_ALIGN == 0 && "xp is not SIMD aligned");
    // assert(reinterpret_cast<uintptr_t>(yp) % SIMD_ALIGN == 0 && "yp is not SIMD aligned");

    xp = (const double *)__builtin_assume_aligned(xp, SIMD_ALIGN);
    yp = (const double *)__builtin_assume_aligned(yp, SIMD_ALIGN);

    const usize ydim = y_degree + 1;
    double total = 0.0;
    for (usize j = 0; j <= x_degree; ++j) {
        const double xj = xp[j];
        const double *__restrict row = &grid[j * ydim];
        double dot = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(enable) interleave(enable)
#endif
        for (usize k = 0; k <= y_degree; ++k) {
            dot += row[k] * yp[k];
        }
        total += xj * dot;
    }
    return total;
}

SparseOssifiedPolynomial::SparseOssifiedPolynomial() {}

SparseOssifiedPolynomial::SparseOssifiedPolynomial(Polynomial<double, 2> polynomial) {
    num_coefficients = polynomial.coefficients.size();
    coefficients = SimdHeapArray<double, SIMD_ALIGN>(num_coefficients);
    x_exponents = SimdHeapArray<exp_t, SIMD_ALIGN>(num_coefficients);
    y_exponents = SimdHeapArray<exp_t, SIMD_ALIGN>(num_coefficients);
    usize i = 0;
    for (const auto& [monomial, coefficient] : polynomial.coefficients) {
        coefficients[i] = coefficient;
        x_exponents[i] = monomial.exponents[0];
        y_exponents[i] = monomial.exponents[1];
        i++;
    }
}

double SparseOssifiedPolynomial::eval(std::span<const double> xpowers, std::span<const double> ypowers) const {
    const double *__restrict xp = xpowers.data();
    const double *__restrict yp = ypowers.data();
    const double *__restrict cof = coefficients.data;
    const exp_t *__restrict xe = x_exponents.data;
    const exp_t *__restrict ye = y_exponents.data;

    // assert(reinterpret_cast<uintptr_t>(xp) % SIMD_ALIGN == 0 && "xp is not SIMD aligned");
    // assert(reinterpret_cast<uintptr_t>(yp) % SIMD_ALIGN == 0 && "yp is not SIMD aligned");
    // assert(reinterpret_cast<uintptr_t>(cof) % SIMD_ALIGN == 0 && "coefficients is not SIMD aligned");
    // assert(reinterpret_cast<uintptr_t>(xe) % SIMD_ALIGN == 0 && "x_exponents is not SIMD aligned");
    // assert(reinterpret_cast<uintptr_t>(ye) % SIMD_ALIGN == 0 && "y_exponents is not SIMD aligned");

    xp = (const double *)__builtin_assume_aligned(xp, SIMD_ALIGN);
    yp = (const double *)__builtin_assume_aligned(yp, SIMD_ALIGN);
    cof = (const double *)__builtin_assume_aligned(cof, SIMD_ALIGN);
    xe = (const exp_t *)__builtin_assume_aligned(xe, SIMD_ALIGN);
    ye = (const exp_t *)__builtin_assume_aligned(ye, SIMD_ALIGN);

    double total = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(enable) interleave(enable)
#endif
    for (usize i = 0; i < num_coefficients; ++i) {
        total += cof[i] * xp[xe[i]] * yp[ye[i]];
    }
    return total;
}
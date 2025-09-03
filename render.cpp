#include "render.hpp"
#include "paper.hpp"
#include "thread.hpp"

template <typename T> using Point = std::pair<T, T>;

OffsetPolynomial to_offset_polynomial(const HashmapPolynomial<double, 2>& poly) {
    // Promote the 2D polynomial p(x, y) to a 4D polynomial p(x, y, z, w)
    HashmapPolynomial<double, 4> p4;
    for (const auto& [mon2, cof] : poly.coefficients) {
        Monomial<4> mon4;
        mon4.exponents[0] = mon2.exponents[0];
        mon4.exponents[1] = mon2.exponents[1];
        mon4.exponents[2] = 0;
        mon4.exponents[3] = 0;
        p4.coefficients[mon4] += cof;
    }

    // Build substitutions: x -> x - z, y -> y - w (z = dx, w = dy)
    HashmapPolynomial<double, 4> xsub;
    xsub.set(Monomial<4>(1, 0, 0, 0), 1.0);  // +x
    xsub.set(Monomial<4>(0, 0, 1, 0), -1.0); // -z
    HashmapPolynomial<double, 4> ysub;
    ysub.set(Monomial<4>(0, 1, 0, 0), 1.0);  // +y
    ysub.set(Monomial<4>(0, 0, 0, 1), -1.0); // -w

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
        usize xdeg = inner.degree(0);
        usize ydeg = inner.degree(1);
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
            double value = n_values <= 1 ? min : (min + i * (max - min) / (double)(n_values - 1));
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

HashmapPolynomial<double, 2>
PreparedLattices::eval(usize granularity, Point<usize> at, const OssifiedOffsetPolynomial& poly) {
    auto& p = powers[granularity];
    auto xpowers = p.get(at.first);
    auto ypowers = p.get(at.second);
    return poly.map<double>([&](const OssifiedPolynomial<double, 2>& pcoeff) -> double {
        return std::visit([&](const auto& p) { return p.eval_with_precalculated_powers({xpowers, ypowers}); }, pcoeff);
    });
}

std::vector<u8> are_points_viable(std::span<Point<usize>> points,
                                  usize granularity,
                                  PreparedLattices& lattices,
                                  const OssifiedOffsetPolynomial& poly) {
    auto box_width = lattices.width / (1 << granularity);
    auto delta = box_width / 4;
    auto num_powers = lattices.max_degree * 2 + 4;
    auto delta_powers = std::vector<double>(num_powers);
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

SparseTexture<2, BlackWhite>
render_image(const OssifiedOffsetPolynomial& poly, PreparedLattices& lattices, ImageParams& params) {
    assert(params.width == params.height);
    auto img = SparseTexture<2, BlackWhite>(params.width, params.height);
    auto max_granularity = lattices.powers.size() - 1;
    assert((1 << max_granularity) == params.width);
    usize granularity = 0;
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

std::vector<SparseTexture<2, BlackWhite>> render_images(PreparedLattices& lattices, AnimationParams& params) {
    std::cout << "Rendering " << params.lambdas.size() << " images" << std::endl;
    return parallel_map<SparseTexture<2, BlackWhite>>(
        [&](double lambda) {
            auto base = params.p1 * lambda + params.p2 * (1 - lambda);
            auto offset_poly = to_offset_polynomial(base);
            auto type = decide_ossified_polynomial_type(offset_poly);
            if (type == DENSE) {
                auto dense_poly = offset_poly.map<OssifiedPolynomial<double, 2>>(
                    [](const HashmapPolynomial<double, 2>& p) -> OssifiedPolynomial<double, 2> {
                        return DenseOssifiedPolynomial<double, 2>(p);
                    });
                return render_image(dense_poly, lattices, params.image);
            } else {
                auto sparse_poly = offset_poly.map<OssifiedPolynomial<double, 2>>(
                    [](const HashmapPolynomial<double, 2>& p) -> OssifiedPolynomial<double, 2> {
                        return SparseOssifiedPolynomial<double, 2>(p);
                    });
                return render_image(sparse_poly, lattices, params.image);
            }
        },
        std::span(params.lambdas));
}
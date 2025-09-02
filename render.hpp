#include "image.hpp"
#include "polynomial.hpp"
#include <cstdlib>

struct ImageParams {
    usize width, height;
};

struct AnimationParams {
    ImageParams image;
    HashmapPolynomial<double, 2> p1, p2;
    std::span<double> lambdas;
};

usize align_to_simd(usize x);
void get_powers(std::span<double> powers, double x, usize n, bool aligned_to_simd = false);

/// A container for precalculated powers of a list of values.
/// Essentially a thin wrapper around a 2D array of doubles (a list of rows),
/// where each row is aligned to `SIMD_ALIGN`.
struct PreparedPowers {
    usize num_rows;
    usize row_len;
    SimdHeapArray<double, SIMD_ALIGN> powers;
    PreparedPowers();
    static PreparedPowers from_values(std::span<double> values, usize max_power);
    static PreparedPowers from_range(double min, double max, usize n_values, usize max_power);
    PreparedPowers(double min, double max, usize n_values, usize max_power);
    PreparedPowers(std::span<double> values, usize max_power);
    std::span<const double> get(usize index) const;
};

template <typename T, usize NVARS>
using OssifiedPolynomialVariant = std::variant<SparseOssifiedPolynomial<T, NVARS>, DenseOssifiedPolynomial<T, NVARS>>;

template <typename T, usize NVARS>
using OssifiedPolynomial = std::variant<SparseOssifiedPolynomial<T, NVARS>, DenseOssifiedPolynomial<T, NVARS>>;

using OssifiedOffsetPolynomial = HashmapPolynomial<OssifiedPolynomial<double, 2>, 2>;
using OffsetPolynomial = HashmapPolynomial<HashmapPolynomial<double, 2>, 2>;

/// A container for precalculated data for offsetting a polynomial (calculating the substitution p(x, y) -> p(x - dx, y
/// - dy)). The specific calculation required for checking for roots of a polynomial (`may_have_root`) is done by
/// checking points on a lattice (a square grid of points with a side length of 2^granularity). This struct stores a
/// list of PreparedPowers, one for each granularity.
struct PreparedLattices {
    std::vector<PreparedPowers> powers;
    usize max_degree;
    double width;
    PreparedLattices(double min, double max, usize max_degree, usize max_granularity);
    HashmapPolynomial<double, 2>
    eval(usize granularity, std::pair<usize, usize> at, const OssifiedOffsetPolynomial& poly);
};

std::vector<Texture2D<BlackWhite>> render_images(PreparedLattices& lattices, AnimationParams& params);

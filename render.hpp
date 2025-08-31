#include "image.hpp"
#include "polynomial.hpp"
#include <cstdlib>

#define SIMD_SIZE (512 / 8)
#define SIMD_ALIGN SIMD_SIZE
#define SIMD_NUM_VALUES (SIMD_SIZE / sizeof(double))

struct ImageParams {
    usize width, height;
};

struct AnimationParams {
    ImageParams image;
    Polynomial<double, 2> p1, p2;
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

struct OssifiedPolynomial {
    virtual double eval(std::span<const double> xpowers, std::span<const double> ypowers) const = 0;
};

struct SparseOssifiedPolynomial : OssifiedPolynomial {
    SimdHeapArray<double, SIMD_ALIGN> coefficients;
    SimdHeapArray<exp_t, SIMD_ALIGN> x_exponents;
    SimdHeapArray<exp_t, SIMD_ALIGN> y_exponents;
    usize num_coefficients;
    SparseOssifiedPolynomial();
    SparseOssifiedPolynomial(Polynomial<double, 2> polynomial);

    double eval(std::span<const double> xpowers, std::span<const double> ypowers) const override;
};

struct DenseOssifiedPolynomial : OssifiedPolynomial {
    usize x_degree;
    usize y_degree;
    SimdHeapArray<double, SIMD_ALIGN> grid;

    DenseOssifiedPolynomial();
    DenseOssifiedPolynomial(Polynomial<double, 2> polynomial);

    double eval(std::span<const double> xpowers, std::span<const double> ypowers) const override;
};

template <typename P>
    requires std::is_base_of_v<OssifiedPolynomial, P>
using OssifiedOffsetPolynomial = Polynomial<P, 2>;
using OffsetPolynomial = Polynomial<Polynomial<double, 2>, 2>;

/// A container for precalculated data for offsetting a polynomial (calculating the substitution p(x, y) -> p(x - dx, y
/// - dy)). The specific calculation required for checking for roots of a polynomial (`may_have_root`) is done by
/// checking points on a lattice (a square grid of points with a side length of 2^granularity). This struct stores a
/// list of PreparedPowers, one for each granularity.
struct PreparedLattices {
    std::vector<PreparedPowers> powers;
    usize max_degree;
    double width;
    PreparedLattices(double min, double max, usize max_degree, usize max_granularity);
    template <typename P>
        requires std::is_base_of_v<OssifiedPolynomial, P>
    Polynomial<double, 2> eval(usize granularity, std::pair<usize, usize> at, OssifiedOffsetPolynomial<P>& poly);
};

std::vector<Texture2D<BlackWhite>> render_images(PreparedLattices& lattices, AnimationParams& params);

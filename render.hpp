#include "image.hpp"
#include "polynomial.hpp"
#include <cstdlib>

#define SIMD_SIZE 512 / 8
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

using OffsetPolynomial = Polynomial<Polynomial<double, 2>, 2>;

usize align_to_simd(usize x);
void get_powers(std::span<double> powers, double x, usize n, bool aligned_to_simd = false);

/// A container for precalculated powers of a list of values.
/// Essentially a thin wrapper around a 2D array of doubles (a list of rows),
/// where each row is aligned to `SIMD_ALIGN`.
struct PreparedPowers {
    usize num_rows;
    usize row_len;
    double *powers;
    PreparedPowers();
    PreparedPowers(const PreparedPowers&) = delete;
    PreparedPowers& operator=(const PreparedPowers&) = delete;
    PreparedPowers(PreparedPowers&& other) noexcept;
    PreparedPowers& operator=(PreparedPowers&& other) noexcept;
    static PreparedPowers from_values(std::span<double> values, usize max_power);
    static PreparedPowers from_range(double min, double max, usize n_values, usize max_power);
    PreparedPowers(double min, double max, usize n_values, usize max_power);
    PreparedPowers(std::span<double> values, usize max_power);
    ~PreparedPowers();
    std::span<const double> get(usize index) const;
};

/// A container for precalculated data for offsetting a polynomial (calculating the substitution p(x, y) -> p(x - dx, y
/// - dy)). The specific calculation required for checking for roots of a polynomial (`may_have_root`) is done by
/// checking points on a lattice (a square grid of points with a side length of 2^granularity). This struct stores a
/// list of PreparedPowers, one for each granularity.
struct PreparedLattices {
    std::vector<PreparedPowers> powers;
    usize max_degree;
    double width;
    PreparedLattices(double min, double max, usize max_degree, usize max_granularity);
    Polynomial<double, 2> eval(usize granularity, std::pair<usize, usize> at, OffsetPolynomial& poly);
};

OffsetPolynomial to_offset_polynomial(const Polynomial<double, 2>& poly);
Texture2D<BlackWhite> render_image(OffsetPolynomial& poly, PreparedLattices& lattices, ImageParams& params);
std::vector<Texture2D<BlackWhite>> render_images(PreparedLattices& lattices, AnimationParams& params);

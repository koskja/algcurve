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

/// Calculate the powers of `x` up to `n - 1`.
void get_powers(std::span<double> powers, double x, usize n, bool aligned_to_simd = false);

/// A container for precalculated powers of a list of values.
/// Essentially a thin wrapper around a 2D array of doubles (a list of rows),
/// where each row is aligned to `SIMD_ALIGN`.
struct PreparedPowers {
    usize num_rows;
    usize row_len;
    SimdHeapArray<double> powers;
    PreparedPowers();
    static PreparedPowers from_values(std::span<double> values, usize max_power);
    static PreparedPowers from_range(double min, double max, usize n_values, usize max_power);
    PreparedPowers(double min, double max, usize n_values, usize max_power);
    PreparedPowers(std::span<double> values, usize max_power);
    std::span<const double> get(usize index) const;
};

/// A polynomial with 2 variables (x, y) whose coefficients are polynomials with 2 variables (dx, dy).
/// Used for representing p(x - dx, y - dy).
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
    Polynomial<double, 2> eval(usize granularity, std::pair<usize, usize> at, const OffsetPolynomial& poly);
};

/// Render a single image of a polynomial.
Image render_image(const OffsetPolynomial& poly, PreparedLattices& lattices, ImageParams& params);
/// Render a sequence of images to create an animation.
std::vector<Image> render_images(PreparedLattices& lattices, AnimationParams& params);

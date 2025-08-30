#include "polynomial.hpp"
#include "image.hpp"

#define SIMD_SIZE 512 / 8
#define SIMD_ALIGN SIMD_SIZE
#define SIMD_NUM_VALUES (SIMD_SIZE / sizeof(double))

struct ImageParams {
    std::function<std::array<double, 2>(usize, usize)> to_plane;
    std::function<double(double)> soft_clamp;
    usize width, height;
};

struct AnimationParams {
    ImageParams image;
    Polynomial<double, 2> p1, p2;
    std::span<double> lambdas;
};

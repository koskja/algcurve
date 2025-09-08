#pragma once

#include "polynomial.hpp"

// #define DESINGULARIZED_ROOT_CHECK
#define MAGIC_CONSTANT 20.0

bool default_may_have_root(const Polynomial<double, 2>& poly, std::span<double> delta_powers);
bool desingularized_may_have_root(const Polynomial<double, 2>& poly, std::span<double> delta_powers);

#ifdef DESINGULARIZED_ROOT_CHECK
#define DELTA_FACTOR 0.5
#define may_have_root desingularized_may_have_root
#else
#define DELTA_FACTOR 1.0
#define may_have_root default_may_have_root
#endif
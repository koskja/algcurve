#pragma once

#include "polynomial.hpp"

// #define DESINGULARIZED_ROOT_CHECK
#define MAGIC_CONSTANT 20.0

bool default_may_have_root(const Polynomial<double, 2>& poly, std::span<double> delta_powers);
bool desingularized_may_have_root(const Polynomial<double, 2>& poly, std::span<double> delta_powers);

#ifdef DESINGULARIZED_ROOT_CHECK
#define DELTA_FACTOR 0.5
#define _may_have_root desingularized_may_have_root
#else
#define DELTA_FACTOR 1.0
#define _may_have_root default_may_have_root
#endif

/// Check if the polynomial may have a root in the box `[-delta, delta]^2`.
/// Returns false if it is CERTAIN that the polynomial DOES NOT have a root in the box.
/// Returns true if it is POSSIBLE that the polynomial DOES HAVE a root in the box.
#define may_have_root _may_have_root
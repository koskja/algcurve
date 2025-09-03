#pragma once

#include "polynomial.hpp"

bool default_may_have_root(const HashmapPolynomial<double, 2>& poly, std::span<double> delta_powers);
bool desingularized_may_have_root(const HashmapPolynomial<double, 2>& poly, std::span<double> delta_powers);
#define may_have_root desingularized_may_have_root
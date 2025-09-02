#pragma once

#include "polynomial.hpp"

bool may_have_root(const HashmapPolynomial<double, 2>& poly, std::span<double> delta_powers);
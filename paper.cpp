#include "paper.hpp"

/// Check if the polynomial may have a root in the box `[-delta, delta]^2`.
/// Returns false if it is CERTAIN that the polynomial DOES NOT have a root in the box.
/// Returns true if it is POSSIBLE that the polynomial DOES HAVE a root in the box.
bool may_have_root(const Polynomial<double, 2>& poly, std::span<double> delta_powers) {
    double sum = 0;
    for (const auto& [mon, cof] : poly.coefficients) {
        usize degree = mon.degree();
        if (degree == 0) {
            sum += abs(cof);
            continue;
        }
        sum -= abs(cof) * delta_powers[degree]; // TODO: maybe sum first and multiply afterwards?
    }
    return sum <= 0;
}
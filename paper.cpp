#include "paper.hpp"
#include "polynomial.hpp"

/// Check if the polynomial may have a root in the box `[-delta, delta]^2`.
/// Returns false if it is CERTAIN that the polynomial DOES NOT have a root in the box.
/// Returns true if it is POSSIBLE that the polynomial DOES HAVE a root in the box.
bool may_have_root(const HashmapPolynomial<double, 2>& poly, std::span<double> delta_powers) {
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

bool desingularized_may_have_root(const HashmapPolynomial<double, 2>& f, std::span<double> delta_powers) {
    auto phi = f * f;
    auto d = f.degree();
    auto f_inf = f.map_with_exponent<double>(
        [&](const Monomial<2>& mon, const double& coeff) { return (double)(d - mon.degree()) / d * coeff; });
    assert(f_inf.degree() == d - 1);
    auto fx = f.partial_derivative(0);
    auto fy = f.partial_derivative(1);
    auto psi = f_inf * f_inf + fx * fx + fy * fy;
    auto phi_sparse = SparseOssifiedPolynomial<double, 2>(phi);
    auto psi_sparse = SparseOssifiedPolynomial<double, 2>(psi);
    auto phi_slices = phi_sparse.get_single_degree_slices();
    auto psi_slices = psi_sparse.get_single_degree_slices();
    // phi[i] = sum_{j from 0 to i} psi[j] * g[i-j]
    // g[i] = (phi[i] - sum_{j from 0 to (i-1)} psi[j] * g[i-j]) / psi[0]
    auto d2 = 2 * d;
    assert(phi_slices.size() == d2);
    assert(psi_slices.size() == d2 - 2);
    std::vector<HashmapPolynomial<double, 2>> gs(d2 + 1);
    std::vector<HashmapPolynomial<double, 2>> psis(d2 + 1);
    auto psi0 = psi.coefficients[Monomial<2>(0, 0)];
    if (psi0 == 0) {
        return true;
    }
    for (usize i = 0; i < d2; ++i) {
        auto phi_slice = phi_slices[i];
        auto psi_slice = psi_slices[i];
        for (usize j = 0; j < phi_slice.coefficients.size(); ++j) {
            auto phi_mon = Monomial<2>(phi_slice.exponents[0][j], phi_slice.exponents[1][j]);
            gs[phi_mon.degree()].coefficients[phi_mon] += phi_slice.coefficients[j];
            auto psi_mon = Monomial<2>(psi_slice.exponents[0][j], psi_slice.exponents[1][j]);
            psis[psi_mon.degree()].coefficients[psi_mon] += psi_slice.coefficients[j];
        }
    }
    return false;
}
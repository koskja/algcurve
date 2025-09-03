#include "paper.hpp"
#include "polynomial.hpp"

/// Check if the polynomial may have a root in the box `[-delta, delta]^2`.
/// Returns false if it is CERTAIN that the polynomial DOES NOT have a root in the box.
/// Returns true if it is POSSIBLE that the polynomial DOES HAVE a root in the box.
bool default_may_have_root(const HashmapPolynomial<double, 2>& poly, std::span<double> delta_powers) {
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
    if (f_inf.degree() != d - 1) {
        auto newfinf = HashmapPolynomial<double, 2>();
        for (const auto& [mon, coeff] : f_inf.coefficients) {
            if (mon.degree() < d) {
                newfinf.coefficients[mon] = coeff;
            }
        }
        f_inf = newfinf;
    }
    assert(f_inf.degree() == d - 1);
    auto fx = f.partial_derivative(0);
    auto fy = f.partial_derivative(1);
    auto psi = f_inf * f_inf + fx * fx + fy * fy;
    auto phi_sparse = SparseOssifiedPolynomial<double, 2>(phi);
    auto psi_sparse = SparseOssifiedPolynomial<double, 2>(psi);
    auto phi_slices = phi_sparse.get_single_degree_slices();
    auto psi_slices = psi_sparse.get_single_degree_slices();
    // phi[i] = sum_{j from 0 to i} psi[j] * g[i-j]
    // g[i] = (phi[i] - sum_{j from 1 to i} psi[j] * g[i-j]) / psi[0]
    auto d2 = 2 * d;
    assert(phi_slices.size() == (d2 + 1));
    assert(psi_slices.size() == (d2 - 1));
    std::vector<HashmapPolynomial<double, 2>> gs(d2 + 1);
    std::vector<HashmapPolynomial<double, 2>> psis(d2 - 1);
    auto psi0 = psi.coefficients[Monomial<2>(0, 0)];
    if (psi0 == 0) {
        return true;
    }
    // g = phi / psi, this does not change the function and sets the constant term of psi to 1
    psi *= 1.0 / psi0;
    phi *= 1.0 / psi0;
    for (usize i = 0; i <= d2; ++i) {
        auto phi_slice = phi_slices[i];
        for (usize j = 0; j < phi_slice.coefficients.size(); ++j) {
            auto phi_mon = phi_slice.get_monomial(j);
            assert(i == phi_mon.degree());
            gs[i].coefficients[phi_mon] += phi_slice.coefficients[j];
        }
    }
    for (usize i = 0; i <= d2 - 2; ++i) {
        auto psi_slice = psi_slices[i];
        for (usize j = 0; j < psi_slice.coefficients.size(); ++j) {
            auto psi_mon = psi_slice.get_monomial(j);
            assert(i == psi_mon.degree());
            psis[i].coefficients[psi_mon] += psi_slice.coefficients[j];
        }
    }
    for (usize i = 1; i <= d2 - 2; ++i) {
        for (usize j = 0; j + 1 <= i; ++j) {
            gs[i] -= psis[j] * gs[i - j];
        }
    }
    for (usize i = d2 - 1; i <= d2; ++i) {
        for (usize j = 0; j <= d2 - 2; ++j) {
            gs[i] -= psis[j] * gs[i - j];
        }
    }
    std::vector<double> cap_gs(d2 + 1);
    for (usize i = 0; i <= d2; ++i) {
        for (const auto& [mon, cof] : gs[i].coefficients) {
            cap_gs[i] += abs(cof);
        }
    }
    auto cap_m = 0.0;
    for (usize i = 1; i <= d2 - 2; ++i) {
        for (const auto& [mon, cof] : psis[i].coefficients) {
            cap_m += abs(cof) * delta_powers[mon.degree()];
        }
    }
    if (cap_m >= 1.0) {
        return true;
    }
    auto cap_k1 = 0.0;
    auto accum = cap_gs[0];
    for (usize i = 1; i <= d2; ++i) {
        auto val = cap_gs[i] * delta_powers[i];
        accum -= val;
        if (i <= d2 - 2) {
            cap_k1 = std::max(cap_k1, val);
        }
    }
    auto cap_g2dp1 = d2 * cap_k1 * cap_m * (1 - cap_m);
    accum -= cap_g2dp1;
    if (std::isnan(accum)) {
        __builtin_trap(); // triggers a breakpoint or aborts execution if nan
    }
    return accum <= 0;
}
#pragma once

#include "core.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#define PARENTHESIS_PRINT_MIN_EXPONENT 10

char get_var_name(usize i);

/// A monomial in NVARS variables.
/// The field `exponents` stores the exponents of the variables, i. e.
/// the monomial of `exponents = {1, 2, 3}` is `x^1 * y^2 * z^3`.
template <usize NVARS> struct Monomial {
    std::array<exp_t, NVARS> exponents;
    constexpr Monomial() : exponents() {}
    template <typename... Args> constexpr Monomial(Args... exponents) : exponents({static_cast<exp_t>(exponents)...}) {}
    constexpr exp_t degree() const {
        exp_t d = 0;
        for (exp_t e : exponents) {
            d += e;
        }
        return d;
    }
    constexpr exp_t operator[](usize i) const {
        return exponents[i];
    }
    constexpr bool operator==(const Monomial& other) const {
        return exponents == other.exponents;
    }
    constexpr Monomial<NVARS> operator*(const Monomial& other) const {
        Monomial<NVARS> result;
        for (usize i = 0; i < NVARS; ++i) {
            result.exponents[i] = exponents[i] + other.exponents[i];
        }
        return result;
    }
    template <typename T> constexpr T eval(std::span<const T> values) const {
        T result = 1;
        for (usize i = 0; i < NVARS; ++i) {
            result *= std::pow(values[i], exponents[i]);
        }
        return result;
    }
    constexpr std::pair<Monomial<1>, Monomial<NVARS - 1>> split(usize split_variable) const {
        Monomial<1> outer = {this->exponents[split_variable]};
        Monomial<NVARS - 1> inner;
        usize from, to;
        for (from = 0, to = 0; from < NVARS; ++from) {
            if (from == split_variable) continue;
            inner.exponents[to++] = this->exponents[from];
        }
        return {outer, inner};
    }
    operator std::string() const {
        std::string result;
        for (usize i = 0; i < NVARS; ++i) {
            if (exponents[i] != 0) {
                result += get_var_name(i);
                if (exponents[i] > 1) {
                    result += "^";
                    if (exponents[i] >= PARENTHESIS_PRINT_MIN_EXPONENT) {
                        result += "(" + std::to_string(exponents[i]) + ")";
                    } else {
                        result += std::to_string(exponents[i]);
                    }
                }
            }
        }
        return result;
    }
    constexpr Monomial(const std::array<exp_t, NVARS>& e) : exponents(e) {}
    constexpr Monomial(const exp_t (&e)[NVARS]) : exponents() {
        std::copy(std::begin(e), std::end(e), exponents.begin());
    }
};

template <usize NVARS> struct std::hash<Monomial<NVARS>> {
    usize operator()(const Monomial<NVARS>& m) const {
        if (sizeof(usize) >= NVARS * sizeof(exp_t)) {
            usize h = 0;
            constexpr usize bitshift = sizeof(exp_t) * 8;
            for (exp_t e : m.exponents) {
                h = (h << bitshift) | e;
            }
            return h;
        }
        usize h = 0;
        for (exp_t e : m.exponents) {
            h = h * 31 + std::hash<exp_t>{}(e);
        }
        return h;
    }
};

/// A polynomial represented as a map from monomials to coefficients.
template <typename T, usize NVARS> struct HashmapPolynomial {
    std::unordered_map<Monomial<NVARS>, T> coefficients;
    T *get(const Monomial<NVARS>& monomial) {
        auto it = coefficients.find(monomial);
        return it == coefficients.end() ? nullptr : &it->second;
    }
    T *get_or_insert(const Monomial<NVARS>& monomial) {
        return &coefficients[monomial];
    }
    void set(const Monomial<NVARS>& monomial, const T& value) {
        coefficients[monomial] = value;
    }
    exp_t degree() const {
        exp_t max_deg = 0;
        for (const auto& [mon, _] : coefficients) {
            max_deg = std::max(max_deg, mon.degree());
        }
        return max_deg;
    }
    exp_t degree(usize variable) const {
        exp_t max_degree = 0;
        for (const auto& [monomial, _] : coefficients) {
            max_degree = std::max(max_degree, monomial.exponents[variable]);
        }
        return max_degree;
    }
    T eval(std::array<T, NVARS> values) const {
        return eval(std::span<const T>(values));
    }
    T eval(std::span<const T> values) const {
        T result = 0;
        for (const auto& [monomial, coefficient] : coefficients) {
            result += coefficient * monomial.eval(values);
        }
        return result;
    }
    /// Evaluate the polynomial using precalculated powers.
    /// `powers[i]` is the span of the powers of the `i`-th variable.
    T eval_with_precalculated_powers(std::array<std::span<const T>, NVARS> powers) const {
        T result = 0;
        for (const auto& [monomial, coefficient] : coefficients) {
            T term = coefficient;
            for (usize i = 0; i < NVARS; ++i) {
                term *= powers[i][monomial.exponents[i]];
            }
            result += term;
        }
        return result;
    }
    void iterate(const std::function<void(const Monomial<NVARS>&, const T&)>& func) const {
        for (const auto& [monomial, coefficient] : coefficients) {
            func(monomial, coefficient);
        }
    }
    template <typename U> HashmapPolynomial<U, NVARS> map(const std::function<U(const T&)>& func) const {
        HashmapPolynomial<U, NVARS> result;
        iterate([&](const Monomial<NVARS>& monomial, const T& coefficient) {
            result.coefficients.emplace(std::make_pair(monomial, func(coefficient)));
        });
        return result;
    }
    template <typename U>
    HashmapPolynomial<U, NVARS>
    map_with_exponent(const std::function<U(const Monomial<NVARS>&, const T&)>& func) const {
        HashmapPolynomial<U, NVARS> result;
        iterate([&](const Monomial<NVARS>& monomial, const T& coefficient) {
            result.coefficients.emplace(std::make_pair(monomial, func(monomial, coefficient)));
        });
        return result;
    }
    HashmapPolynomial<T, NVARS> partial_derivative(usize variable) const {
        HashmapPolynomial<T, NVARS> result;
        for (const auto& [monomial, coefficient] : coefficients) {
            auto var_deg = monomial.exponents[variable];
            if (var_deg == 0) continue;
            Monomial<NVARS> nmon = monomial;
            nmon.exponents[variable]--;
            result.coefficients.emplace(std::make_pair(nmon, coefficient * var_deg));
        }
        return result;
    }
    HashmapPolynomial<T, NVARS> operator+(const HashmapPolynomial<T, NVARS>& other) const {
        HashmapPolynomial<T, NVARS> result = *this;
        for (const auto& [monomial, coefficient] : other.coefficients) {
            result.coefficients[monomial] += coefficient;
        }
        return result;
    }
    HashmapPolynomial<T, NVARS>& operator+=(const HashmapPolynomial<T, NVARS>& other) {
        for (const auto& [monomial, coefficient] : other.coefficients) {
            this->coefficients[monomial] += coefficient;
        }
        return *this;
    }
    HashmapPolynomial<T, NVARS> operator-(const HashmapPolynomial<T, NVARS>& other) const {
        HashmapPolynomial<T, NVARS> result = *this;
        for (const auto& [monomial, coefficient] : other.coefficients) {
            result.coefficients[monomial] -= coefficient;
        }
        return result;
    }
    HashmapPolynomial<T, NVARS>& operator-=(const HashmapPolynomial<T, NVARS>& other) {
        for (const auto& [monomial, coefficient] : other.coefficients) {
            this->coefficients[monomial] -= coefficient;
        }
        return *this;
    }
    HashmapPolynomial<T, NVARS> operator*(const HashmapPolynomial<T, NVARS>& other) const {
        HashmapPolynomial<T, NVARS> result;
        for (const auto& [monomial1, coefficient1] : coefficients) {
            for (const auto& [monomial2, coefficient2] : other.coefficients) {
                result.coefficients[monomial1 * monomial2] += coefficient1 * coefficient2;
            }
        }
        return result;
    }
    HashmapPolynomial<T, NVARS> operator*(const T& scalar) const {
        HashmapPolynomial<T, NVARS> result = *this;
        for (auto& [monomial, coefficient] : result.coefficients) {
            coefficient *= scalar;
        }
        return result;
    }
    HashmapPolynomial<T, NVARS>& operator*=(const T& scalar) {
        for (auto& [monomial, coefficient] : coefficients) {
            coefficient *= scalar;
        }
        return *this;
    }
    HashmapPolynomial<T, NVARS> pow(exp_t exponent) const {
        assert(exponent >= 0 && "Exponent must be non-negative");
        if (exponent == 0) {
            HashmapPolynomial<T, NVARS> result;
            result.coefficients[Monomial<NVARS>{}] = T{1};
            return result;
        }
        if (exponent == 1) {
            return *this;
        }
        auto half = this->pow(exponent / 2);
        auto squared = half * half;
        if (exponent % 2 == 0) {
            return squared;
        }
        return squared * (*this);
    }
    std::unordered_map<exp_t, HashmapPolynomial<T, NVARS>> get_powers(std::span<const exp_t> exponents) const {
        assert(std::is_sorted(exponents.begin(), exponents.end()) && "Exponents must be sorted in ascending order");
        std::unordered_map<exp_t, HashmapPolynomial<T, NVARS>> result;
        std::function<HashmapPolynomial<T, NVARS>(exp_t)> calculate_power;
        calculate_power = [&result, this, &calculate_power](exp_t e) {
            assert(e >= 0 && "Exponent must be non-negative");
            if (e == 0) {
                HashmapPolynomial<T, NVARS> power;
                power.coefficients[Monomial<NVARS>()] = T{1};
                return power;
            }
            if (e == 1) {
                return *this;
            }
            if (result.contains(e)) {
                return result[e];
            }
            auto half = calculate_power(e / 2);
            auto squared = half * half;
            auto res = squared;
            if (e % 2 != 0) {
                res = res * (*this);
            }
            result[e] = res;
            return res;
        };
        for (usize i = 0; i < exponents.size(); ++i) {
            auto exponent = exponents[i];
            result[exponent] = calculate_power(exponent);
        }
        return result;
    }
    constexpr HashmapPolynomial<T, NVARS> substitute(usize variable, HashmapPolynomial<T, NVARS> inner_val) const {
        HashmapPolynomial<T, NVARS> result;
        std::vector<exp_t> degrees;
        for (const auto& [monomial, _] : coefficients) {
            degrees.push_back(monomial.exponents[variable]);
        }
        std::sort(degrees.begin(), degrees.end());
        degrees.erase(std::unique(degrees.begin(), degrees.end()), degrees.end());
        auto powers = inner_val.get_powers(degrees);
        for (const auto& [monomial, coefficient] : coefficients) {
            auto var_exp = monomial.exponents[variable];
            auto inner_val_power = powers[var_exp];
            auto monomial_without_var = monomial;
            monomial_without_var.exponents[variable] = 0; // Represent only the variables that are not substituted
            for (const auto& [inner_monomial, inner_coefficient] : inner_val_power.coefficients) {
                result.coefficients[monomial_without_var * inner_monomial] += inner_coefficient * coefficient;
            }
        }
        return result;
    }
    operator std::string() const {
        std::string result;
        auto it = coefficients.begin();
        auto end = coefficients.end();
        while (it != end) {
            const auto& [monomial, coefficient] = *it;
            if (coefficient == 0) {
                ++it;
                continue;
            }
            if (it != coefficients.begin() && coefficient > 0) {
                result += "+";
            }
            if (coefficient == std::floor(coefficient)) {
                result += std::to_string((long long)coefficient);
                if (monomial != Monomial<NVARS>()) {
                    result += "*" + std::string(monomial);
                }
            } else {
                result += std::to_string(coefficient);
                if (monomial != Monomial<NVARS>()) {
                    result += "*" + std::string(monomial);
                }
            }
            ++it;
        }
        return result;
    }
    friend std::ostream& operator<<(std::ostream& os, const HashmapPolynomial<T, NVARS>& poly) {
        return os << std::string(poly);
    }
    T& operator[](const Monomial<NVARS>& monomial) {
        return coefficients[monomial];
    }
    const T& operator[](const Monomial<NVARS>& monomial) const {
        return coefficients.at(monomial);
    }
    HashmapPolynomial<HashmapPolynomial<T, NVARS - 1>, 1> unnest_outer(usize outer_variable) const {
        HashmapPolynomial<HashmapPolynomial<T, NVARS - 1>, 1> result;
        for (const auto& [monomial, coefficient] : coefficients) {
            auto [outer, inner] = monomial.split(outer_variable);
            result.coefficients[outer].coefficients[inner] += coefficient;
        }
        return result;
    }
    HashmapPolynomial<HashmapPolynomial<T, 1>, NVARS - 1> unnest_inner(usize inner_variable) const {
        HashmapPolynomial<HashmapPolynomial<T, 1>, NVARS - 1> result;
        for (const auto& [monomial, coefficient] : coefficients) {
            auto [inner, outer] = monomial.split(inner_variable);
            result.coefficients[outer].coefficients[inner] += coefficient;
        }
        return result;
    }
};

template <typename T, usize NVARS> struct HashmapPolynomial;
template <typename T, usize NVARS> struct SparseOssifiedPolynomial;
template <typename T, usize NVARS> struct DenseOssifiedPolynomial;
template <typename T, usize NVARS> struct SparseOssifiedSlice;

/// A polynomial represented by a list of coefficients and corresponding exponents.
/// This representation is efficient for sparse polynomials (polynomials with few non-zero coefficients).
/// Once created, it is immutable. The exponents are stored as an array of structures.
template <typename T, usize NVARS> struct SparseOssifiedPolynomial {
    SimdHeapArray<T> coefficients;
    std::array<SimdHeapArray<exp_t>, NVARS> exponents;
    usize num_coefficients;
    exp_t _degree;
    bool sorted;
    std::array<exp_t, NVARS> _degrees_per_var;

    SparseOssifiedPolynomial() : num_coefficients(0), _degree(0) {
        _degrees_per_var.fill(0);
        sorted = false;
    }
    SparseOssifiedPolynomial(const HashmapPolynomial<T, NVARS>& polynomial) {
        num_coefficients = polynomial.coefficients.size();
        coefficients = SimdHeapArray<T>(num_coefficients);
        for (size_t i = 0; i < NVARS; ++i) {
            exponents[i] = SimdHeapArray<exp_t>(num_coefficients);
        }

        _degree = 0;
        _degrees_per_var.fill(0);

        usize i = 0;
        for (const auto& [monomial, coefficient] : polynomial.coefficients) {
            coefficients[i] = coefficient;
            exp_t current_degree = 0;
            for (size_t v = 0; v < NVARS; ++v) {
                exp_t e = monomial.exponents[v];
                exponents[v][i] = e;
                _degrees_per_var[v] = std::max(_degrees_per_var[v], e);
                current_degree += e;
            }
            _degree = std::max(_degree, current_degree);
            i++;
        }
        sorted = false;
    }

    void sort_coeffs() {
        if (num_coefficients == 0 || sorted) return;

        std::vector<usize> p(num_coefficients);
        std::iota(p.begin(), p.end(), 0);

        std::sort(p.begin(), p.end(), [&](usize i, usize j) {
            exp_t i_degree = 0;
            exp_t j_degree = 0;
            for (usize v = 0; v < NVARS; ++v) {
                i_degree += exponents[v][i];
                j_degree += exponents[v][j];
            }
            if (i_degree < j_degree) return true;
            if (i_degree > j_degree) return false;
            for (usize v = 0; v < NVARS; ++v) {
                if (exponents[v][i] < exponents[v][j]) return true;
                if (exponents[v][i] > exponents[v][j]) return false;
            }
            return false;
        });

        SimdHeapArray<T> new_coefficients(num_coefficients);
        std::array<SimdHeapArray<exp_t>, NVARS> new_exponents;
        for (usize v = 0; v < NVARS; ++v) {
            new_exponents[v] = SimdHeapArray<exp_t>(num_coefficients);
        }

        for (usize i = 0; i < num_coefficients; ++i) {
            new_coefficients[i] = coefficients[p[i]];
            for (usize v = 0; v < NVARS; ++v) {
                new_exponents[v][i] = exponents[v][p[i]];
            }
        }

        coefficients = std::move(new_coefficients);
        exponents = std::move(new_exponents);
        sorted = true;
    }
    std::vector<SparseOssifiedSlice<T, NVARS>> get_single_degree_slices() {
        sort_coeffs();
        std::vector<SparseOssifiedSlice<T, NVARS>> slices;

        if (num_coefficients == 0) {
            for (exp_t d = 0; d <= _degree; ++d) {
                slices.push_back(get_slice(0, 0));
            }
            return slices;
        }

        exp_t running_degree = 0;
        usize start_index = 0;
        for (usize i = 0; i < num_coefficients; ++i) {
            exp_t current_degree = 0;
            for (usize v = 0; v < NVARS; ++v) {
                current_degree += exponents[v][i];
            }

            if (current_degree > running_degree) {
                slices.push_back(get_slice(start_index, i)); // Slice for running_degree

                // Fill in empty slices for missing degrees
                for (exp_t d = running_degree + 1; d < current_degree; ++d) {
                    slices.push_back(get_slice(i, i));
                }

                start_index = i;
                running_degree = current_degree;
            }
        }

        // Push the very last slice of coefficients
        slices.push_back(get_slice(start_index, num_coefficients));

        // Add empty slices up to the polynomial's degree
        for (exp_t d = running_degree + 1; d <= _degree; ++d) {
            slices.push_back(get_slice(num_coefficients, num_coefficients));
        }

        return slices;
    }
    SparseOssifiedSlice<T, NVARS> get_slice(usize start, usize end) {
        if (!sorted) sort_coeffs();
        std::span<const T> coeffs = coefficients.slice(start, end);
        std::array<std::span<const exp_t>, NVARS> exps;
        for (usize v = 0; v < NVARS; ++v) {
            exps[v] = exponents[v].slice(start, end);
        }
        return SparseOssifiedSlice<T, NVARS>{coeffs, exps};
    }

    T *get(const Monomial<NVARS>&) {
        throw std::runtime_error("Not mutable");
    }
    T *get_or_insert(const Monomial<NVARS>&) {
        throw std::runtime_error("Not mutable");
    }
    void set(const Monomial<NVARS>&, const T&) {
        throw std::runtime_error("Not mutable");
    }
    T eval(std::array<T, NVARS> /*values*/) const {
        throw std::runtime_error("Not implemented, use eval_with_precalculated_powers");
    }

    T eval_with_precalculated_powers(std::array<std::span<const T>, NVARS> powers) const {
        T total = T{0};
        for (usize i = 0; i < num_coefficients; ++i) {
            T term = coefficients[i];
            for (usize v = 0; v < NVARS; ++v) {
                term *= powers[v][exponents[v][i]];
            }
            total += term;
        }
        return total;
    }

    exp_t degree() const {
        return _degree;
    }
    exp_t degree(usize variable) const {
        return _degrees_per_var[variable];
    }

    HashmapPolynomial<T, NVARS> to_hashmap() const {
        HashmapPolynomial<T, NVARS> result;
        for (usize i = 0; i < num_coefficients; ++i) {
            Monomial<NVARS> mon;
            for (usize v = 0; v < NVARS; ++v) {
                mon.exponents[v] = exponents[v][i];
            }
            if (coefficients[i] != T{0}) {
                result.coefficients[mon] = coefficients[i];
            }
        }
        return result;
    }
    void iterate(const std::function<void(const Monomial<NVARS>&, const T&)>& func) const {
        usize num_coefficients = coefficients.size();
        for (usize i = 0; i < num_coefficients; ++i) {
            Monomial<NVARS> mon;
            for (usize v = 0; v < NVARS; ++v) {
                mon.exponents[v] = exponents[v][i];
            }
            func(mon, coefficients[i]);
        }
    }
    template <typename U> SparseOssifiedPolynomial<U, NVARS> map(const std::function<U(const T&)>& func) const {
        SparseOssifiedPolynomial<U, NVARS> result;
        result.coefficients = SimdHeapArray<U>(coefficients.size());
        for (usize i = 0; i < NVARS; ++i) {
            result.exponents[i] = SimdHeapArray<exp_t>(exponents[i].size());
        }
        auto num_coefficients = coefficients.size();
        for (usize i = 0; i < num_coefficients; ++i) {
            Monomial<NVARS> mon;
            for (usize v = 0; v < NVARS; ++v) {
                mon.exponents[v] = exponents[v][i];
            }
            auto val = func(coefficients[i]);
            result.coefficients[i] = val;
            for (usize v = 0; v < NVARS; ++v) {
                result.exponents[v][i] = exponents[v][i];
            }
        }
        result._degree = this->_degree;
        result._degrees_per_var = this->_degrees_per_var;
        result.sorted = false;
        return result;
    }
};

template <typename T, usize NVARS> struct SparseOssifiedSlice {
    std::span<const T> coefficients;
    std::array<std::span<const exp_t>, NVARS> exponents;
    Monomial<NVARS> get_monomial(usize i) const {
        Monomial<NVARS> mon;
        for (usize v = 0; v < NVARS; ++v) {
            mon.exponents[v] = exponents[v][i];
        }
        return mon;
    }
};

/// A polynomial represented by a dense grid of coefficients.
/// This representation is efficient for dense polynomials (polynomials with many non-zero coefficients up to a certain
/// degree). Once created, it is immutable.
template <typename T, usize NVARS> struct DenseOssifiedPolynomial {
    std::array<exp_t, NVARS> _degrees_per_var;
    exp_t _degree;
    SimdHeapArray<T> grid;

    DenseOssifiedPolynomial() : _degree(0) {
        _degrees_per_var.fill(0);
    }
    DenseOssifiedPolynomial(const HashmapPolynomial<T, NVARS>& polynomial) {
        _degree = polynomial.degree();
        usize total_size = 1;
        for (usize i = 0; i < NVARS; ++i) {
            _degrees_per_var[i] = polynomial.degree(i);
            total_size *= (_degrees_per_var[i] + 1);
        }

        grid = SimdHeapArray<T>(total_size);

        std::array<usize, NVARS> strides;
        if constexpr (NVARS > 0) {
            strides[NVARS - 1] = 1;
            for (int i = NVARS - 2; i >= 0; --i) {
                strides[i] = strides[i + 1] * (_degrees_per_var[i + 1] + 1);
            }
        }

        for (const auto& [monomial, coefficient] : polynomial.coefficients) {
            usize index = 0;
            for (usize i = 0; i < NVARS; ++i) {
                index += monomial.exponents[i] * strides[i];
            }
            grid[index] = coefficient;
        }
    }

    T *get(const Monomial<NVARS>&) {
        throw std::runtime_error("Not mutable");
    }
    T *get_or_insert(const Monomial<NVARS>&) {
        throw std::runtime_error("Not mutable");
    }
    void set(const Monomial<NVARS>&, const T&) {
        throw std::runtime_error("Not mutable");
    }
    T eval(std::array<T, NVARS> /*values*/) const {
        throw std::runtime_error("Not implemented, use eval_with_precalculated_powers");
    }

    T eval_with_precalculated_powers(std::array<std::span<const T>, NVARS> powers) const {
        if (grid.size() == 0) return T{0};

        if constexpr (NVARS == 0) return grid.size() > 0 ? grid[0] : T{0};

        if constexpr (NVARS == 1) {
            T total = T{0};
            for (usize i = 0; i <= _degrees_per_var[0]; ++i) {
                total += grid[i] * powers[0][i];
            }
            return total;
        }

        if constexpr (NVARS == 2) {
            const T *__restrict p0 = powers[0].data();
            const T *__restrict p1 = powers[1].data();
            p0 = (const T *)__builtin_assume_aligned(p0, SIMD_ALIGN);
            p1 = (const T *)__builtin_assume_aligned(p1, SIMD_ALIGN);

            const usize dim1 = _degrees_per_var[1] + 1;
            T total = T{0};
            for (usize j = 0; j <= _degrees_per_var[0]; ++j) {
                const T *__restrict row = &grid[j * dim1];
                T dot = T{0};
#if defined(__clang__)
#pragma clang loop vectorize(enable) interleave(enable)
#endif
                for (usize k = 0; k <= _degrees_per_var[1]; ++k) {
                    dot += row[k] * p1[k];
                }
                total += p0[j] * dot;
            }
            return total;
        }

        if constexpr (NVARS > 2) {
            throw std::runtime_error("DenseOssifiedPolynomial::eval not implemented for NVARS > 2");
        }
    }

    exp_t degree() const {
        return _degree;
    }
    exp_t degree(usize variable) const {
        return _degrees_per_var[variable];
    }

    HashmapPolynomial<T, NVARS> to_hashmap() const {
        HashmapPolynomial<T, NVARS> result;
        if (grid.size() == 0) return result;
        std::array<exp_t, NVARS> current_exponents;
        current_exponents.fill(0);

        std::function<void(usize)> recurse = [&](usize var_idx) {
            if (var_idx == NVARS) {
                std::array<usize, NVARS> strides;
                if constexpr (NVARS > 0) {
                    strides[NVARS - 1] = 1;
                    for (int i = NVARS - 2; i >= 0; --i) {
                        strides[i] = strides[i + 1] * (_degrees_per_var[i + 1] + 1);
                    }
                }
                usize index = 0;
                for (usize i = 0; i < NVARS; ++i) {
                    index += current_exponents[i] * strides[i];
                }

                if (grid[index] != T{0}) {
                    result.coefficients[Monomial<NVARS>(current_exponents)] = grid[index];
                }
                return;
            }

            for (exp_t e = 0; e <= _degrees_per_var[var_idx]; ++e) {
                current_exponents[var_idx] = e;
                recurse(var_idx + 1);
            }
        };
        recurse(0);
        return result;
    }
    usize to_index(std::array<exp_t, NVARS> indices) const {
        usize index = 0;
        for (usize i = 0; i < NVARS; ++i) {
            index += indices[i] * (this->_degrees_per_var[i] + 1);
        }
        return index;
    }
    void iterate(const std::function<void(const Monomial<NVARS>&, const T&)>& func) const {
        std::array<exp_t, NVARS> indices;
        indices.fill(0);
        while (true) {
            func(Monomial<NVARS>(indices), grid[to_index(indices)]);
            indices[0]++;
            for (usize i = 0; i < NVARS; ++i) {
                if (indices[i] > _degrees_per_var[i]) {
                    indices[i] = 0;
                    if (i == NVARS - 1) {
                        return;
                    }
                    indices[i + 1]++;
                }
            }
        }
    }
    template <typename U> DenseOssifiedPolynomial<U, NVARS> map(const std::function<U(const T&)>& func) const {
        DenseOssifiedPolynomial<U, NVARS> result;
        result.grid = SimdHeapArray<U>(this->grid.size());
        for (usize i = 0; i < this->grid.size(); ++i) {
            result.grid[i] = func(this->grid[i]);
        }
        result._degree = this->_degree;
        result._degrees_per_var = this->_degrees_per_var;
        return result;
    }
};

template <typename T, usize NDIM>
HashmapPolynomial<HashmapPolynomial<T, 2>, NDIM>
merge_coeffs(HashmapPolynomial<HashmapPolynomial<HashmapPolynomial<T, 1>, 1>, NDIM> nested_poly) {
    HashmapPolynomial<HashmapPolynomial<T, 2>, NDIM> result;

    for (const auto& [outer_monomial, inner_poly] : nested_poly.coefficients) {
        for (const auto& [mid_monomial, innermost_poly] : inner_poly.coefficients) {
            for (const auto& [inner_monomial, coeff] : innermost_poly.coefficients) {
                // Create a 2D monomial from the two 1D monomials
                Monomial<2> combined_monomial;
                combined_monomial.exponents[0] = inner_monomial[0];
                combined_monomial.exponents[1] = mid_monomial[0];

                result.coefficients[outer_monomial][combined_monomial] += coeff;
            }
        }
    }

    return result;
}

/// A variant that can hold a `HashmapPolynomial`, `SparseOssifiedPolynomial`, or `DenseOssifiedPolynomial`.
/// It provides a unified interface for polynomial operations.
template <typename T, usize NVARS> struct Polynomial {
    std::variant<HashmapPolynomial<T, NVARS>, SparseOssifiedPolynomial<T, NVARS>, DenseOssifiedPolynomial<T, NVARS>>
        poly;

    Polynomial() : poly(HashmapPolynomial<T, NVARS>{}) {}
    Polynomial(HashmapPolynomial<T, NVARS> p) : poly(std::move(p)) {}
    Polynomial(SparseOssifiedPolynomial<T, NVARS> p) : poly(std::move(p)) {}
    Polynomial(DenseOssifiedPolynomial<T, NVARS> p) : poly(std::move(p)) {}

    HashmapPolynomial<T, NVARS> to_hashmap() const {
        return std::visit(
            [](const auto& p) -> HashmapPolynomial<T, NVARS> {
                using P = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<P, HashmapPolynomial<T, NVARS>>) {
                    return p;
                } else {
                    return p.to_hashmap();
                }
            },
            poly);
    }

    exp_t degree() const {
        return std::visit([](const auto& p) { return p.degree(); }, poly);
    }
    exp_t degree(usize variable) const {
        return std::visit([=](const auto& p) { return p.degree(variable); }, poly);
    }
    T eval(std::array<T, NVARS> values) const {
        if (auto *p = std::get_if<HashmapPolynomial<T, NVARS>>(&poly)) {
            return p->eval(values);
        }
        throw std::runtime_error("eval not supported for ossified polynomials, use eval_with_precalculated_powers");
    }
    T eval_with_precalculated_powers(std::array<std::span<const T>, NVARS> powers) const {
        return std::visit([&](const auto& p) { return p.eval_with_precalculated_powers(powers); }, poly);
    }
    T *get(const Monomial<NVARS>& monomial) {
        if (auto *p = std::get_if<HashmapPolynomial<T, NVARS>>(&poly)) {
            return p->get(monomial);
        }
        throw std::runtime_error("get not supported for ossified polynomials");
    }
    T *get_or_insert(const Monomial<NVARS>& monomial) {
        if (auto *p = std::get_if<HashmapPolynomial<T, NVARS>>(&poly)) {
            return p->get_or_insert(monomial);
        }
        throw std::runtime_error("get_or_insert not supported for ossified polynomials");
    }
    void set(const Monomial<NVARS>& monomial, const T& value) {
        if (auto *p = std::get_if<HashmapPolynomial<T, NVARS>>(&poly)) {
            p->set(monomial, value);
            return;
        }
        throw std::runtime_error("set not supported for ossified polynomials");
    }
    Polynomial<T, NVARS> partial_derivative(usize variable) const {
        return Polynomial(to_hashmap().partial_derivative(variable));
    }
    Polynomial<T, NVARS> operator+(const Polynomial<T, NVARS>& other) const {
        return Polynomial(to_hashmap() + other.to_hashmap());
    }
    Polynomial<T, NVARS> operator-(const Polynomial<T, NVARS>& other) const {
        return Polynomial(to_hashmap() - other.to_hashmap());
    }
    Polynomial<T, NVARS> operator*(const Polynomial<T, NVARS>& other) const {
        return Polynomial(to_hashmap() * other.to_hashmap());
    }
    Polynomial<T, NVARS> operator*(const T& scalar) const {
        return Polynomial(to_hashmap() * scalar);
    }
    Polynomial<T, NVARS> pow(exp_t exponent) const {
        return Polynomial(to_hashmap().pow(exponent));
    }
    Polynomial<T, NVARS> substitute(usize variable, Polynomial<T, NVARS> inner_val) const {
        return Polynomial(to_hashmap().substitute(variable, inner_val.to_hashmap()));
    }
    operator std::string() const {
        return std::visit([](const auto& p) { return std::string(p); }, poly);
    }
    void iterate(const std::function<void(const Monomial<NVARS>&, const T&)>& func) const {
        return std::visit([&](const auto& p) { p.iterate(func); }, poly);
    }
    template <typename U> Polynomial<U, NVARS> map(const std::function<U(const T&)>& func) const {
        return std::visit([&](const auto& p) { return Polynomial<U, NVARS>(p.map(func)); }, poly);
    }
};
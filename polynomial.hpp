#pragma once

#include "core.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <span>
#include <string>
#include <unordered_map>

#define PARENTHESIS_PRINT_MIN_EXPONENT 10

char get_var_name(usize i);

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
    std::size_t operator()(const Monomial<NVARS>& m) const {
        std::size_t h = 0;
        for (exp_t e : m.exponents) {
            h = h * 31 + std::hash<exp_t>{}(e);
        }
        return h;
    }
};

template <typename T, usize NVARS> struct Polynomial {
    std::unordered_map<Monomial<NVARS>, T> coefficients;
    usize max_degree() const {
        usize max_degree = 0;
        for (const auto& [monomial, _] : coefficients) {
            max_degree = std::max(max_degree, monomial.degree());
        }
        return max_degree;
    }
    usize max_degree(usize variable) const {
        usize max_degree = 0;
        for (const auto& [monomial, _] : coefficients) {
            max_degree = std::max(max_degree, monomial.exponents[variable]);
        }
        return max_degree;
    }
    T eval(std::span<const T> values) const {
        T result = 0;
        for (const auto& [monomial, coefficient] : coefficients) {
            result += coefficient * monomial.eval(values);
        }
        return result;
    }
    template <typename U> Polynomial<U, NVARS> map(const std::function<U(T)>& func) const {
        Polynomial<U, NVARS> result;
        for (const auto& [monomial, coefficient] : coefficients) {
            result.coefficients[monomial] = func(coefficient);
        }
        return result;
    }
    Polynomial<T, NVARS> operator+(const Polynomial<T, NVARS>& other) const {
        Polynomial<T, NVARS> result = *this;
        for (const auto& [monomial, coefficient] : other.coefficients) {
            result.coefficients[monomial] += coefficient;
        }
        return result;
    }
    Polynomial<T, NVARS> operator-(const Polynomial<T, NVARS>& other) const {
        Polynomial<T, NVARS> result = *this;
        for (const auto& [monomial, coefficient] : other.coefficients) {
            result.coefficients[monomial] -= coefficient;
        }
        return result;
    }
    Polynomial<T, NVARS> operator*(const Polynomial<T, NVARS>& other) const {
        Polynomial<T, NVARS> result;
        for (const auto& [monomial1, coefficient1] : coefficients) {
            for (const auto& [monomial2, coefficient2] : other.coefficients) {
                result.coefficients[monomial1 * monomial2] += coefficient1 * coefficient2;
            }
        }
        return result;
    }
    Polynomial<T, NVARS> operator*(const T& scalar) const {
        Polynomial<T, NVARS> result = *this;
        for (auto& [monomial, coefficient] : result.coefficients) {
            coefficient *= scalar;
        }
        return result;
    }
    Polynomial<T, NVARS> pow(exp_t exponent) const {
        assert(exponent >= 0 && "Exponent must be non-negative");
        if (exponent == 0) {
            Polynomial<T, NVARS> result;
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
    std::unordered_map<exp_t, Polynomial<T, NVARS>> get_powers(std::span<const exp_t> exponents) const {
        assert(std::is_sorted(exponents.begin(), exponents.end()) && "Exponents must be sorted in ascending order");
        std::unordered_map<exp_t, Polynomial<T, NVARS>> result;
        std::function<Polynomial<T, NVARS>(exp_t)> calculate_power;
        calculate_power = [&result, this, &calculate_power](exp_t e) {
            assert(e >= 0 && "Exponent must be non-negative");
            if (e == 0) {
                Polynomial<T, NVARS> power;
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
    constexpr Polynomial<T, NVARS> substitute(usize variable, Polynomial<T, NVARS> inner_val) const {
        Polynomial<T, NVARS> result;
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
    friend std::ostream& operator<<(std::ostream& os, const Polynomial<T, NVARS>& poly) {
        return os << std::string(poly);
    }
    T& operator[](const Monomial<NVARS>& monomial) {
        return coefficients[monomial];
    }
    const T& operator[](const Monomial<NVARS>& monomial) const {
        return coefficients.at(monomial);
    }
    template <typename... Args> T& set(Args... exponents) {
        static_assert(sizeof...(Args) == NVARS, "Number of exponents must match number of variables");
        return coefficients[Monomial<NVARS>({static_cast<exp_t>(exponents)...})];
    }
    template <typename... Args> const T& get(Args... exponents) const {
        static_assert(sizeof...(Args) == NVARS, "Number of exponents must match number of variables");
        return coefficients.at(Monomial<NVARS>({static_cast<exp_t>(exponents)...}));
    }
    Polynomial<Polynomial<T, NVARS - 1>, 1> unnest_outer(usize outer_variable) const {
        Polynomial<Polynomial<T, NVARS - 1>, 1> result;
        for (const auto& [monomial, coefficient] : coefficients) {
            auto [outer, inner] = monomial.split(outer_variable);
            result.coefficients[outer].coefficients[inner] += coefficient;
        }
        return result;
    }
    Polynomial<Polynomial<T, 1>, NVARS - 1> unnest_inner(usize inner_variable) const {
        Polynomial<Polynomial<T, 1>, NVARS - 1> result;
        for (const auto& [monomial, coefficient] : coefficients) {
            auto [inner, outer] = monomial.split(inner_variable);
            result.coefficients[outer].coefficients[inner] += coefficient;
        }
        return result;
    }
};

template <typename T, usize NDIM>
Polynomial<Polynomial<T, 2>, NDIM> merge_coeffs(Polynomial<Polynomial<Polynomial<T, 1>, 1>, NDIM> nested_poly) {
    Polynomial<Polynomial<T, 2>, NDIM> result;

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

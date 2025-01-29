#include <vector>
#include <iostream>
#include <unordered_map>
#include <array>
#include <span>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <functional>

typedef size_t usize;
typedef short int exp_t;

const char *VAR_NAMES[4] = {"x", "y", "z", "w"};
const char *get_var_name(usize i) {
    assert(i < 4 && "Index out of bounds");
    return VAR_NAMES[i];
}

template <usize NVARS>
struct Monomial {
    std::array<exp_t, NVARS> exponents;
    constexpr Monomial() : exponents() {}
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
    template <typename T>
    constexpr T eval(std::span<const T> values) const {
        T result = 1;
        for (usize i = 0; i < NVARS; ++i) {
            result *= std::pow(values[i], exponents[i]);
        }
        return result;
    }
    operator std::string() const {
        std::string result;
        for (usize i = 0; i < NVARS; ++i) {
            if (exponents[i] != 0) {
                result += std::string(get_var_name(i));
                if (exponents[i] > 1) {
                    result += "^" + std::to_string(exponents[i]);
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

template<usize NVARS>
struct std::hash<Monomial<NVARS>> {
    std::size_t operator()(const Monomial<NVARS>& m) const {
        std::size_t h = 0;
        for (exp_t e : m.exponents) {
            h = h * 31 + std::hash<exp_t>{}(e);
        }
        return h;
    }
};

template <typename T, usize NVARS>
struct Polynomial {
    std::unordered_map<Monomial<NVARS>, T> coefficients;
    constexpr T eval(std::span<const T> values) const {
        T result = 0;
        for (const auto& [monomial, coefficient] : coefficients) {
            result += coefficient * monomial.eval(values);
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
    Polynomial<T, NVARS> operator*(const Polynomial<T, NVARS>& other) const {
        Polynomial<T, NVARS> result;
        for (const auto& [monomial1, coefficient1] : coefficients) {
            for (const auto& [monomial2, coefficient2] : other.coefficients) {
                result.coefficients[monomial1 * monomial2] += coefficient1 * coefficient2;
            }
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
            if (coefficient == std::floor(coefficient)) {
                result += std::to_string((int)coefficient);
                if (monomial != Monomial<NVARS>()) {
                    result += " * " + std::string(monomial);
                }
            } else {
                result += std::to_string(coefficient);
                if (monomial != Monomial<NVARS>()) {
                    result += " * " + std::string(monomial);
                }
            }
            ++it;
            if (it != end) {
                result += " + ";
            }
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
    template<typename... Args>
    T& set(Args... exponents) {
        static_assert(sizeof...(Args) == NVARS, "Number of exponents must match number of variables");
        return coefficients[Monomial<NVARS>({static_cast<exp_t>(exponents)...})];
    }
    template<typename... Args>
    const T& get(Args... exponents) const {
        static_assert(sizeof...(Args) == NVARS, "Number of exponents must match number of variables");
        return coefficients.at(Monomial<NVARS>({static_cast<exp_t>(exponents)...}));
    }
};

int main() {
    Polynomial<double, 2> p, p1, p2;
    p.set(2, 0) = 7;
    p.set(0, 2) = 11;
    p1.set(1, 0) = 2;
    p1.set(0, 0) = 3;
    p2.set(0, 1) = 6;
    p2.set(0, 0) = 5;
    std::cout << "(" << p << ") o_x (" << p1 << ") = (" << p.substitute(0, p1) << ")" << std::endl;
    std::cout << "(" << p << ") o_y (" << p1 << ") = (" << p.substitute(1, p1) << ")" << std::endl;
    std::cout << "(" << p << ") o_x (" << p2 << ") = (" << p.substitute(0, p2) << ")" << std::endl;
    std::cout << "(" << p << ") o_y (" << p2 << ") = (" << p.substitute(1, p2) << ")" << std::endl;
    return 0;
}



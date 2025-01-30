#include <vector>
#include <iostream>
#include <unordered_map>
#include <array>
#include <span>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <functional>

#include "polynomial.hpp"
#include "input.hpp"

static std::vector<char> var_names = {'x', 'y', 'z', 'w'};
char get_var_name(usize i) {
    if (i < var_names.size()) {
        return var_names[i];
    }
    throw std::runtime_error("Variable index out of bounds");
}
usize get_var_index(char c) {
    for (usize i = 0; i < var_names.size(); ++i) {
        if (var_names[i] == c) {
            return i;
        }
    }
    throw std::runtime_error("Variable not found");
}

template<usize NVARS>
Polynomial<double, NVARS> parse_ast(ASTNodePtrType node) {
    auto result = Polynomial<double, NVARS>();
    if (node->type == ASTNodeType::Const) {
        result[Monomial<NVARS>()] = std::get<ASTNodeConstType>(node->value);
        return result;
    }
    if (node->type == ASTNodeType::Var) {
        auto monomial = Monomial<NVARS>();
        monomial.exponents[get_var_index(std::get<ASTNodeVarType>(node->value))] = 1;
        result[monomial] = 1;
        return result;
    }
    auto&& [left, right] = std::get<ASTNodeForkType>(node->value);
    auto left_node = parse_ast<NVARS>(std::move(left));
    auto right_node = parse_ast<NVARS>(std::move(right));
    if (node->type == ASTNodeType::Add) {
        return left_node + right_node;
    }
    if (node->type == ASTNodeType::Sub) {
        return left_node - right_node;
    }
    if (node->type == ASTNodeType::Mul) {
        return left_node * right_node;
    }
    if (node->type == ASTNodeType::Pow) {
        auto base = std::move(left_node);
        auto exponent = std::move(right_node);
        auto exp = exponent.coefficients[Monomial<NVARS>()];
        for (const auto& [monomial, coefficient] : exponent.coefficients) {
            assert(monomial == Monomial<NVARS>() && "Exponent must be a single monomial");
        }
        assert((exp_t)exp == exp && "Exponent must be an integer");
        assert(exp >= 0 && "Exponent must be non-negative");
        return base.pow(exp);
    }
    throw std::runtime_error("Invalid AST node type");
}

template <usize NVARS>
Polynomial<double, NVARS> parse_expression(std::string_view expression) {
    auto tokens = tokenize(expression);
    auto root = parse_tokens(tokens);
    return parse_ast<NVARS>(std::move(root));
}

int main() {
    std::string expression = "x(x^2+y^2-2y*1)(x^2+y^2-2y*2)(x^2+y^2-2y*4)(x^2+y^2-2y*8)(x^2+y^2-2y*16)(x^2+y^2-2y*32)(x^2+y^2-2y*64)(x^2+y^2-2y*128)";
    auto p0 = parse_expression<2>(expression);
    std::cout << p0 << std::endl;
    return 0;
}



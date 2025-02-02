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
#include "image.hpp"
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
    // std::string expression = "(x^2+y^2-1)xy(x^2-y^2-1)(y^2-x^2-1)(x^2-y^2)";
    // std::string expression = "(y^2 + x^2 - 1)^3 - (x^2)*(y^3)";
    // std::string expression = "(x^2+y^2-1)^3-x^2y^3";
    // std::string expression = "(x^2+y^2)^5-(x^4-6x^2y^2+y^4)^2";
    // std::string expression = "2.8x^2(x^2(2.5x^2+y^2-2)+1.2y^2(y(3y-0.75)-6.0311)+3.09)+0.98y^2((y^2-3.01)y^2+3)-1.005";
    std::string expression1, expression2;
    std::cin >> expression1 >> expression2;
    auto p1 = parse_expression<2>(expression1);
    auto p2 = parse_expression<2>(expression2);
    auto num_steps = 100;
    auto width = 1280;
    auto height = 1280;
    auto [x0, x1, y0, y1] = std::array{-2.0, 2.0, -2.0, 2.0};
    for (usize i = 0; i < num_steps; ++i) {
        auto lambda = (double)i / (num_steps - 1);
        auto p = p1 * lambda + p2 * (1 - lambda);
        auto img = Image<GrayscalePixel>(width, height);
        auto to_plane = [&](usize px, usize py) {
            double x = x0 + (x1 - x0) * px / (width - 1);
            double y = y0 + (y1 - y0) * py / (height - 1);
            return std::array{x, y};
        };
        auto soft_clamp = [](double v) {
            return 1 - std::abs(std::atan(v * 100) / M_PI * 2);
        };
        for (usize y = 0; y < img.height; ++y) {
            for (usize x = 0; x < img.width; ++x) {
                auto plane_coords = to_plane(x, y);
                auto pval = p.eval(plane_coords);
                img(x, y).set_v(soft_clamp(pval));
            }
        }
        std::string filename_fwd = "test" + std::format("{:04}", i) + ".bmp";
        std::string filename_bwd = "test" + std::format("{:04}", num_steps * 2 - 1 - i) + ".bmp";
        img.save_bmp(filename_fwd);
        img.save_bmp(filename_bwd);
    }
    return 0;
}



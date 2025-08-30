#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <span>
#include <vector>

#include "image.hpp"
#include "input.hpp"
#include "polynomial.hpp"
#include "render.hpp"
#include "thread.hpp"

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

template <usize NVARS> Polynomial<double, NVARS> parse_ast(ASTNodePtrType node) {
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

template <usize NVARS> Polynomial<double, NVARS> parse_expression(std::string_view expression) {
    auto tokens = tokenize(expression);
    auto root = parse_tokens(tokens);
    return parse_ast<NVARS>(std::move(root));
}

Texture2D<Real> render_image(Polynomial<double, 2> p, ImageParams params) {
    auto img = Texture<2, Real>(params.width, params.height);
    for (usize y = 0; y < img.height; ++y) {
        for (usize x = 0; x < img.width; ++x) {
            auto plane_coords = params.to_plane(x, y);
            auto pval = p.eval(plane_coords);
            img(x, y).set_v(params.soft_clamp(pval));
        }
    }
    return img;
}

std::vector<Texture2D<Real>> render_images(AnimationParams& params) {
    return parallel_map<Texture2D<Real>>(
        [&](double lambda) {
            auto p = params.p1 * lambda + params.p2 * (1 - lambda);
            return render_image(p, params.image);
        },
        std::span(params.lambdas));
}

template <typename P> double image_difference(Texture2D<P>& left, Texture2D<P>& right) {
    assert(left.width == right.width && left.height == right.height);
    double sum = 0.0;
    for (usize y = 0; y < left.height; ++y) {
        for (usize x = 0; x < left.width; ++x) {
            auto& l = left(x, y);
            auto& r = right(x, y);
            double dr = l.r() - r.r();
            double dg = l.g() - r.g();
            double db = l.b() - r.b();
            sum += dr * dr + dg * dg + db * db;
        }
    }
    return sum / left.width / left.height;
}

int main() {
    auto p = parse_expression<4>("(x^2+y^2-1)xy(x^2-y^2-1)(y^2-x^2-1)(x^2-y^2)");
    auto xsub = parse_expression<4>("x+z");
    auto ysub = parse_expression<4>("y+w");
    p = p.substitute(0, xsub);
    p = p.substitute(1, ysub);
    std::cout << p << std::endl;
    return 1;
    // std::string expression = "(x^2+y^2-1)xy(x^2-y^2-1)(y^2-x^2-1)(x^2-y^2)";
    // std::string expression = "(y^2 + x^2 - 1)^3 - (x^2)*(y^3)";
    // std::string expression = "(x^2+y^2-1)^3-x^2y^3";
    // std::string expression = "(x^2+y^2)^5-(x^4-6x^2y^2+y^4)^2";
    // std::string expression =
    // "2.8x^2(x^2(2.5x^2+y^2-2)+1.2y^2(y(3y-0.75)-6.0311)+3.09)+0.98y^2((y^2-3.01)y^2+3)-1.005";
    std::string expression1, expression2;
    std::cin >> expression1 >> expression2;
    auto p1 = parse_expression<2>(expression1);
    auto p2 = parse_expression<2>(expression2);
    usize width = 1440;
    usize height = 1440;
    auto plane_height = 4.0;
    auto aspect_ratio = (double)width / height;
    auto [x0, x1, y0, y1] = std::array{
        -plane_height * aspect_ratio / 2, plane_height * aspect_ratio / 2, -plane_height / 2, plane_height / 2};
    auto to_plane = [&](usize px, usize py) {
        double x = x0 + (x1 - x0) * px / (width - 1);
        double y = y0 + (y1 - y0) * py / (height - 1);
        return std::array{x, y};
    };
    auto soft_clamp = [](double v) -> double { return 1 - std::abs(std::atan(v * 100) / M_PI * 2); };
    auto img_params = ImageParams{to_plane, soft_clamp, width, height};
    std::map<double, Texture2D<Real>> interpolation_steps = std::map<double, Texture2D<Real>>();
    interpolation_steps.insert({0.0, render_image(p2, img_params)});
    interpolation_steps.insert({1.0, render_image(p1, img_params)});
    std::vector<std::pair<double, double>> candidate_intervals;
    candidate_intervals.push_back({0.0, 1.0});
    const auto max_distance = 0.005;
    while (!candidate_intervals.empty()) {
        auto new_candidates = std::vector<std::pair<double, double>>();
        auto to_render = std::vector<double>();
        auto distances = parallel_map<double>(
            [&](auto interval) {
                auto [left, right] = interval;
                auto l_img_iter = interpolation_steps.find(left);
                auto r_img_iter = interpolation_steps.find(right);
                assert(l_img_iter != interpolation_steps.end() && r_img_iter != interpolation_steps.end());
                auto& l_img = l_img_iter->second;
                auto& r_img = r_img_iter->second;
                auto distance = image_difference(l_img, r_img);
                return distance;
            },
            std::span(candidate_intervals));
        for (usize i = 0; i < candidate_intervals.size(); ++i) {
            auto [left, right] = candidate_intervals[i];
            auto distance = distances[i];
            if (distance > max_distance && right - left > 1e-14) {
                auto midpoint = (left + right) / 2;
                to_render.push_back(midpoint);
                new_candidates.emplace_back(left, midpoint);
                new_candidates.emplace_back(midpoint, right);
            }
        }
        auto ani_params = AnimationParams{img_params, p1, p2, std::span(to_render)};
        auto new_images = render_images(ani_params);
        for (usize i = 0; i < to_render.size(); ++i) {
            std::pair<double, Texture2D<Real>> kv = std::pair(to_render[i], std::move(new_images[i]));
            interpolation_steps.insert(kv);
        }
        candidate_intervals = new_candidates;
    }
    auto num_images = interpolation_steps.size();
    auto num_end_reps = 4;
    std::vector<std::pair<std::string, Texture2D<Real>&>> work_queue;
    usize i = 0;
    for (auto& [t, img] : interpolation_steps) {
        auto num_reps = 1;
        if (i == 0 || i == num_images + num_end_reps - 2) num_reps = num_end_reps;
        for (usize j = 0; j < num_reps; ++j) {
            std::string filename_fwd = "zzz" + std::format("{:04}", i) + ".bmp";
            std::string filename_bwd = "zzz" + std::format("{:04}", num_images * 2 + 4 * num_end_reps - i - 1) + ".bmp";
            work_queue.push_back({filename_fwd, img});
            work_queue.push_back({filename_bwd, img});
            ++i;
        }
    }
    parallel_for(work_queue.size(), [&](usize i) {
        auto& [filename, img] = work_queue[i];
        img.save_bmp(filename);
    });
    return 0;
}

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <span>
#include <unordered_set>
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

template <usize NVARS> HashmapPolynomial<double, NVARS> parse_ast(ASTNodePtrType node) {
    auto result = HashmapPolynomial<double, NVARS>();
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

template <usize NVARS> HashmapPolynomial<double, NVARS> parse_expression(std::string_view expression) {
    auto tokens = tokenize(expression);
    auto root = parse_tokens(tokens);
    return parse_ast<NVARS>(std::move(root));
}

double image_difference(const SparseTexture<2, BlackWhite>& left, const SparseTexture<2, BlackWhite>& right) {
    assert(left.dims == right.dims);

    usize left_not_right = 0;
    usize right_not_left = 0;

    std::unordered_set<usize> keys;
    for (const auto& [idx, _] : left.data) {
        keys.insert(idx);
        if (right.data.find(idx) == right.data.end()) {
            left_not_right++;
        }
    }
    for (const auto& [idx, _] : right.data) {
        keys.insert(idx);
        if (left.data.find(idx) == left.data.end()) {
            right_not_left++;
        }
    }
    if (keys.empty()) {
        return 0.0;
    }
    return (left_not_right + right_not_left) / (double)(keys.size());
}

int main() {
    auto max_distance = 0.4;
    usize width = 2048;
    auto plane_height = 4.0;
    std::string intermediate_dir = "intermediate_images";

    try {
        std::filesystem::create_directories(intermediate_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cout << "Warning: Could not create intermediate directory. Images will be saved in current directory."
                  << std::endl;
        intermediate_dir = ".";
    }
    std::string expression1, expression2;
    std::cin >> expression1 >> expression2;
    auto p1 = parse_expression<2>(expression1);
    auto p2 = parse_expression<2>(expression2);
    auto img_params = ImageParams{width, width};

    usize degree = std::max(p1.degree(), p2.degree());
    usize max_granularity = (usize)std::log2((double)width);
    PreparedLattices lattices(-plane_height / 2, plane_height / 2, degree, max_granularity);

    std::cout << "Using " << num_threads() << " threads" << std::endl;
    std::map<double, SparseTexture<2, BlackWhite>> interpolation_steps;
    auto zeroone = std::vector<double>({0.0, 1.0});
    auto ani_params0 = AnimationParams{img_params, p1, p2, std::span<double>(zeroone)};
    auto images = render_images(lattices, ani_params0);
    interpolation_steps.insert({0.0, std::move(images[0])});
    interpolation_steps.insert({1.0, std::move(images[1])});
    std::vector<std::pair<double, double>> candidate_intervals;
    candidate_intervals.push_back({0.0, 1.0});
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
        auto max_current_distance = 0.0;
        for (usize i = 0; i < candidate_intervals.size(); ++i) {
            auto [left, right] = candidate_intervals[i];
            auto distance = distances[i];
            if (distance > max_current_distance) max_current_distance = distance;
            if (distance > max_distance && right - left > 4 * std::numeric_limits<double>::epsilon()) {
                auto midpoint = (left + right) / 2;
                to_render.push_back(midpoint);
                new_candidates.emplace_back(left, midpoint);
                new_candidates.emplace_back(midpoint, right);
            }
        }
        std::cout << "Max distance: " << max_current_distance << std::endl;
        auto ani_params = AnimationParams{img_params, p1, p2, std::span(to_render)};
        auto new_images = render_images(lattices, ani_params);
        for (usize i = 0; i < to_render.size(); ++i) {
            std::pair<double, SparseTexture<2, BlackWhite>> kv = std::pair(to_render[i], std::move(new_images[i]));
            interpolation_steps.insert(kv);
        }
        candidate_intervals = new_candidates;
    }
    auto num_images = interpolation_steps.size();
    auto num_end_reps = 12;
    usize forward_frame_count = (num_images - 2) + 2 * num_end_reps;
    std::vector<std::pair<std::string, double>> work_queue;
    usize i = 0;
    for (auto const& [t, img] : interpolation_steps) {
        auto num_reps = 1;
        if (i == 0 || i == num_images + num_end_reps - 2) num_reps = num_end_reps;
        for (usize j = 0; j < num_reps; ++j) {
            std::string filename_fwd = intermediate_dir + "/zzz" + std::format("{:04}", i) + ".bmp";
            std::string filename_bwd =
                intermediate_dir + "/zzz" + std::format("{:04}", 2 * forward_frame_count - i - 1) + ".bmp";
            work_queue.push_back({filename_fwd, t});
            work_queue.push_back({filename_bwd, t});
            ++i;
        }
    }
    std::cout << "Saving " << work_queue.size() << " images" << std::endl;
    parallel_for(work_queue.size(), [&](usize i) {
        auto& [filename, t] = work_queue[i];
        interpolation_steps.at(t).to_dense().save_bmp(filename);
    });

    // Run ffmpeg to create video from images
    std::cout << "Running ffmpeg to create video..." << std::endl;

    // Check if ffmpeg is available
#ifdef _WIN32
    int ffmpeg_check = std::system("ffmpeg -version >nul 2>&1");
#else
    int ffmpeg_check = std::system("ffmpeg -version >/dev/null 2>&1");
#endif
    if (ffmpeg_check != 0) {
        std::cout << "Warning: ffmpeg not found in PATH. Please install ffmpeg to create video." << std::endl;
        std::cout << "Images have been saved in the '" << intermediate_dir << "' directory." << std::endl;
        return 0;
    }

    std::string ffmpeg_cmd;
#ifdef _WIN32
    ffmpeg_cmd = "ffmpeg -y -framerate 30 -i \"" + intermediate_dir +
        "\\zzz%04d.bmp\" -c:v libx264 -pix_fmt yuv420p -crf 23 output_video.mp4";
#else
    ffmpeg_cmd = "ffmpeg -y -framerate 30 -i " + intermediate_dir +
        "/zzz%04d.bmp -c:v libx264 -pix_fmt yuv420p -crf 23 output_video.mp4";
#endif

    std::cout << "Executing: " << ffmpeg_cmd << std::endl;
    int result = std::system(ffmpeg_cmd.c_str());
    if (result == 0) {
        std::cout << "Video created successfully: output_video.mp4" << std::endl;
    } else {
        std::cout << "Error creating video with ffmpeg. Exit code: " << result << std::endl;
        std::cout << "Images have been saved in the '" << intermediate_dir << "' directory." << std::endl;
    }
    if (intermediate_dir != ".") {
        std::filesystem::remove_all(intermediate_dir);
    } else {
        for (usize i = 0; i < work_queue.size(); ++i) {
            std::filesystem::remove(work_queue[i].first);
        }
    }

    return 0;
}

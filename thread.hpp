#pragma once
#include <functional>
#include <thread>
#include <vector>
#include "core.hpp"
#include <span>

usize num_threads();

std::thread spawn_thread(const std::function<void()>& func, usize _allowed_threads);

/// Run the function `func` for each index from 0 to `n - 1`.
/// If there are threads available (i.e. `num_threads() > 1`), they will be used to run the function in parallel.
/// Each thread will be assigned at least `min_work_per_thread` indices, with the exception of the last thread.
void parallel_for(usize n, const std::function<void(usize)>& func, usize min_work_per_thread = 1);

/// Apply the function `func` to each element of `args` in parallel.
/// The results are collected into a vector.
template<typename RETVAL, typename ARG, typename T>
std::vector<RETVAL> parallel_map(const T& func, std::span<ARG> args) {
    std::vector<RETVAL> result;
    result.resize(args.size());
    std::function<RETVAL(ARG)> f = [func](ARG arg) -> RETVAL {
        return func(arg);
    };
    parallel_for(args.size(), [&](usize i) {
        result[i] = f(args[i]);
    });
    return result;
}

/// Apply the function `func` to each index from 0 to `n - 1` in parallel.
/// The results are collected into a vector.
template<typename RETVAL, typename T>
std::vector<RETVAL> parallel_nat_map(usize n, const T& func) {
    std::vector<RETVAL> results;
    results.resize(n);
    std::function<RETVAL(usize)> f = [func](usize i) -> RETVAL {
        return func(i);
    };
    parallel_for(n, [&](usize i) {
        results[i] = f(i);
    });
    return results;
}
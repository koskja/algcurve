#pragma once
#include <functional>
#include <thread>
#include <vector>
#include "core.hpp"
#include <iostream>

thread_local usize allowed_threads = std::thread::hardware_concurrency();

usize num_threads() {
    return allowed_threads;
}

std::thread spawn_thread(const std::function<void()>& func, usize _allowed_threads) {
    std::thread thread([func, _allowed_threads]() {
        allowed_threads = _allowed_threads;
        func();
    });
    return thread;
}

void parallel_for(usize n, const std::function<void(usize)>& func, usize min_work_per_thread = 1) {
    if (n == 0) return;
    if (allowed_threads == 1) {
        for (usize i = 0; i < n; ++i) {
            func(i);
        }
        return;
    }
    usize rem_threads = n / min_work_per_thread;
    if (allowed_threads < rem_threads) rem_threads = allowed_threads;
    std::vector<std::thread> threads;
    usize i = 0;
    while (i < n) {
        usize to_take;
        if (rem_threads > 1) {
            to_take = (n - i + rem_threads - 1) / rem_threads;
        } else {
            to_take = n - i;
        }
        threads.push_back(spawn_thread([i, to_take, &func]() {
            for (usize j = i; j < i + to_take; ++j) {
                func(j);
            }
        }, 1));
        i += to_take;
        rem_threads -= 1;
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

template<typename RETVAL, typename ARG, typename T>
std::vector<RETVAL> parallel_map(T func, std::span<ARG> args) {
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

template<typename RETVAL, typename T>
std::vector<RETVAL> parallel_nat_map(usize n, T func) {
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
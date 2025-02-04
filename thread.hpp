#pragma once
#include <functional>
#include <thread>
#include <vector>
#include "core.hpp"
#include <iostream>

usize num_threads() {
    return std::thread::hardware_concurrency();
    //return 1;
}

void parallel_for(usize n, std::function<void(usize)> func) {
    usize rem_threads = num_threads();
    std::vector<std::thread> threads;
    auto tfunc = [&](usize start, usize end) {
        for (usize i = start; i < end; ++i) {
            func(i);
        }
    };
    usize i = 0;
    while (i < n) {
        usize to_take;
        if (rem_threads > 1) {
            to_take = (n - i + rem_threads - 1) / rem_threads;
        } else {
            to_take = n - i;
        }
        threads.push_back(std::thread(tfunc, i, i + to_take));
        i += to_take;
        rem_threads -= 1;
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

template<typename RETVAL, typename ARG>
std::vector<RETVAL> parallel_map(std::function<RETVAL(ARG)> func, std::span<ARG> args) {
    std::vector<RETVAL> result;
    result.resize(args.size());
    parallel_for(args.size(), [&](usize i) {
        result[i] = func(args[i]);
    });
    return result;
}

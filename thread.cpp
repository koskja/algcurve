#include "thread.hpp"

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
void parallel_for(usize n, const std::function<void(usize)>& func, usize min_work_per_thread) {
    if (n == 0) return;
    if (allowed_threads == 1) {
        for (usize i = 0; i < n; ++i) {
            func(i);
        }
        return;
    }
    usize rem_threads = n / min_work_per_thread;
    if (allowed_threads < rem_threads) rem_threads = allowed_threads;
    usize inherited_threads = allowed_threads - rem_threads;
    std::vector<std::thread> threads;
    usize i = 0;
    while (i < n) {
        usize to_take;
        usize to_inherit;
        if (rem_threads > 1) {
            to_take = (n - i + rem_threads - 1) / rem_threads;
            to_inherit = inherited_threads / rem_threads;
        } else {
            to_take = n - i;
            to_inherit = inherited_threads;
        }
        threads.push_back(spawn_thread(
            [i, to_take, &func]() {
                for (usize j = i; j < i + to_take; ++j) {
                    func(j);
                }
            },
            1 + to_inherit));
        i += to_take;
        inherited_threads -= to_inherit;
        rem_threads -= 1;
    }
    for (auto& thread : threads) {
        thread.join();
    }
}
#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <tuple>
#include <chrono>
#include <vector>
#include <numeric>
#include <utility>

struct MyTime {
    MyTime() : time(std::chrono::steady_clock::now()) {}
    long long operator-(const MyTime& other) {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(time - other.time).count();
    }
private:
    std::chrono::steady_clock::time_point time;
};
MyTime getTime() {
    return MyTime();
}

// run once, return pair { functionReturnValue, time }
template<typename ReturnType, typename Function, typename... Args>
std::pair<ReturnType, long long> benchmark(Function function, Args&&... args) {
    auto start = getTime();
    ReturnType result = function(std::forward<Args>(args)...);
    auto end = getTime();
    return std::pair(result, end - start);
}

// run once, return time
template<typename Function, typename... Args> long long benchmark(Function function, Args&&... args) {
    auto start = getTime();
    function(std::forward<Args>(args)...);
    auto end = getTime();
    return end - start;
}


namespace details {
    template<typename Function, typename... Args> std::vector<long long> runMultipleTimes(int repetitions, Function function, Args&&... args) {
        repetitions += 1;
        std::vector<long long> times(repetitions);
        for (int i = repetitions-1; i >= 0; --i) {
            times[i] = benchmark(function, args...);
        }
        times.pop_back();
        return times;
    }
}

// run "repetitions" times, return avarage time
template<typename Function, typename... Args> double benchmark(int repetitions, Function function, Args&&... args) {
    auto times = details::runMultipleTimes(repetitions, function, std::forward<Args>(args)...);
    auto value = std::accumulate(times.begin(), times.end(), 0.0);
    return static_cast<double>(value) / repetitions;
}

// run "repetitions" times, return pair { functionReturnValue, avarage time }
template<typename ReturnType, typename Function, typename... Args>
std::pair<ReturnType, double> benchmark(int repetitions, Function function, Args&&... args) {
    return { function(args...), benchmark(repetitions, function, args...) };
}
#endif

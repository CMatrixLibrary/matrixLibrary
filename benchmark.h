#pragma once
#include <tuple>
#include <chrono>
#include <vector>
#include <numeric>

struct MyTime {
    MyTime() : time(std::chrono::steady_clock::now()) {}
    double operator-(const MyTime& other) {
        return std::chrono::duration<double>(time - other.time).count();
    }
private:
    std::chrono::steady_clock::time_point time;
};
MyTime getTime() {
    return MyTime();
}

// run once, return pair { functionReturnValue, time }
template<typename ReturnType, typename Function, typename... Args>
std::pair<ReturnType, double> benchmark(Function function, Args... args) {
    auto start = getTime();
    ReturnType result = function(args...);
    auto end = getTime();
    return std::pair(result, end - start);
}

// run once, return time
template<typename Function, typename... Args> double benchmark(Function function, Args... args) {
    auto start = getTime();
    function(args...);
    auto end = getTime();
    return end - start;
}


namespace details {
    template<typename Function, typename... Args> std::vector<double> runMultipleTimes(int repetitions, Function function, Args... args) {
        std::vector<double> times(repetitions);
        for (int i = 0; i < repetitions; ++i) {
            times[i] = benchmark(function, args...);
        }
        return times;
    }
}

// run "repetitions" times, return avarage time
template<typename Function, typename... Args> double benchmark(int repetitions, Function function, Args... args) {
    auto times = details::runMultipleTimes(repetitions, function, args...);
    auto value = std::accumulate(times.begin(), times.end(), 0.0);
    return value / repetitions;
}

// run "repetitions" times, return pair { functionReturnValue, avarage time }
template<typename ReturnType, typename Function, typename... Args>
std::pair<ReturnType, double> benchmark(int repetitions, Function function, Args... args) {
    return { function(args...), benchmark(repetitions, function, args...) };
}
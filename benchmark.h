#pragma once
#include <tuple>
#include <chrono>

struct MyTime {
    MyTime() : time(std::chrono::high_resolution_clock::now()) {}
    double operator-(const MyTime& other) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(time - other.time).count() / 1000.0;
    }
private:
    std::chrono::high_resolution_clock::time_point time;
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

// run "repetitions" times, return avarage time
template<typename Function, typename... Args> double benchmark(int repetitions, Function function, Args... args) {
    double time = 0;
    for (int i = 0; i < repetitions; ++i) {
        time += benchmark(function, args...);
    }
    return time / repetitions;
}

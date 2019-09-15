#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <tuple>
#include <chrono>
#include <vector>
#include <numeric>
#include <utility>

struct TimeResult {
    double value;

    TimeResult(double value=0) : value(value) {}
    operator double() { return value; }

    double getSec()   { return value; }
    double getMilli() { return value * 1'000; }
    double getMicro() { return value * 1'000'000; }
    double getNano()  { return value * 1'000'000'000; }
};

struct MyTime {
    MyTime() : time(std::chrono::steady_clock::now()) {}
    TimeResult operator-(const MyTime& other) {
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
std::pair<ReturnType, TimeResult> benchmark(Function function, Args&&... args) {
    auto start = getTime();
    ReturnType result = function(std::forward<Args>(args)...);
    auto end = getTime();
    return std::pair(result, end - start);
}

// run once, return time
template<typename Function, typename... Args> TimeResult benchmark(Function function, Args&&... args) {
    auto start = getTime();
    function(std::forward<Args>(args)...);
    auto end = getTime();
    return end - start;
}


namespace details {
    template<typename Function, typename... Args> std::vector<TimeResult> runMultipleTimes(int repetitions, Function function, Args&&... args) {
        repetitions += 1;
        std::vector<TimeResult> times(repetitions);
        for (int i = repetitions-1; i >= 0; --i) {
            times[i] = benchmark(function, args...);
        }
        times.pop_back();
        return times;
    }
}

// run "repetitions" times, return avarage time
template<typename Function, typename... Args> TimeResult benchmark(int repetitions, Function function, Args&&... args) {
    auto times = details::runMultipleTimes(repetitions, function, std::forward<Args>(args)...);
    auto value = std::accumulate(times.begin(), times.end(), 0.0);
    return value / repetitions;
}

// run "repetitions" times, return median time
template<typename Function, typename... Args> TimeResult benchmarkMedian(int repetitions, Function function, Args&&... args) {
    auto times = details::runMultipleTimes(repetitions, function, std::forward<Args>(args)...);
    if (repetitions % 2 == 0) {
        return (times[repetitions / 2] + times[repetitions / 2 - 1]) / 2;
    } else {
        return times[repetitions / 2];
    }
}

// run "repetitions" times, return pair { functionReturnValue, avarage time }
template<typename ReturnType, typename Function, typename... Args>
std::pair<ReturnType, TimeResult> benchmark(int repetitions, Function function, Args&&... args) {
    return { function(args...), benchmark(repetitions, function, args...) };
}
#endif

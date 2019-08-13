#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "Matrix.h"
#include "MatrixView.h"
#include "naiveBasicOperations.h"
#include "matrixOperators.h"
#include "Range.h"
#include "RangeZip.h"
#include "MatrixExtendedFunctions.h"
#include "benchmark.h"
#include "strassen.h"

template<typename T, typename Function> void strassenDynamicBenchmark(std::string label, int n, int steps, Function fun) {
    HeapMatrix<T> a(n, n);
    HeapMatrix<T> b(n, n);
    for (int i = 0; i < a.size(); ++i) {
        a.data()[i] = i;
        b.data()[i] = i * 2;
    }
    auto start = MyTime();
    auto c = fun(a, b, steps);
    auto end = MyTime();
    std::cout << "Dynamic " << label << ' ' << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
}
template<typename T, int N, typename Function> void strassenStaticBenchmark(std::string label, Function fun) {
    StaticHeapMatrix<T, N, N> a;
    StaticHeapMatrix<T, N, N> b;
    for (int i = 0; i < a.size(); ++i) {
        a.data()[i] = i;
        b.data()[i] = i * 2;
    }
    auto start = MyTime();
    auto c = fun(a, b);
    auto end = MyTime();
    std::cout << "Static  " << label << ' ' << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
}
template<typename T, int N, typename Function> void strassenStaticBenchmark(std::string label, int steps, Function fun) {
    StaticHeapMatrix<T, N, N> a;
    StaticHeapMatrix<T, N, N> b;
    for (int i = 0; i < a.size(); ++i) {
        a.data()[i] = i;
        b.data()[i] = i * 2;
    }
    auto start = MyTime();
    auto c = fun(a, b, steps);
    auto end = MyTime();
    std::cout << "Static  " << label << ' ' << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
}
template<typename T, int N, int Steps, bool OnlyAvx=false> void strassenTest(int n) {
    using StaticMatrix = StaticHeapMatrix<T, N, N>;
    using DynamicMatrix = HeapMatrix<T>;
    
    if (n != N) {
        std::cout << "incorrect stdin\n";
        return;
    }
    std::cout << std::left;
    std::cout << "=================================================================\n";
    std::cout << "N               : " << N << '\n';
    std::cout << "Recursive Steps : " << Steps << '\n';
    std::cout << "-------------------------------\n";
    if (!OnlyAvx) {
        strassenDynamicBenchmark<T>("High-level noAVX :", n, Steps, [](const auto& a, const auto& b, int steps) {
            return strassen(a, b, steps);
        });
        strassenStaticBenchmark<T, N>("High-level noAVX :", [](const auto& a, const auto& b) {
            return strassen<Steps>(a, b);
        });
        strassenDynamicBenchmark<T>("Low-level  noAVX :", n, Steps, [](const auto& a, const auto& b, int steps) {
            return lowLevelStrassen(a, b, steps);
        });
        strassenStaticBenchmark<T, N>("Low-level  noAVX :", Steps, [](const auto& a, const auto& b, int steps) {
            return lowLevelStrassen(a, b, steps);
        });
        strassenDynamicBenchmark<T>("Min-Space  noAVX :", n, Steps, [](const auto& a, const auto& b, int steps) {
            return minSpaceStrassen(a, b, steps);
        });
        strassenStaticBenchmark<T, N>("Min-Space  noAVX :", [](const auto& a, const auto& b) {
            return minSpaceStrassen<Steps>(a, b);
        });
    }
    strassenDynamicBenchmark<T>("High-level AVX   :", n, Steps, [](const auto& a, const auto& b, int steps) {
        return strassenAvx(a, b, steps);
    });
    strassenStaticBenchmark<T, N>("High-level AVX   :", [](const auto& a, const auto& b) {
        return strassenAvx<Steps>(a, b);
    });
    strassenDynamicBenchmark<T>("Low-level  AVX   :", n, Steps, [](const auto& a, const auto& b, int steps) {
        return lowLevelAvxStrassen(a, b, steps);
    });
    strassenStaticBenchmark<T, N>("Low-level  AVX   :", Steps, [](const auto& a, const auto& b, int steps) {
        return lowLevelAvxStrassen(a, b, steps);
    });
    strassenDynamicBenchmark<T>("Min-Space  AVX   :", n, Steps, [](const auto& a, const auto& b, int steps) {
        return minSpaceAvxStrassen(a, b, steps);
    });
    strassenStaticBenchmark<T, N>("Min-Space  AVX   :", [](const auto& a, const auto& b) {
        return minSpaceAvxStrassen<Steps>(a, b);
    });
    std::cout << "=================================================================\n\n";
}

template<typename T> void strassenTestForType() {
    std::array<int, 6> vec;
    std::cout << "copy line below and click enter: \n";
    std::cout << "480 992 1024 1920 2048 3840\n";
    std::cin >> vec[0] >> vec[1] >> vec[2] >> vec[3] >> vec[4] >> vec[5];
    strassenTest<T, 480, 0>(vec[0]);
    strassenTest<T, 480, 1>(vec[0]);
    strassenTest<T, 480, 2>(vec[0]);
    strassenTest<T, 992, 0>(vec[1]);
    strassenTest<T, 992, 1>(vec[1]);
    strassenTest<T, 992, 2>(vec[1]);
    strassenTest<T, 992, 3>(vec[1]);
    strassenTest<T, 992, 4>(vec[1]);
    strassenTest<T, 992, 5>(vec[1]);
    strassenTest<T, 1024, 0>(vec[2]);
    strassenTest<T, 1024, 1>(vec[2]);
    strassenTest<T, 1024, 2>(vec[2]);
    strassenTest<T, 1920, 0>(vec[3]);
    strassenTest<T, 1920, 1>(vec[3]);
    strassenTest<T, 1920, 2>(vec[3]);
    strassenTest<T, 1920, 3>(vec[3]);
    strassenTest<T, 1920, 4>(vec[3]);
    strassenTest<T, 1920, 5>(vec[3]);
    strassenTest<T, 2048, 2>(vec[4]);
    strassenTest<T, 2048, 3>(vec[4]);
    strassenTest<T, 3840, 0, true>(vec[5]);
    strassenTest<T, 3840, 1, true>(vec[5]);
    strassenTest<T, 3840, 2, true>(vec[5]);
}

int main() {
    strassenTestForType<float>();
    std::cin.get();
    return 0;
}

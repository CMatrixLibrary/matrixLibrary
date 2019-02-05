#include <vector>
#include <iostream>
#include "benchmark.h"
#include "FullMatrix.h"
#include "FullSubMatrix.h"
#include "FullMatrixView.h"
#include "strassenMultiply.h"
#include "fastMultiply3x3.h"

void strassenVsNaiveMulSpeedTest() {
    constexpr int N = 1024;

    FullMatrix<int> a(N, N);
    FullMatrix<int> b(N, N);
    for (int i = 0; i < N*N; ++i) {
        a[i] = i;
        b[i] = i * i;
    }

    auto[c, strassenTime] = benchmark<FullMatrix<int>>([&]() -> auto {
        return strassenMul(a, b);
    });
    auto[d, naiveTime] = benchmark<FullMatrix<int>>([&]() -> auto {
        return naiveMul(a, b);
    });

    std::cout << "strassenTime = " << strassenTime << '\n';
    std::cout << "naiveTime    = " << naiveTime << '\n';
}

void fast3x3VsNaiveMulSpeedTest() {
    constexpr int N = 729;

    FullMatrix<int> a(N, N);
    FullMatrix<int> b(N, N);
    for (int i = 0; i < N*N; ++i) {
        a[i] = i;
        b[i] = i * i;
    }

    auto[c, fast3x3Time] = benchmark<FullMatrix<int>>([&]() -> auto {
        return fastMul3x3(a, b);
    });
    auto[d, naiveTime] = benchmark<FullMatrix<int>>([&]() -> auto {
        return naiveMul(a, b);
    });

    std::cout << "fast3x3Time  = " << fast3x3Time << '\n';
    std::cout << "naiveTime    = " << naiveTime << '\n';
}


int main() {
    strassenVsNaiveMulSpeedTest();
    fast3x3VsNaiveMulSpeedTest();

    std::cin.get();
    return 0;
}
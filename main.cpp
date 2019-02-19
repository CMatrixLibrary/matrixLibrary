#include <vector>
#include <iostream>
#include "benchmark.h"
#include "FullMatrix.h"
#include "FullMatrixView.h"
#include "FullMatrixConstView.h"
#include "strassenMultiply.h"
#include "fastMultiply3x3.h"
#include "FastMatrixMultiplyGenerator.h"
#include "fastMatrixMultiplyAlgorithms/genStrassen.h"

void strassenVsNaiveMulSpeedTest() {
    constexpr int N = 1024;

    FullMatrix<int> a(N, N);
    FullMatrix<int> b(N, N);
    for (int i = 0; i < N*N; ++i) {
        a.data()[i] = i;
        b.data()[i] = i * i;
    }

    auto[c, strassenTime] = benchmark<FullMatrix<int>>([&]() -> auto {
        return strassenMul(a, b);
    });
    auto[d, naiveTime] = benchmark<FullMatrix<int>>([&]() -> auto {
        return naiveMul(a, b);
    });
    auto[e, genStrassenTime] = benchmark<FullMatrix<int>>([&]() -> auto {
        return genStrassen(a, b);
    });

    std::cout << "strassenTime    = " << strassenTime << '\n';
    std::cout << "genStrassenTime = " << genStrassenTime << '\n';
    std::cout << "naiveTime       = " << naiveTime << '\n';
}

void fast3x3VsNaiveMulSpeedTest() {
    constexpr int N = 729;

    FullMatrix<int> a(N, N);
    FullMatrix<int> b(N, N);
    for (int i = 0; i < N*N; ++i) {
        a.data()[i] = i;
        b.data()[i] = i * i;
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

/*
    all iteration methods for T=int compile the most inner loop in x86-64 architecture to:

    when using gcc 8.2:
    Loop:
        mov     DWORD PTR [rax], 10
        add     rax, 4
        cmp     rdx, rax
        jne     Loop

    when using clang 7.0.0:
    Loop:
        movups  xmmword ptr [rsi + 4*rdx], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 16], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 32], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 48], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 64], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 80], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 96], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 112], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 128], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 144], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 160], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 176], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 192], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 208], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 224], xmm0
        movups  xmmword ptr [rsi + 4*rdx + 240], xmm0
        add     rdx, 64
        add     rbp, 8
        jne     Loop

    In case of clang matrixIndexLoop was compiled to less efficient code (lower xmm usage by half).

    In practice all methods seem to be equal in speed using msvc 19.16 compiler.
*/

template<typename T> void matrixDirectIndexing(T* data, int n, int m) {
    // note that this method is not usable for FullMatrixView and FullMatrixConstView classes
    int size = n * m;
    for (int i = 0; i < size; ++i) {
        data[i] = 10;
    }
}
template<typename T> void matrixRangeBasedForLoop(FullMatrixView<T> matrix) {
    for (auto& row : matrix) {
        for (auto& value : row) {
            value = 10;
        }
    }
}
template<typename T> void matrixIteratorLoop(T* data, int n, int m) {
    for (auto iter = data, end = data + n * m; iter != end; iter += m) {
        for (auto vIter = iter, vEnd = iter + m; vIter != vEnd; ++vIter) {
            *vIter = 10;
        }
    }
}
template<typename T> void matrixIndexLoop(T* data, int n, int m) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            data[j + i * n] = 10;
        }
    }
}
template<typename T> void matrixIndexIteratorHybridLoop(T* data, int n, int m) {
    T* iter = data;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            iter[j] = 10;
        }
        iter += n;
    }
}

template<typename T> void matrixIterationMethodsTest() {
    int n = 2000;
    int m = 2000;
    int times = 300;
    FullMatrix<T> v(n, m);

    auto matrixDirectIndexingTime = benchmark(times, matrixDirectIndexing<T>, v.data(), v.rowCount(), v.columnCount());
    auto matrixRangeBasedForLoopTime = benchmark(times, matrixRangeBasedForLoop<T>, v);
    auto matrixIteratorLoopTime = benchmark(times, matrixIteratorLoop<T>, v.data(), v.rowCount(), v.columnCount());
    auto matrixIndexLoopTime = benchmark(times, matrixIndexLoop<T>, v.data(), v.rowCount(), v.columnCount());
    auto matrixIndexIteratorHybridLoopTime = benchmark(times, matrixIndexIteratorHybridLoop<T>, v.data(), v.rowCount(), v.columnCount());

    matrixDirectIndexingTime = benchmark(times, matrixDirectIndexing<T>, v.data(), v.rowCount(), v.columnCount());
    matrixRangeBasedForLoopTime = benchmark(times, matrixRangeBasedForLoop<T>, v);
    matrixIteratorLoopTime = benchmark(times, matrixIteratorLoop<T>, v.data(), v.rowCount(), v.columnCount());
    matrixIndexLoopTime = benchmark(times, matrixIndexLoop<T>, v.data(), v.rowCount(), v.columnCount());
    matrixIndexIteratorHybridLoopTime = benchmark(times, matrixIndexIteratorHybridLoop<T>, v.data(), v.rowCount(), v.columnCount());

    std::cout << "matrixDirectIndexingTime          : " << matrixDirectIndexingTime << '\n';
    std::cout << "matrixRangeBasedForLoopTime       : " << matrixRangeBasedForLoopTime << '\n';
    std::cout << "matrixIteratorLoopTime            : " << matrixIteratorLoopTime << '\n';
    std::cout << "matrixIndexLoopTime               : " << matrixIndexLoopTime << '\n';
    std::cout << "matrixIndexIteratorHybridLoopTime : " << matrixIndexIteratorHybridLoopTime << '\n';
}

int main() {
    strassenVsNaiveMulSpeedTest();
    fast3x3VsNaiveMulSpeedTest();
    matrixIterationMethodsTest<int>();

    std::cin.get();
    return 0;
}
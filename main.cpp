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

template<int N, int Steps, bool OnlyAvx=false> void strassenTest(int n, int steps) {
    if (n != N || steps != Steps) {
        std::cout << "incorrect stdin\n";
        return;
    }
    std::cout << std::left;
    std::cout << "=================================================================\n";
    std::cout << "N               : " << N << '\n';
    std::cout << "Recursive Steps : " << Steps << '\n';
    std::cout << "-------------------------------\n";
    if (!OnlyAvx) {
        {
            HeapMatrix<int> a(n, n);
            HeapMatrix<int> b(n, n);
            for (int i = 0; i < a.size(); ++i) {
                a.data()[i] = i;
                b.data()[i] = i * 2;
            }
            auto start = MyTime();
            auto c = strassen(a, b, steps);
            auto end = MyTime();
            std::cout << "Dynamic High-level noAVX : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
        }
        {
            HeapMatrix<int> a(n, n);
            HeapMatrix<int> b(n, n);
            for (int i = 0; i < a.size(); ++i) {
                a.data()[i] = i;
                b.data()[i] = i * 2;
            }
            auto start = MyTime();
            auto c = lowLevelStrassen(a, b, steps);
            auto end = MyTime();
            std::cout << "Static  High-level noAVX : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
        }
        {
            StaticHeapMatrix<int, N, N> a;
            StaticHeapMatrix<int, N, N> b;
            for (int i = 0; i < a.size(); ++i) {
                a.data()[i] = i;
                b.data()[i] = i * 2;
            }
            auto start = MyTime();
            auto c = strassen<Steps>(a, b);
            auto end = MyTime();
            std::cout << "Dynamic Low-level  noAVX : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
        }
        {
            StaticHeapMatrix<int, N, N> a;
            StaticHeapMatrix<int, N, N> b;
            for (int i = 0; i < a.size(); ++i) {
                a.data()[i] = i;
                b.data()[i] = i * 2;
            }
            auto start = MyTime();
            auto c = lowLevelStrassen(a, b, Steps);
            auto end = MyTime();
            std::cout << "Static  Low-level  noAVX : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
        }
    }
    {
        HeapMatrix<int> a(n, n);
        HeapMatrix<int> b(n, n);
        for (int i = 0; i < a.size(); ++i) {
            a.data()[i] = i;
            b.data()[i] = i*2;
        }
        auto start = MyTime();
        auto c = strassenAvx(a, b, steps);
        auto end = MyTime();
        std::cout << "Dynamic High-level AVX   : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
    }
    {
        HeapMatrix<int> a(n, n);
        HeapMatrix<int> b(n, n);
        for (int i = 0; i < a.size(); ++i) {
            a.data()[i] = i;
            b.data()[i] = i * 2;
        }
        auto start = MyTime();
        auto c = lowLevelAvxStrassen(a, b, steps);
        auto end = MyTime();
        std::cout << "Static  High-level AVX   : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
    }
    {
        StaticHeapMatrix<int, N, N> a;
        StaticHeapMatrix<int, N, N> b;
        for (int i = 0; i < a.size(); ++i) {
            a.data()[i] = i;
            b.data()[i] = i * 2;
        }
        auto start = MyTime();
        auto c = strassenAvx<Steps>(a, b);
        auto end = MyTime();
        std::cout << "Dynamic Low-level  AVX   : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
    }
    {
        StaticHeapMatrix<int, N, N> a;
        StaticHeapMatrix<int, N, N> b;
        for (int i = 0; i < a.size(); ++i) {
            a.data()[i] = i;
            b.data()[i] = i * 2;
        }
        auto start = MyTime();
        auto c = lowLevelAvxStrassen(a, b, Steps);
        auto end = MyTime();
        std::cout << "Static  Low-level  AVX   : " << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
    }
    std::cout << "=================================================================\n\n";
}

int main() {
    std::array<int, 4> vec;
    std::cout << "copy line below and click enter: \n";
    std::cout << "480 992 1024 1920 2048 4096\n";
    std::cin >> vec[0] >> vec[1] >> vec[2] >> vec[3] >> vec[4] >> vec[5];
    strassenTest<480, 0>(vec[0], 0);
    strassenTest<480, 1>(vec[0], 1);
    strassenTest<480, 2>(vec[0], 2);
    strassenTest<992, 0>(vec[1], 0);
    strassenTest<992, 1>(vec[1], 1);
    strassenTest<992, 2>(vec[1], 2);
    strassenTest<992, 3>(vec[1], 3);
    strassenTest<992, 4>(vec[1], 4);
    strassenTest<992, 5>(vec[1], 5);
    strassenTest<1024, 0>(vec[2], 0);
    strassenTest<1024, 1>(vec[2], 1);
    strassenTest<1024, 2>(vec[2], 2);
    strassenTest<1920, 0>(vec[3], 0);
    strassenTest<1920, 1>(vec[3], 1);
    strassenTest<1920, 2>(vec[3], 2);
    strassenTest<1920, 3>(vec[3], 3);
    strassenTest<1920, 4>(vec[3], 4);
    strassenTest<1920, 5>(vec[3], 5);
    strassenTest<2048, 2>(vec[4], 2);
    strassenTest<2048, 3>(vec[4], 3);
    strassenTest<4096, 0, true>(vec[5], 0);
    strassenTest<4096, 1, true>(vec[5], 1);
    strassenTest<4096, 2, true>(vec[5], 2);

    std::cin.get();
    return 0;
}

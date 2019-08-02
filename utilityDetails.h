/*#pragma once
#include <cstdint>
#include <cmath>
#include <array>
#include "debugAssert.h"

namespace details {
    double log(double base, double value) {
        return ::log(value) / ::log(base);
    }
    int numberOfSteps(int size, int base, int cutOffPoint) {
        if (size <= cutOffPoint) return 0;
        return static_cast<int>(ceil(log(base, static_cast<double>(size) / cutOffPoint)));
    }
    int staticPaddingNewSize(int size, int base, int steps) {
        if (steps <= 0) return size;
        int divisor = pow(base, steps);
        auto remainder = size % divisor;
        if (remainder == 0) return size;
        else                return size + divisor - remainder;
    }

    template<
        int BaseSize1, 
        int BaseSize2,
        int BaseSize3,
        typename T, 
        template<typename> typename MatrixA, 
        template<typename> typename MatrixB, 
        typename Function
    >
    FullMatrix<T> fastMul(const MatrixA<T>& a, const MatrixB<T>& b, int cutOffPoint, Function recursiveFunction) {
        debugAssertOp(a.columnCount(), ==, b.rowCount());

        int steps = details::numberOfSteps(a.columnCount(), BaseSize2, cutOffPoint);

        std::array<int, 3> sizes = { a.rowCount(), a.columnCount(), b.columnCount() };
        std::array<int, 3> paddedSizes = {
            details::staticPaddingNewSize(sizes[0], BaseSize1, steps),
            details::staticPaddingNewSize(sizes[1], BaseSize2, steps),
            details::staticPaddingNewSize(sizes[2], BaseSize3, steps)
        };

        if (paddedSizes[0] > sizes[0] || paddedSizes[1] > sizes[1] || paddedSizes[2] > sizes[2]) {
            FullMatrix<T> newA(paddedSizes[0], paddedSizes[1]);
            FullMatrix<T> newB(paddedSizes[1], paddedSizes[2]);
            newA.copy(a);
            newB.copy(b);
            auto result = recursiveFunction(newA, newB, steps);
            result.shrink(sizes[0], sizes[2]);
            return result;
        } else {
            return recursiveFunction(a, b, steps);
        }
    }

    int nextPowerOf2(uint32_t n) {
        n -= 1;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        return n + 1;
    }

    bool isPowerOf2(int n) {
        return n && !(n & (n - 1));
    }

    int nextPowerOf(int number, int power) {
        int nextPower = 1;
        while (nextPower < number) {
            nextPower *= power;
        }
        return nextPower;
    }
    bool isPowerOf(int number, int power) {
        double value = std::log(power) / std::log(number);
        return trunc(value) == value;
    }
}*/
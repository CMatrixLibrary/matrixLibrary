#pragma once
#include <tuple>
#include "FullMatrix.h"
#include "FullSubMatrix.h"
#include "FullMatrixView.h"
#include "naiveOperations.h"
#include "operators.h"
#include "matrixDivision.h"

namespace details {
    template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
    FullMatrix<T> strassenMulRecursive(const MatrixA<T>& a, const MatrixB<T>& b) {
        if (a.rowsCount() <= 32) {
            return naiveMul(a, b);
        }

        auto dA = matrixDivide<2, 2>(a);
        auto dB = matrixDivide<2, 2>(b);

        auto m1 = strassenMulRecursive(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1]);
        auto m2 = strassenMulRecursive(dA[1][0] + dA[1][1], dB[0][0]);
        auto m3 = strassenMulRecursive(dA[0][0], dB[0][1] - dB[1][1]);
        auto m4 = strassenMulRecursive(dA[1][1], dB[1][0] - dB[0][0]);
        auto m5 = strassenMulRecursive(dA[0][0] + dA[0][1], dB[1][1]);
        auto m6 = strassenMulRecursive(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1]);
        auto m7 = strassenMulRecursive(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1]);

        auto c11 = m1 + m4 - m5 + m7;
        auto c12 = m3 + m5;
        auto c21 = m2 + m4;
        auto c22 = m1 + m3 - m2 + m6;

        auto n = c11.rowsCount();
        FullMatrix<T> c(n * 2, n * 2);
        for (int y = 0; y < n; y++) {
            for (int x = 0; x < n; x++) {
                c.at(x, y) = c11.at(x, y);
                c.at(x + n, y) = c12.at(x, y);
                c.at(x, y + n) = c21.at(x, y);
                c.at(x + n, y + n) = c22.at(x, y);
            }
        }
        return c;
    }
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB> 
FullMatrix<T> strassenMul(MatrixA<T> a, MatrixB<T> b) {
    if (!details::isPowerOf2(a.rowsCount())) {
        auto newSize = details::nextPowerOf2(a.rowsCount());
        FullMatrix<T> newA(newSize, newSize);
        FullMatrix<T> newB(newSize, newSize);
        newA.insert(a);
        newB.insert(b);
        auto result = details::strassenMulRecursive(newA, newB);
        result.shrink(a.rowsCount(), a.rowsCount());
        return result;
    } else {
        return details::strassenMulRecursive(a, b);
    }
}
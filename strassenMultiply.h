#pragma once
#include <tuple>
#include "FullMatrix.h"
#include "FullMatrixView.h"
#include "FullMatrixConstView.h"
#include "utilityDetails.h"
#include "naiveOperations.h"
#include "operators.h"
#include "matrixDivision.h"

/*
    "Gaussian Elimination is not Optimal"
    Volker Strassen
    https://link.springer.com/article/10.1007%2FBF02165411
*/


namespace details {
    template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
    FullMatrix<T> strassenMulRecursive(const MatrixA<T>& a, const MatrixB<T>& b) {
        if (a.rowCount() <= 32) {
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

        FullMatrix<T> c(a.rowCount(), b.columnCount());
        auto dC = matrixDivide<2, 2>(c);
        dC[0][0].copy(m1 + m4 - m5 + m7);
        dC[0][1].copy(m3 + m5);
        dC[1][0].copy(m2 + m4);
        dC[1][1].copy(m1 + m3 - m2 + m6);

        return c;
    }
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB> 
FullMatrix<T> strassenMul(MatrixA<T> a, MatrixB<T> b) {
    if (!details::isPowerOf2(a.rowCount())) {
        auto newSize = details::nextPowerOf2(a.rowCount());
        FullMatrix<T> newA(newSize, newSize);
        FullMatrix<T> newB(newSize, newSize);
        newA.copy(a);
        newB.copy(b);
        auto result = details::strassenMulRecursive(newA, newB);
        result.shrink(a.rowCount(), a.rowCount());
        return result;
    } else {
        return details::strassenMulRecursive(a, b);
    }
}
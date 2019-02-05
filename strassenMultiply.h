#pragma once
#include <tuple>
#include "FullMatrix.h"
#include "FullSubMatrix.h"
#include "FullMatrixView.h"
#include "naiveOperations.h"
#include "operators.h"

namespace details {

    template<typename T, template<typename> typename Matrix>
    std::tuple<FullMatrixView<T>, FullMatrixView<T>, FullMatrixView<T>, FullMatrixView<T>> 
    strassenDivide(const Matrix<T>& m) {
        auto halfRows = m.rowsCount() / 2;
        auto halfCols = m.columnsCount() / 2;
        return {
            FullMatrixView(m, 0,        0,        halfRows, halfCols),
            FullMatrixView(m, 0,        halfCols, halfRows, halfCols),
            FullMatrixView(m, halfRows, 0,        halfRows, halfCols),
            FullMatrixView(m, halfRows, halfCols, halfRows, halfCols)
        };
    }

    template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
    FullMatrix<T> strassenMulRecursive(const MatrixA<T>& a, const MatrixB<T>& b) {
        if (a.rowsCount() <= 32) {
            return naiveMul(a, b);
        }

        auto[a11, a12, a21, a22] = strassenDivide(a);
        auto[b11, b12, b21, b22] = strassenDivide(b);

        auto m1 = strassenMulRecursive(a11 + a22, b11 + b22);
        auto m2 = strassenMulRecursive(a21 + a22, b11);
        auto m3 = strassenMulRecursive(a11, b12 - b22);
        auto m4 = strassenMulRecursive(a22, b21 - b11);
        auto m5 = strassenMulRecursive(a11 + a12, b22);
        auto m6 = strassenMulRecursive(a21 - a11, b11 + b12);
        auto m7 = strassenMulRecursive(a12 - a22, b21 + b22);

        auto c11 = m1 + m4 - m5 + m7;
        auto c12 = m3 + m5;
        auto c21 = m2 + m4;
        auto c22 = m1 - m2 + m3 + m6;

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
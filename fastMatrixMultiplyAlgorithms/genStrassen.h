// this file was generated using fastMatrixMultiplyAlgorithms/generator
/*#pragma once
#include <array>
#include "../FullMatrix.h"
#include "../FullMatrixView.h"
#include "../FullMatrixConstView.h"
#include "../naiveOperations.h"
#include "../operators.h"
#include "../utilityDetails.h"

namespace details {
    template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
    FullMatrix<T> genStrassen(const MatrixA<T>& A, const MatrixB<T>& B, int steps) {
        if (steps <= 0) {
            return naiveMul(A, B);
        }

        auto a = matrixDivide<2, 2>(A);
        auto b = matrixDivide<2, 2>(B);

        std::array<FullMatrix<T>, 7> m;
        m[0] = genStrassen(a[0][0] + a[1][1], b[0][0] + b[1][1], steps - 1);
        m[1] = genStrassen(a[1][0] + a[1][1], b[0][0], steps - 1);
        m[2] = genStrassen(a[0][0], b[0][1] - b[1][1], steps - 1);
        m[3] = genStrassen(a[1][1], -b[0][0] + b[1][0], steps - 1);
        m[4] = genStrassen(a[0][0] + a[0][1], b[1][1], steps - 1);
        m[5] = genStrassen(-a[0][0] + a[1][0], b[0][0] + b[0][1], steps - 1);
        m[6] = genStrassen(a[0][1] - a[1][1], b[1][0] + b[1][1], steps - 1);

        FullMatrix<T> C(A.rowCount(), B.columnCount());
        auto c = matrixDivide<2, 2>(C);

        c[0][0].copy(m[0] + m[3] - m[4] + m[6]);
        c[0][1].copy(m[2] + m[4]);
        c[1][0].copy(m[1] + m[3]);
        c[1][1].copy(m[0] - m[1] + m[2] + m[5]);

        return C;
    }
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> genStrassen(const MatrixA<T>& A, const MatrixB<T>& B) {
    return details::fastMul<2, 2, 2>(A, B, 50, details::genStrassen<T, MatrixA, MatrixB>);
}
*/
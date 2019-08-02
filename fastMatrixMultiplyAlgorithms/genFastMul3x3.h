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
    FullMatrix<T> genFastMul3x3(const MatrixA<T>& A, const MatrixB<T>& B, int steps) {
        if (steps <= 0) {
            return naiveMul(A, B);
        }

        auto a = matrixDivide<3, 3>(A);
        auto b = matrixDivide<3, 3>(B);

        std::array<FullMatrix<T>, 23> m;
        auto a_11_15 = a[2][1] + a[2][2];
        auto a_8_10 = a[2][0] - a[0][0];
        auto a_0_1 = a[0][0] + a[0][1];
        auto a_2_19 = a[0][2] - a_0_1;
        auto a_5_20 = -a[2][1] - a_2_19;
        auto a_7_16 = a[1][1] + a[1][2];
        auto a_4_21 = -a[1][1] - a_5_20;
        auto a_9_7 = a[1][0] + a[1][1];
        auto b_12_17 = b[1][2] - b[2][2];
        auto b_2_7 = -b[0][0] - b[2][0];
        auto b_4_19 = b[1][0] - b_2_7;
        auto b_0_16 = b[1][1] - b[2][1];
        auto b_6_20 = -b[1][2] - b_4_19;
        auto b_5_22 = -b[1][1] - b_6_20;
        m[0] = genFastMul3x3(-a[1][0] - a[2][2] - a_4_21, b[1][1], steps - 1);
        m[1] = genFastMul3x3(a[0][0] - a[1][0], -b[1][0] + b[1][1], steps - 1);
        m[2] = genFastMul3x3(a[1][1], b[0][1] + b[2][2] - b_5_22, steps - 1);
        m[3] = genFastMul3x3(-a[0][0] - a_9_7, b[0][0] - b[0][1] + b[1][1], steps - 1);
        m[4] = genFastMul3x3(a_9_7, -b[0][0] + b[0][1], steps - 1);
        m[5] = genFastMul3x3(a[0][0], b[0][0], steps - 1);
        m[6] = genFastMul3x3(a[2][1] - a_8_10, b[0][0] - b[0][2] + b[1][2], steps - 1);
        m[7] = genFastMul3x3(a_8_10, b[0][2] - b[1][2], steps - 1);
        m[8] = genFastMul3x3(a[2][0] + a[2][1], -b[0][0] + b[0][2], steps - 1);
        m[9] = genFastMul3x3(-a[1][2] - a[2][0] - a_4_21, b[1][2], steps - 1);
        m[10] = genFastMul3x3(a[2][1], b[0][2] + b[2][1] - b_5_22, steps - 1);
        m[11] = genFastMul3x3(-a[0][2] - a_11_15, b[2][0] - b_0_16, steps - 1);
        m[12] = genFastMul3x3(a[0][2] - a[2][2], b_0_16, steps - 1);
        m[13] = genFastMul3x3(a[0][2], b[2][0], steps - 1);
        m[14] = genFastMul3x3(a_11_15, -b[2][0] + b[2][1], steps - 1);
        m[15] = genFastMul3x3(-a[0][2] - a_7_16, b[2][0] - b_12_17, steps - 1);
        m[16] = genFastMul3x3(a[0][2] - a[1][2], b_12_17, steps - 1);
        m[17] = genFastMul3x3(a_7_16, -b[2][0] + b[2][2], steps - 1);
        m[18] = genFastMul3x3(a[0][1], b[1][0], steps - 1);
        m[19] = genFastMul3x3(a[1][2], b[2][1], steps - 1);
        m[20] = genFastMul3x3(a[1][0], b[0][2], steps - 1);
        m[21] = genFastMul3x3(a[2][0], b[0][1], steps - 1);
        m[22] = genFastMul3x3(a[2][2], b[2][2], steps - 1);

        FullMatrix<T> C(A.rowCount(), B.columnCount());
        auto c = matrixDivide<3, 3>(C);

        c[0][0].copy(m[5] + m[13] + m[18]);
        c[0][1].copy(m[0] + m[3] + m[4] + m[5] + m[11] + m[13] + m[14]);
        c[0][2].copy(m[5] + m[6] + m[8] + m[9] + m[13] + m[15] + m[17]);
        c[1][0].copy(m[1] + m[2] + m[3] + m[5] + m[13] + m[15] + m[16]);
        c[1][1].copy(m[1] + m[2] + m[3] + m[5] + m[13] + m[15] + m[16]);
        c[1][2].copy(m[13] + m[15] + m[16] + m[17] + m[20]);
        c[2][0].copy(m[5] + m[6] + m[7] + m[10] + m[11] + m[12] + m[13]);
        c[2][1].copy(m[11] + m[12] + m[13] + m[14] + m[21]);
        c[2][2].copy(m[5] + m[6] + m[7] + m[8] + m[22]);

        return C;
    }
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> genFastMul3x3(const MatrixA<T>& A, const MatrixB<T>& B) {
    return details::fastMul<3, 3, 3>(A, B, 50, details::genFastMul3x3<T, MatrixA, MatrixB>);
}
*/
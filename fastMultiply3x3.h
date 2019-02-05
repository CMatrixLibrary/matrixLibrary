#pragma once
#include <tuple>
#include "FullMatrix.h"
#include "FullMatrixView.h"
#include "naiveOperations.h"
#include "operators.h"
#include "matrixDivision.h"

/*
    A noncommutative algorithm for multiplying 3 x 3 matrices using 23 multiplications
    Julian D. Laderman 
    http://www.ams.org/journals/bull/1976-82-01/S0002-9904-1976-13988-2/S0002-9904-1976-13988-2.pdf
*/

namespace details {
    template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
    FullMatrix<T> fastMul3x3Recursive(const MatrixA<T>& a, const MatrixB<T>& b) {
        if (a.rowsCount() <= 27) {
            return naiveMul(a, b);
        }

        auto dA = matrixDivide<3, 3>(a);
        auto dB = matrixDivide<3, 3>(b);
        auto[a11, a12, a13, a21, a22, a23, a31, a32, a33] = std::make_tuple(
            dA[0][0], dA[0][1], dA[0][2],
            dA[1][0], dA[1][1], dA[1][2],
            dA[2][0], dA[2][1], dA[2][2]
        );
        auto[b11, b12, b13, b21, b22, b23, b31, b32, b33] = std::make_tuple(
            dB[0][0], dB[0][1], dB[0][2],
            dB[1][0], dB[1][1], dB[1][2],
            dB[2][0], dB[2][1], dB[2][2]
        );

        auto m1  = fastMul3x3Recursive(a11 + a12 + a13 - a21 - a22 - a32 - a33, b22);
        auto m2  = fastMul3x3Recursive(a11 - a21, b22 - b12);
        auto m3  = fastMul3x3Recursive(a22, b12 - b11 + b21 - b22 - b23 - b31 + b33);
        auto m4  = fastMul3x3Recursive(a21 - a11 + a22, b11 - b12 + b22);
        auto m5  = fastMul3x3Recursive(a21 + a22, b12 - b11);
        auto m6  = fastMul3x3Recursive(a11, b11);
        auto m7  = fastMul3x3Recursive(a31 - a11 + a32, b11 - b13 + b23);
        auto m8  = fastMul3x3Recursive(a31 - a11, b13 - b23);
        auto m9  = fastMul3x3Recursive(a31 + a32, b13 - b11);
        auto m10 = fastMul3x3Recursive(a11 + a12 + a13 - a22 - a23 - a31 - a32, b23);
        auto m11 = fastMul3x3Recursive(a32, b13 - b11 + b21 - b22 - b23 - b31 + b32);
        auto m12 = fastMul3x3Recursive(a32 - a13 + a33, b22 + b31 - b32);
        auto m13 = fastMul3x3Recursive(a13 - a33, b22 - b32);
        auto m14 = fastMul3x3Recursive(a13, b31);
        auto m15 = fastMul3x3Recursive(a32 + a33, b32 - b31);
        auto m16 = fastMul3x3Recursive(a22 - a13 + a23, b23 + b31 - b33);
        auto m17 = fastMul3x3Recursive(a13 - a23, b23 - b33);
        auto m18 = fastMul3x3Recursive(a22 + a23, b33 - b31);
        auto m19 = fastMul3x3Recursive(a12, b21);
        auto m20 = fastMul3x3Recursive(a23, b32);
        auto m21 = fastMul3x3Recursive(a21, b13);
        auto m22 = fastMul3x3Recursive(a31, b12);
        auto m23 = fastMul3x3Recursive(a33, b33);

        auto c11 = m6 + m14 + m19;
        auto c12 = m1 + m4 + m5 + m6 + m12 + m14 + m15;
        auto c13 = m6 + m7 + m9 + m10 + m14 + m16 + m18;
        auto c21 = m2 + m3 + m4 + m6 + m14 + m16 + m17;
        auto c22 = m2 + m4 + m5 + m6 + m20;
        auto c23 = m14 + m16 + m17 + m18 + m21;
        auto c31 = m6 + m7 + m8 + m11 + m12 + m13 + m14;
        auto c32 = m12 + m13 + m14 + m15 + m22;
        auto c33 = m6 + m7 + m8 + m9 + m23;

        auto n = c11.rowsCount();
        FullMatrix<T> c(n * 3, n * 3);
        for (int y = 0; y < n; y++) {
            for (int x = 0; x < n; x++) {
                c.at(x,       y)       = c11.at(x, y);
                c.at(x + n,   y)       = c12.at(x, y);
                c.at(x + 2*n, y)       = c13.at(x, y);
                c.at(x,       y + n)   = c21.at(x, y);
                c.at(x + n,   y + n)   = c22.at(x, y);
                c.at(x + 2*n, y + n)   = c23.at(x, y);
                c.at(x + 2*n, y + 2*n) = c31.at(x, y);
                c.at(x + 2*n, y + 2*n) = c32.at(x, y);
                c.at(x + 2*n, y + 2*n) = c33.at(x, y);
            }
        }

        return c;
    }
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB> 
FullMatrix<T> fastMul3x3(MatrixA<T> a, MatrixB<T> b) {
    if (!details::isPowerOf(a.rowsCount(), 3)) {
        auto newSize = details::nextPowerOf(a.rowsCount(), 3);
        FullMatrix<T> newA(newSize, newSize);
        FullMatrix<T> newB(newSize, newSize);
        newA.insert(a);
        newB.insert(b);
        auto result = details::fastMul3x3Recursive(newA, newB);
        result.shrink(a.rowsCount(), a.rowsCount());
        return result;
    } else {
        return details::fastMul3x3Recursive(a, b);
    }
}
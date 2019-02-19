#pragma once
#include "FullMatrix.h"
#include "debugAssert.h"

template<typename T> class FullMatrix;

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> naiveMul(const MatrixA<T>& a, const MatrixB<T>& b) {
    debugAssertOp(a.columnCount(), == , b.rowCount());

    FullMatrix<T> result(a.rowCount(), b.columnCount());
    for (int i = 0; i < a.rowCount(); ++i) {
        for (int j = 0; j < b.columnCount(); ++j) {
            for (int k = 0; k < a.columnCount(); ++k) {
                result.at(i, j) += a.at(i, k) * b.at(k, j);
            }
        }
    }
    return result;
}

template<typename T, template<typename> typename Matrix>
FullMatrix<T> naiveMul(const Matrix<T>& m, T scalar) {
    FullMatrix<T> result(m.rowCount(), m.columnCount());
    for (int row = 0; row < result.rowCount(); ++row) {
        for (int column = 0; column < result.columnCount(); ++column) {
            result.at(row, column) += scalar * m.at(row, column);
        }
    }
    return result;
}
template<typename T, template<typename> typename Matrix>
FullMatrix<T> naiveMul(T scalar, const Matrix<T>& m) {
    return naiveMul(m, scalar);
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> naiveAdd(const MatrixA<T>& a, const MatrixB<T>& b) {
    debugAssertOp(a.rowCount(), ==, b.rowCount());
    debugAssertOp(a.columnCount(), ==, b.columnCount());

    FullMatrix<T> result(a.rowCount(), a.columnCount());
    for (int row = 0; row < result.rowCount(); ++row) {
        for (int column = 0; column < result.columnCount(); ++column) {
            result.at(row, column) = a.at(row, column) + b.at(row, column);
        }
    }
    return result;
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> naiveSub(const MatrixA<T>& a, const MatrixB<T>& b) {
    debugAssertOp(a.rowCount(), ==, b.rowCount());
    debugAssertOp(a.columnCount(), ==, b.columnCount());

    FullMatrix<T> result(a.rowCount(), a.columnCount());
    for (int row = 0; row < result.rowCount(); ++row) {
        for (int column = 0; column < result.columnCount(); ++column) {
            result.at(row, column) = a.at(row, column) - b.at(row, column);
        }
    }
    return result;
}

template<typename T, template<typename> typename Matrix>
FullMatrix<T> naiveNeg(const Matrix<T>& m) {
    FullMatrix<T> result(m.rowCount(), m.columnCount());
    for (int row = 0; row < result.rowCount(); ++row) {
        for (int column = 0; column < result.columnCount(); ++column) {
            result.at(row, column) = -m.at(row, column);
        }
    }
    return result;
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
void naiveAddAssign(MatrixA<T>& a, const MatrixB<T>& b) {
    debugAssertOp(a.rowCount(), ==, b.rowCount());
    debugAssertOp(a.columnCount(), ==, b.columnCount());

    for (int row = 0; row < a.rowCount(); ++row) {
        for (int column = 0; column < a.columnCount(); ++column) {
            a.at(row, column) += b.at(row, column);
        }
    }
}
template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
void naiveSubAssign(MatrixA<T>& a, const MatrixB<T>& b) {
    debugAssertOp(a.rowCount(), ==, b.rowCount());
    debugAssertOp(a.columnCount(), ==, b.columnCount());

    for (int row = 0; row < a.rowCount(); ++row) {
        for (int column = 0; column < a.columnCount(); ++column) {
            a.at(row, column) -= b.at(row, column);
        }
    }
}
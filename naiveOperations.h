#pragma once
#include "FullMatrix.h"

template<typename T> class FullMatrix;

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> naiveMul(const MatrixA<T>& a, const MatrixB<T>& b) {
    FullMatrix<T> result(a.rowsCount(), b.columnsCount());
    for (int i = 0; i < a.rowsCount(); ++i) {
        for (int j = 0; j < b.columnsCount(); ++j) {
            for (int k = 0; k < a.columnsCount(); ++k) {
                result.at(j, i) += a.at(k, i) * b.at(j, k);
            }
        }
    }
    return result;
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> naiveAdd(const MatrixA<T>& a, const MatrixB<T>& b) {
    FullMatrix<T> result(a.rowsCount(), a.columnsCount());
    for (int i = 0; i < result.data().size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> naiveSub(const MatrixA<T>& a, const MatrixB<T>& b) {
    FullMatrix<T> result(a.rowsCount(), a.columnsCount());
    for (int i = 0; i < result.data().size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}
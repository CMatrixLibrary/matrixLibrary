#ifndef PARALLEL_MUL_H
#define PARALLEL_MUL_H

#include "MatrixInterface.h"
#include "Matrix.h"

template<typename MR, typename M1, typename M2>
void parallelMul(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
#pragma omp parallel for
    for (mtl::size_t i = 0; i < a.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < b.columnCount(); ++j) {
            typename MR::ValueType sum{};
            for (mtl::size_t k = 0; k < a.columnCount(); ++k) {
                sum += a.at(i, k) * b.at(k, j);
            }
            r.at(i, j) = sum;
        }
    }
}

template<typename M1, typename M2>
auto parallelMul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        parallelMul(result, a, b);
        return result;
    } else {
        Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
        parallelMul(result, a, b);
        return result;
    }
}

#endif

#pragma once
#include "Matrix.h"
#include "blasUtility.h"

template<typename M1, typename M2>
auto blasMul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (blas::IsAvailable) {
        if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
            Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
            for (int i = 0; i < a.size(); ++i) result.data()[i] = typename M1::ValueType{};
            blas::mul(result.data(), a.data(), b.data(), M1::CRow(), M1::CCol(), M2::CCol(), result.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
            return result;
        } else {
            Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
            for (int i = 0; i < a.size(); ++i) result.data()[i] = typename M1::ValueType{};
            blas::mul(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), result.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
            return result;
        }
    } else {
        static_assert(blas::IsAvailable && always_false_v<M1>, "blas is not available");
    }
}
#ifndef MATRIX_OPERATORS_H
#define MATRIX_OPERATORS_H
#include "naiveBasicOperations.h"
#include "avxMul.h"
#include "blasMul.h"
#include "parallelBlockMul.h"
#include <iostream>
#include <iomanip>

template<typename M1, typename M2> auto operator+(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return naiveAdd(a, b);
}
template<typename M1, typename M2> auto& operator+=(MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return naiveAddAssign(a, b);
}
template<typename M1, typename M2> auto operator-(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return naiveSub(a, b);
}
template<typename M1, typename M2> auto& operator-=(MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return naiveSubAssign(a, b);
}
template<typename M> auto operator-(const MatrixInterface<M>& m) {
    return naiveNeg(m);
}

template<typename M1, typename M2> auto operator*(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (blas::IsAvailable && blas::IsCompatible<typename M1::ValueType> && blas::IsCompatible<typename M2::ValueType>)
        return blas::mul(a, b);
    if constexpr (avx::IsAvailable && avx::IsCompatible<typename M1::ValueType> && avx::IsCompatible<typename M2::ValueType>)
        return avx::parallelMul(a, b);
    return parallelBlockMul(a, b);
}
template<typename M> auto operator*(const MatrixInterface<M>& m, const typename M::ValueType& scalar) {
    return naiveScalarMul(m, scalar);
}
template<typename M> auto& operator*=(MatrixInterface<M>& m, const typename M::ValueType& scalar) {
    return naiveScalarMulAssign(m, scalar);
}

template<typename MT> std::ostream& operator<<(std::ostream& out, const MatrixInterface<MT>& m) {
    for (auto row : m) {
        for (auto value : row) {
            out << std::setw(6) << value << ' ';
        }
        out << '\n';
    }
    return out;
}
#endif

#ifndef NAIVE_MUL_H
#define NAIVE_MUL_H
#include "MatrixInterface.h"
#include "Matrix.h"

template<typename T> void naiveMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            T sum{};
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * effA] * b[j + k * effB];
            }
            c[j + i * effC] = sum;
        }
    }
}

template<typename T> void naiveMul(T* c, const T* a, const T* b, int n, int m, int p) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            T sum{};
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * m] * b[j + k * p];
            }
            c[j + i * p] = sum;
        }
    }
}
template<int n, int m, int p, int effC, int effA, int effB, typename T> void naiveMul(T* c, const T* a, const T* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            T sum{};
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * effA] * b[j + k * effB];
            }
            c[j + i * effC] = sum;
        }
    }
}
template<int n, int m, int p, typename T> void naiveMul(T* c, const T* a, const T* b) {
    return naiveMul<n, m, p, p, m, p>(c, a, b);
}

template<typename MR, typename M1, typename M2>
void naiveMul(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
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
auto naiveMul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        naiveMul(result, a, b);
        return result;
    } else {
        Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
        naiveMul(result, a, b);
        return result;
    }
}

#endif

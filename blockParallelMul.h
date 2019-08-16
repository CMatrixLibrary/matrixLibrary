#ifndef BLOCK_PARALLEL_MUL_H
#define BLOCK_PARALLEL_MUL_H

#include "MatrixInterface.h"
#include "Matrix.h"

template<typename MR, typename M1, typename M2>
void blockParallelMul(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    constexpr mtl::size_t ib = 16;
    constexpr mtl::size_t jb = 16;
    constexpr mtl::size_t kb = 16;
    auto n = a.rowCount();
    auto m = a.columnCount();
    auto p = b.columnCount();

    for (mtl::size_t i = 0; i < n; ++i) {
        for (mtl::size_t j = 0; j < p; ++j) {
            r.at(i, j) = typename MR::ValueType{};
        }
    }

#pragma omp parallel for
    for (mtl::size_t ii = 0; ii < n; ii += ib) {
        mtl::size_t iEnd = std::min(ii + ib, n);
        for (mtl::size_t jj = 0; jj < m; jj += jb) {
            mtl::size_t jEnd = std::min(jj + jb, m);
            for (mtl::size_t kk = 0; kk < p; kk += kb) {
                mtl::size_t kEnd = std::min(kk + kb, p);
                for (auto i = ii; i < iEnd; ++i) {
                    for (auto j = jj; j < jEnd; ++j) {
                        for (auto k = kk; k < kEnd; ++k) {
                            r.at(i, j) += a.at(i, k) * b.at(k, j);
                        }
                    }
                }
            }
        }
    }
}

template<typename M1, typename M2>
auto blockParallelMul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        blockParallelMul(result, a, b);
        return result;
    } else {
        Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
        blockParallelMul(result, a, b);
        return result;
    }
}

#endif

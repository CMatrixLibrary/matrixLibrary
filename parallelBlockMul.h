#ifndef PARALLEL_BLOCK_MUL_H
#define PARALLEL_BLOCK_MUL_H

#include "MatrixInterface.h"
#include "Matrix.h"
#include "blockMul.h"
#include "ThreadPool.h"

template<typename T> void parallelBlockMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (int i = 0; i < n; i += chunkSize) {
        int nPart = std::min(chunkSize, n - i);
        pool.addTask([=]() { blockMul(c + i, a + i, b + i, nPart, m, p, effC, effA, effB); });
    }
}

template<typename T> void parallelBlockMul(T* c, const T* a, const T* b, int n, int m, int p) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (int i = 0; i < n; i += chunkSize) {
        int nPart = std::min(chunkSize, n - i);
        pool.addTask([=]() { blockMul(c + i, a + i, b + i, nPart, m, p); });
    }
}

template<int n, int m, int p, int effC, int effA, int effB, typename T> void parallelBlockMul(T* c, const T* a, const T* b) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (int i = 0; i < n; i += chunkSize) {
        constexpr int nPart = (chunkSize < n - i) ? chunkSize : n - i;
        pool.addTask([=]() { blockMul<nPart, m, p, effC, effA, effB>(c + i, a + i, b + i); });
    }
}
template<int n, int m, int p, typename T> void parallelBlockMul(T* c, const T* a, const T* b) {
    return parallelBlockMul<n, m, p, p, m, p>(c, a, b);
}


namespace detail {
    template<typename MR, typename M1, typename M2>
    void parallelBlockMulTask(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, mtl::size_t iStart, mtl::size_t iEnd) {
        constexpr mtl::size_t ib = 16;
        constexpr mtl::size_t jb = 16;
        constexpr mtl::size_t kb = 16;
        auto n = iEnd;
        auto m = a.columnCount();
        auto p = b.columnCount();

        for (mtl::size_t i = iStart; i < n; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                r.at(i, j) = typename MR::ValueType{};
            }
        }

        for (mtl::size_t ii = iStart; ii < n; ii += ib) {
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
}

template<typename MR, typename M1, typename M2>
void parallelBlockMul(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (mtl::size_t i = 0; i < a.rowCount(); i += chunkSize) {
        mtl::size_t partSize = std::min(chunkSize, a.rowCount() - i);
        pool.addTask([=]() { detail::parallelBlockMulTask(r, a, b, i, i + partSize); });
    }
}

template<typename M1, typename M2>
auto parallelBlockMul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        parallelBlockMul(result, a, b);
        return result;
    } else {
        Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
        parallelBlockMul(result, a, b);
        return result;
    }
}

#endif

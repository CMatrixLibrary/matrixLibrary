#ifndef PARALLEL_MUL_H
#define PARALLEL_MUL_H

#include "MatrixInterface.h"
#include "Matrix.h"
#include "ThreadPool.h"
#include "naiveMul.h"

template<typename T> void parallelMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (int i = 0; i < n; i += chunkSize) {
        int nPart = std::min(chunkSize, n - i);
        pool.addTask([=]() { naiveMul(c + i, a + i, b + i, nPart, m, p, effC, effA, effB); });
    }
}

template<typename T> void parallelMul(T* c, const T* a, const T* b, int n, int m, int p) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (int i = 0; i < n; i += chunkSize) {
        int nPart = std::min(chunkSize, n - i);
        pool.addTask([=]() { naiveMul(c + i, a + i, b + i, nPart, m, p); });
    }
}

template<int n, int m, int p, int effC, int effA, int effB, typename T> void parallelMul(T* c, const T* a, const T* b) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (int i = 0; i < n; i += chunkSize) {
        constexpr int nPart = (chunkSize < n - i) ? chunkSize : n - i;
        pool.addTask([=]() { naiveMul<nPart, m, p, effC, effA, effB>(c + i, a + i, b + i); });
    }
}
template<int n, int m, int p, typename T> void parallelMul(T* c, const T* a, const T* b) {
    return parallelMul<n, m, p, p, m, p>(c, a, b);
}


namespace detail {
    template<typename MR, typename M1, typename M2>
    void parallelMulTask(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, mtl::size_t iStart, mtl::size_t iEnd) {
        for (mtl::size_t i = iStart; i < iEnd; ++i) {
            for (mtl::size_t j = 0; j < b.columnCount(); ++j) {
                typename MR::ValueType sum{};
                for (mtl::size_t k = 0; k < a.columnCount(); ++k) {
                    sum += a.at(i, k) * b.at(k, j);
                }
                r.at(i, j) = sum;
            }
        }
    }
}

template<typename MR, typename M1, typename M2>
void parallelMul(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    constexpr int chunkSize = 128;

    ThreadPool pool;
    for (mtl::size_t i = 0; i < a.rowCount(); i += chunkSize) {
        mtl::size_t partSize = std::min(chunkSize, a.rowCount() - i);
        pool.addTask([=]() { detail::parallelMulTask(r, a, b, i, i + partSize); });
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

#ifndef PARALLEL_MUL_H
#define PARALLEL_MUL_H

#include "MatrixInterface.h"
#include "baseMulUtility.h"
#include "ThreadPool.h"

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
        pool.addTask([&r, &a, &b, i, partSize]() { detail::parallelMulTask(r, a, b, i, i + partSize); });
    }
}

CREATE_BASE_MUL_WRAPPERS(parallelMul)

#endif

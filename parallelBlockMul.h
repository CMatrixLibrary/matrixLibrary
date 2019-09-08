#ifndef PARALLEL_BLOCK_MUL_H
#define PARALLEL_BLOCK_MUL_H

#include "MatrixInterface.h"
#include "baseMulUtility.h"
#include "ThreadPool.h"

namespace detail {
    template<typename MR, typename M1, typename M2>
    void parallelBlockMulTask(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, mtl::size_t iStart, mtl::size_t iEnd) {
        constexpr mtl::size_t jb = 16;
        constexpr mtl::size_t kb = 256;

        auto n = a.rowCount();
        auto m = a.columnCount();
        auto p = b.columnCount();

        for (mtl::size_t i = iStart; i < iEnd; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                r.at(i, j) = typename MR::ValueType{};
            }
        }

        for (mtl::size_t kk = 0; kk < p; kk += kb) {
            mtl::size_t kEnd = std::min(kk + kb, p);
            for (mtl::size_t jj = 0; jj < m; jj += jb) {
                mtl::size_t jEnd = std::min(jj + jb, m);
                for (auto i = iStart; i < iEnd; ++i) {
                    for (auto k = kk; k < kEnd; ++k) {
                        typename MR::ValueType sum{};
                        for (auto j = jj; j < jEnd; ++j) {
                            sum += a.at(i, j) * b.at(j, k);
                        }
                        r.at(i, k) += sum;
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
        pool.addTask([&r, &a, &b, i, partSize]() { detail::parallelBlockMulTask(r, a, b, i, i + partSize); });
    }
}

CREATE_BASE_MUL_WRAPPERS(parallelBlockMul)

#endif

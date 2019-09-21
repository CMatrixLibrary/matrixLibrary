#ifndef BLOCK_MUL_H
#define BLOCK_MUL_H

#include "MatrixInterface.h"
#include "baseMulUtility.h"

// algorithm
template<typename MR, typename M1, typename M2>
void blockMul(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) { 
    constexpr mtl::size_t ib = 128;
    constexpr mtl::size_t jb = 16;
    constexpr mtl::size_t kb = 256;

    auto n = a.rowCount();
    auto m = a.columnCount();
    auto p = b.columnCount();

    for (mtl::size_t i = 0; i < n; ++i) {
        for (mtl::size_t j = 0; j < p; ++j) {
            r.at(i, j) = typename MR::ValueType{};
        }
    }

    for (mtl::size_t ii = 0; ii < n; ii += ib) {
        mtl::size_t iEnd = std::min(ii + ib, n);
        for (mtl::size_t kk = 0; kk < p; kk += kb) {
            mtl::size_t kEnd = std::min(kk + kb, p);
            for (mtl::size_t jj = 0; jj < m; jj += jb) {
                mtl::size_t jEnd = std::min(jj + jb, m);
                for (auto i = ii; i < iEnd; ++i) {
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

CREATE_BASE_MUL_WRAPPERS(blockMul)



template<typename MR, typename M1, typename M2>
void blockMul2(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    constexpr mtl::size_t ib = 128;
    constexpr mtl::size_t jb = 16;
    constexpr mtl::size_t kb = 256;

    auto n = a.rowCount();
    auto m = a.columnCount();
    auto p = b.columnCount();

    for (mtl::size_t i = 0; i < n; ++i) {
        for (mtl::size_t j = 0; j < p; ++j) {
            r.at(i, j) = typename MR::ValueType{};
        }
    }

    auto newB = new typename MR::ValueType[m*p];
    auto newBPtr = newB;
    for (mtl::size_t kk = 0; kk < p; kk += kb) {
        mtl::size_t kEnd = std::min(kk + kb, p);
        for (mtl::size_t jj = 0; jj < m; jj += jb) {
            mtl::size_t jEnd = std::min(jj + jb, m);
            for (auto k = kk; k < kEnd; ++k) {
                for (auto j = jj; j < jEnd; ++j) {
                    *newBPtr = b.at(j, k);
                    newBPtr += 1;
                }
            }
        }
    }

    for (mtl::size_t ii = 0; ii < n; ii += ib) {
        mtl::size_t iEnd = std::min(ii + ib, n);
        for (mtl::size_t kk = 0; kk < p; kk += kb) {
            mtl::size_t kEnd = std::min(kk + kb, p);
            for (mtl::size_t jj = 0; jj < m; jj += jb) {
                mtl::size_t jEnd = std::min(jj + jb, m);
                for (auto i = ii; i < iEnd; ++i) {
                    auto aVector = a[i];
                    auto rVector = r[i];
                    auto bPtr = newB + kk * m;
                    for (auto k = kk; k < kEnd; ++k) {
                        typename MR::ValueType sum{};
                        for (auto j = jj; j < jEnd; ++j) {
                            sum += aVector[j] * *bPtr;
                        }
                        bPtr += (jEnd - jj);
                        rVector[k] += sum;
                    }
                }
            }
        }
    }
}

CREATE_BASE_MUL_WRAPPERS(blockMul2)

#endif

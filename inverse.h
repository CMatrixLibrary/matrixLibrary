#ifndef INVERSE_H
#define INVERSE_H

#include <iostream>
#include <string>
#include <numeric>
#include "MatrixInterface.h"
#include "Matrix.h"
#include "MatrixExtendedFunctions.h"
#include "naiveBasicOperations.h"
#include "avxSimd.h"
#include "avxMul.h"
#include "ThreadPool.h"

namespace detail {
    template<bool IsLower, bool IsUnit = false, typename T> void solveTriangularForIdentity(T* x, const T* M, int n, int effX, int effM) {
        int j = 0;
        for (; j <= n - avx::packedCount<T>() * 2; j += avx::packedCount<T>() * 2) {
            for (int i = IsLower ? 0 : n - 1; IsLower ? (i < n) : (i >= 0); IsLower ? ++i : --i) {
                auto startK = IsLower ? j     : i + 1;
                auto endK   = IsLower ? i - 1 : j + avx::packedCount<T>() * 2 - 1;
                auto sum1 = avx::zero<T>();
                auto sum2 = avx::zero<T>();
                for (int k = startK; k <= endK; ++k) {
                    auto lVec = avx::setAllElements(M[k + i * effM]);
                    auto bVec1 = avx::loadUnaligned(&x[j + k * effX]);
                    sum1 += lVec * bVec1;
                    auto bVec2 = avx::loadUnaligned(&x[j + avx::packedCount<T>() + k * effX]);
                    sum2 += lVec * bVec2;
                }
                if constexpr (IsUnit) {
                    avx::storeUnaligned(&x[j + i * effX], sum1);
                    avx::storeUnaligned(&x[avx::packedCount<T>() + j + i * effX], sum2);
                } else {
                    avx::storeUnaligned(&x[j + i * effX], sum1 / avx::setAllElements(M[i*(effM + 1)]));
                    avx::storeUnaligned(&x[avx::packedCount<T>() + j + i * effX], sum2 / avx::setAllElements(M[i*(effM + 1)]));
                }
                for (int jj = j; jj < j + avx::packedCount<T>() * 2; ++jj) {
                    auto temp = (i == jj) / (IsUnit ? 1 : M[i*(effM + 1)]);
                    x[jj + i * effX] = temp - x[jj + i * effX];
                }
            }
        }
        for (; j < n; ++j) {
            for (int i = IsLower ? 0 : n - 1; IsLower ? (i < n) : (i >= 0); IsLower ? ++i : --i) {
                auto startK = IsLower ? j     : i + 1;
                auto endK   = IsLower ? i - 1 : j;
                auto sum = T{};
                for (int k = startK; k <= endK; ++k) {
                    sum += M[k + i * effM] * x[j + k * effX];
                }
                x[j + i * effX] = ((i == j) - sum) / (IsUnit ? 1 : M[i*(effM + 1)]);
            }
        }
    }

    template<bool IsLower, bool IsUnit = false, typename T> void solveTriangularForIdentityParallel(T* x, const T* M, int n, int effX, int effM) {
        ThreadPool pool;
        int j = 0;
        for (; j <= n - avx::packedCount<T>() * 2; j += avx::packedCount<T>() * 2) {
            pool.addTask([x, M, n, effX, effM, j] {
                for (int i = IsLower ? 0 : n - 1; IsLower ? (i < n) : (i >= 0); IsLower ? ++i : --i) {
                    auto startK = IsLower ? j     : i + 1;
                    auto endK   = IsLower ? i - 1 : j + avx::packedCount<T>() * 2 - 1;
                    auto sum1 = avx::zero<T>();
                    auto sum2 = avx::zero<T>();
                    for (int k = startK; k <= endK; ++k) {
                        auto lVec = avx::setAllElements(M[k + i * effM]);
                        auto bVec1 = avx::loadUnaligned(&x[j + k * effX]);
                        sum1 += lVec * bVec1;
                        auto bVec2 = avx::loadUnaligned(&x[j + avx::packedCount<T>() + k * effX]);
                        sum2 += lVec * bVec2;
                    }
                    if constexpr (IsUnit) {
                        avx::storeUnaligned(&x[j + i * effX], sum1);
                        avx::storeUnaligned(&x[avx::packedCount<T>() + j + i * effX], sum2);
                    } else {
                        avx::storeUnaligned(&x[j + i * effX], sum1 / avx::setAllElements(M[i*(effM + 1)]));
                        avx::storeUnaligned(&x[avx::packedCount<T>() + j + i * effX], sum2 / avx::setAllElements(M[i*(effM + 1)]));
                    }
                    for (int jj = j; jj < j + avx::packedCount<T>() * 2; ++jj) {
                        auto temp = (i == jj) / (IsUnit ? 1 : M[i*(effM + 1)]);
                        x[jj + i * effX] = temp - x[jj + i * effX];
                    }
                }
            });
        }
        for (; j < n; ++j) {
            pool.addTask([x, M, n, effX, effM, j] {
                for (int i = IsLower ? 0 : n - 1; IsLower ? (i < n) : (i >= 0); IsLower ? ++i : --i) {
                    auto startK = IsLower ? j     : i + 1;
                    auto endK   = IsLower ? i - 1 : j;
                    auto sum = T{};
                    for (int k = startK; k <= endK; ++k) {
                        sum += M[k + i * effM] * x[j + k * effX];
                    }
                    x[j + i * effX] = ((i == j) - sum) / (IsUnit ? 1 : M[i*(effM + 1)]);
                }
            });
        }
    }


    template<typename T> void avxInPlaceMultiplySubstractForLU(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        static constexpr mtl::size_t ib = 128;
        static constexpr mtl::size_t jb = 128;
        static constexpr mtl::size_t kb = 128;

        auto newB = avx::detail::getNewB<ib, jb, kb>(b, n, m, p, effB);

        ThreadPool pool;
        for (mtl::size_t ii = 0; ii < n; ii += ib) {
            mtl::size_t iEnd = std::min(ii + ib, n);
            for (mtl::size_t kk = 0; kk < p; kk += kb) {
                mtl::size_t kEnd = std::min(kk + kb, p);
                pool.addTask([=] {
                    for (mtl::size_t jj = 0; jj < m; jj += jb) {
                        mtl::size_t jEnd = std::min(jj + jb, m);
                        auto i = ii;
                        for (; i <= iEnd - 4; i += 4) {
                            auto k = kk;
                            auto aPtr = a + i * effA;
                            auto cPtr = c + i * effC;
                            auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                            for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                                auto sum1i1 = avx::zero<T>();
                                auto sum2i1 = avx::zero<T>();
                                auto sum1i2 = avx::zero<T>();
                                auto sum2i2 = avx::zero<T>();
                                auto sum1i3 = avx::zero<T>();
                                auto sum2i3 = avx::zero<T>();
                                auto sum1i4 = avx::zero<T>();
                                auto sum2i4 = avx::zero<T>();
                                for (auto j = jj; j < jEnd; ++j) {
                                    auto bVectr1 = avx::loadUnaligned(bPtr);
                                    auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                                    auto aValue1 = avx::setAllElements(aPtr[j]);
                                    auto aValue2 = avx::setAllElements(aPtr[j + effA]);
                                    auto aValue3 = avx::setAllElements(aPtr[j + 2 * effA]);
                                    auto aValue4 = avx::setAllElements(aPtr[j + 3 * effA]);
                                    sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                                    sum2i1 = avx::fma(aValue1, bVectr2, sum2i1);
                                    sum1i2 = avx::fma(aValue2, bVectr1, sum1i2);
                                    sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                                    sum1i3 = avx::fma(aValue3, bVectr1, sum1i3);
                                    sum2i3 = avx::fma(aValue3, bVectr2, sum2i3);
                                    sum1i4 = avx::fma(aValue4, bVectr1, sum1i4);
                                    sum2i4 = avx::fma(aValue4, bVectr2, sum2i4);
                                    bPtr += 2 * avx::packedCount<T>();
                                }
                                avx::storeUnaligned(&cPtr[k], avx::loadUnaligned(&cPtr[k]) - sum1i1);
                                avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]) - sum2i1);
                                avx::storeUnaligned(&cPtr[k + effC], avx::loadUnaligned(&cPtr[k + effC]) - sum1i2);
                                avx::storeUnaligned(&cPtr[k + avx::packedCount<T>() + effC], avx::loadUnaligned(&cPtr[k + avx::packedCount<T>() + effC]) - sum2i2);
                                avx::storeUnaligned(&cPtr[k + 2 * effC], avx::loadUnaligned(&cPtr[k + 2 * effC]) - sum1i3);
                                avx::storeUnaligned(&cPtr[k + avx::packedCount<T>() + 2 * effC], avx::loadUnaligned(&cPtr[k + avx::packedCount<T>() + 2 * effC]) - sum2i3);
                                avx::storeUnaligned(&cPtr[k + 3 * effC], avx::loadUnaligned(&cPtr[k + 3 * effC]) - sum1i4);
                                avx::storeUnaligned(&cPtr[k + avx::packedCount<T>() + 3 * effC], avx::loadUnaligned(&cPtr[k + avx::packedCount<T>() + 3 * effC]) - sum2i4);
                            }
                            for (; k < kEnd; ++k) {
                                auto sum1 = T{};
                                auto sum2 = T{};
                                auto sum3 = T{};
                                auto sum4 = T{};
                                for (auto j = jj; j < jEnd; ++j) {
                                    auto bVal = b[k + j * effB];
                                    sum1 += aPtr[j] * bVal;
                                    sum2 += aPtr[j + effA] * bVal;
                                    sum3 += aPtr[j + 2 * effA] * bVal;
                                    sum4 += aPtr[j + 3 * effA] * bVal;
                                }
                                cPtr[k] -= sum1;
                                cPtr[k + effC] -= sum2;
                                cPtr[k + 2 * effC] -= sum3;
                                cPtr[k + 3 * effC] -= sum4;
                            }
                        }
                        for (; i < iEnd; ++i) {
                            auto k = kk;
                            auto aPtr = a + i * effA;
                            auto cPtr = c + i * effC;
                            auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                            for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                                auto sum1 = avx::zero<T>();
                                auto sum2 = avx::zero<T>();
                                for (auto j = jj; j < jEnd; ++j) {
                                    auto bVectr1 = avx::loadUnaligned(bPtr);
                                    auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                                    auto aValue = avx::setAllElements(aPtr[j]);
                                    sum1 = avx::fma(aValue, bVectr1, sum1);
                                    sum2 = avx::fma(aValue, bVectr2, sum2);
                                    bPtr += 2 * avx::packedCount<T>();
                                }
                                avx::storeUnaligned(&cPtr[k], avx::loadUnaligned(&cPtr[k]) - sum1);
                                avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]) - sum2);
                            }
                            for (; k < kEnd; ++k) {
                                auto sum = T{};
                                for (auto j = jj; j < jEnd; ++j) {
                                    sum += aPtr[j] * b[k + j * effB];
                                }
                                cPtr[k] -= sum;
                            }
                        }
                    }
                    });
            }
        }

        pool.completeTasksAndStop();
        delete[] newB;
    }

    template<typename MR, typename M1, typename M2>
    void InPlaceMultiplySubstractForLU(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        avxInPlaceMultiplySubstractForLU(r.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), r.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
    }


    template<typename M> void blockLUPForBlockSize(MatrixInterface<M>& remainingRows, int remIndex, ArrayView<int> P, int blockSize) {
        for (int i = 0; i < blockSize; ++i) {
            int maxRowIndex = i;
            auto max = std::abs(remainingRows.at(i, remIndex + i));
            for (int j = i + 1; j < remainingRows.rowCount(); ++j) {
                if (auto newMax = std::abs(remainingRows.at(j, remIndex + i)); max < newMax) {
                    max = newMax;
                    maxRowIndex = j;
                }
            }
            if (maxRowIndex != i) {
                std::swap(P[i], P[maxRowIndex]);
                for (int j = 0; j < remainingRows.columnCount(); ++j) {
                    std::swap(remainingRows.at(i, j), remainingRows.at(maxRowIndex, j));
                }
            }
            for (int j = i + 1; j < remainingRows.rowCount(); ++j) {
                remainingRows.at(j, remIndex + i) /= remainingRows.at(i, remIndex + i);
                for (int k = i + 1; k < blockSize; ++k) {
                    remainingRows.at(j, remIndex + k) -= remainingRows.at(i, remIndex + k) * remainingRows.at(j, remIndex + i);
                }
            }
        }
    }
    
    template<typename M1, typename M2> void solveLowerTriangularInPlace(const MatrixInterface<M1>& L, MatrixInterface<M2>& X) {
        for (int i = 0; i < L.rowCount(); ++i) {
            for (int j = 0; j < X.columnCount(); ++j) {
                typename M1::ValueType sum{};
                for (int k = 0; k < i; ++k) {
                    sum += L.at(i, k) * X.at(k, j);
                }
                X.at(i, j) -= sum;
            }
        }
    }

    template<typename M> void blockLUPImpl(MatrixInterface<M>& remainingRowss, int remIndex, ArrayView<int> P, int blockSize) {
        auto remainingRows = remainingRowss.subMatrixView(0, 0, remainingRowss.rowCount(), remainingRowss.columnCount());
        while (remainingRows.rowCount() > blockSize) {
            blockLUPForBlockSize(remainingRows, remIndex, P, blockSize);

            auto LU = remainingRows.subMatrixView(0, remIndex, blockSize, blockSize);
            auto U1 = remainingRows.subMatrixView(0, remIndex + blockSize, blockSize, remainingRows.columnCount() - remIndex - blockSize);
            solveLowerTriangularInPlace(LU, U1);

            auto L1 = remainingRows.subMatrixView(blockSize, remIndex, remainingRows.rowCount() - blockSize, blockSize);
            auto E = remainingRows.subMatrixView(blockSize, remIndex + blockSize, remainingRows.rowCount() - blockSize, remainingRows.columnCount() - remIndex - blockSize);
            InPlaceMultiplySubstractForLU(E, L1, U1);

            remIndex += blockSize;
            P = ArrayView(P.data() + blockSize, P.size() - blockSize);
            remainingRows = remainingRows.subMatrixView(blockSize, 0, remainingRows.rowCount() - blockSize, remainingRows.columnCount());
        }
        blockLUPForBlockSize(remainingRows, remIndex, P, remainingRows.rowCount());
    }
}

/*
    triangular solvers
*/
template<typename M> auto solveLowerUnitTriangularForIdentity(const MatrixInterface<M>& L) {
    auto X = L.createNew();
    detail::solveTriangularForIdentityParallel<true, true>(X.data(), L.data(), L.rowCount(), X.effectiveColumnCount(), L.effectiveColumnCount());
    return X;
}
template<typename M> auto solveUpperTriangularForIdentity(const MatrixInterface<M>& U) {
    auto X = U.createNew();
    detail::solveTriangularForIdentityParallel<false>(X.data(), U.data(), U.rowCount(), X.effectiveColumnCount(), U.effectiveColumnCount());
    return X;
}

/*
    LUP + permutation apply
*/
template<typename M> void applyLUPPermutations(MatrixInterface<M>& A, std::vector<int>& P) {
    for (int i = 0; i < A.columnCount(); ++i) {
        while (i != P[i]) {
            auto p = P[i];
            std::swap(P[i], P[p]);
            for (int j = 0; j < A.rowCount(); ++j) {
                std::swap(A.at(j, i), A.at(j, p));
            }
        }
    }
}

template<typename M> auto LUP(const MatrixInterface<M>& A) {
    auto LU = A.createNew();
    LU.copy(A);
    std::vector<int> P(LU.rowCount());
    std::iota(P.begin(), P.end(), 0);
    detail::blockLUPForBlockSize(LU, 0, ArrayView(P.data(), P.size()), P.size());
    return std::pair(LU, P);
}

template<typename M> auto blockLUP(const MatrixInterface<M>& A, int cutOff=16) {
    auto LU = A.createNew();
    LU.copy(A);
    std::vector<int> P(LU.rowCount());
    std::iota(P.begin(), P.end(), 0);
    detail::blockLUPImpl(LU, 0, ArrayView(P.data(), P.size()), cutOff);
    return std::pair(LU, P);
}


/*
    Inverse Algorithms
*/

template<typename M> auto inverse(const MatrixInterface<M>& A) {
    auto[LU, P] = LUP(A);
    auto LInverse = solveLowerUnitTriangularForIdentity(LU);
    auto UInverse = solveUpperTriangularForIdentity(LU);
    auto result = avx::parallelMul(UInverse, LInverse);
    applyLUPPermutations(result, P);
    return result;
}

template<typename M> auto blockInverse(const MatrixInterface<M>& A, int cutOff=16) {
    auto[LU, P] = blockLUP(A, cutOff);
    auto LInverse = solveLowerUnitTriangularForIdentity(LU);
    auto UInverse = solveUpperTriangularForIdentity(LU);
    auto result = avx::parallelMul(UInverse, LInverse);
    applyLUPPermutations(result, P);
    return result;
}

template<typename M> auto fastUnstableBlockInverse(const MatrixInterface<M>& AMatrix, int cutOff=64) {
    if (AMatrix.rowCount() <= cutOff) {
        return inverse(AMatrix);
    }

    auto ACopy = AMatrix.createNew();
    ACopy.copy(AMatrix);

    auto A = matrixDivideView<2, 2>(ACopy);
    auto A00Inv = fastUnstableBlockInverse(A[0][0], cutOff);
    auto Temp1 = fastUnstableBlockInverse(A[1][1] - avx::parallelMul(avx::parallelMul(A[1][0], A00Inv), A[0][1]), cutOff);
    auto Temp2 = avx::parallelMul(avx::parallelMul(A00Inv, A[0][1]), Temp1);
    
    auto Cmatrix = ACopy.createNew();
    auto C = matrixDivideView<2, 2>(Cmatrix);
    naiveAdd(C[0][0], A00Inv, avx::parallelMul(avx::parallelMul(Temp2, A[1][0]), A00Inv));
    naiveNeg(C[0][1], Temp2);
    naiveNeg(C[1][0], avx::parallelMul(avx::parallelMul(Temp1, A[1][0]), A00Inv));
    C[1][1].copy(Temp1);

    return Cmatrix;
}

#ifdef USE_BLAS
namespace lapack {
    namespace detail {
        void LUP(float* LU, int* P, int n, int effLU) {
            LAPACKE_sgetrf(LAPACK_ROW_MAJOR, n, n, LU, effLU, P);
        }

        void inverse(float* Result, int n, int effR) {
            auto P = new int[n];
            LAPACKE_sgetrf(LAPACK_ROW_MAJOR, n, n, Result, effR, P);
            LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, Result, effR, P);
            delete[] P;
        }
    }
    
    template<typename M> auto inverse(const MatrixInterface<M>& A) {
        auto Result = A.createNew();
        Result.copy(A);
        detail::inverse(Result.data(), Result.rowCount(), Result.effectiveColumnCount());
        return Result;
    }

    template<typename M> auto LUP(const MatrixInterface<M>& A) {
        auto LU = A.createNew();
        LU.copy(A);
        std::vector<int> P(LU.rowCount());
        detail::LUP(LU.data(), P.data(), P.size(), LU.effectiveColumnCount());
        return std::pair(LU, P);
    }
}
#endif

#endif

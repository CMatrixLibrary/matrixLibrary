#ifndef AVX_MUL_H
#define AVX_MUL_H
#include "avxSimd.h"
#include "MatrixInterface.h"
#include "Matrix.h"
#include "ThreadPool.h"

namespace avx::detail {

    void avxMul2(int* result, const int* a, const int* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n; ++i) {
            int j = 0;
            for (; j <= q - 8; j += 8) {
                auto sum = avx::zero<int>();
                for (int k = 0; k < m; ++k) {
                    auto aValue = avx::setAllElements(a[k + i * m]);
                    auto bVectr = avx::loadUnaligned(&b[j + k * q]);
                    sum += aValue * bVectr;
                }
                avx::storeUnaligned(&result[j + i * q], sum);
            }
            for (; j < q; ++j) {
                int sum = 0;
                for (int k = 0; k < m; ++k) {
                    sum += a[k + i * m] * b[j + k * q];
                }
                result[j + i * q] = sum;
            }
        }
#endif
    }
    void avxMul3(int* result, const int* a, const int* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n; ++i) {
            int j = 0;
            for (; j <= q - 16; j += 16) {
                auto sum1 = avx::zero<int>();
                auto sum2 = avx::zero<int>();
                for (int k = 0; k < m; ++k) {
                    auto aValue = avx::setAllElements(a[k + i * m]);
                    auto bVectr1 = avx::loadUnaligned(&b[j + k * q]);
                    sum1 += aValue * bVectr1;
                    auto bVectr2 = avx::loadUnaligned(&b[8 + j + k * q]);
                    sum2 += aValue * bVectr2;
                }
                avx::storeUnaligned(&result[j + i * q], sum1);
                avx::storeUnaligned(&result[8 + j + i * q], sum2);
            }
            for (; j < q; ++j) {
                int sum = 0;
                for (int k = 0; k < m; ++k) {
                    sum += a[k + i * m] * b[j + k * q];
                }
                result[j + i * q] = sum;
            }
        }
#endif
    }
    void avxMul4(int* result, const int* a, const int* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n; ++i) {
            int j = 0;
            for (; j <= q - 32; j += 32) {
                auto sum1 = avx::zero<int>();
                auto sum2 = avx::zero<int>();
                auto sum3 = avx::zero<int>();
                auto sum4 = avx::zero<int>();
                for (int k = 0; k < m; ++k) {
                    auto aValue = avx::setAllElements(a[k + i * m]);
                    auto bVectr1 = avx::loadUnaligned(&b[j + k * q]);
                    sum1 += aValue * bVectr1;
                    auto bVectr2 = avx::loadUnaligned(&b[8 + j + k * q]);
                    sum2 += aValue * bVectr2;
                    auto bVectr3 = avx::loadUnaligned(&b[16 + j + k * q]);
                    sum3 += aValue * bVectr3;
                    auto bVectr4 = avx::loadUnaligned(&b[24 + j + k * q]);
                    sum4 += aValue * bVectr4;
                }
                avx::storeUnaligned(&result[j + i * q], sum1);
                avx::storeUnaligned(&result[8 + j + i * q], sum2);
                avx::storeUnaligned(&result[16 + j + i * q], sum3);
                avx::storeUnaligned(&result[24 + j + i * q], sum4);
            }
            for (; j < q; ++j) {
                int sum = 0;
                for (int k = 0; k < m; ++k) {
                    sum += a[k + i * m] * b[j + k * q];
                }
                result[j + i * q] = sum;
            }
        }
#endif
    }
    void avxMul5(int* result, const int* a, const int* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n*q; ++i) {
            result[i] = 0;
        }

        int jb = 256;
        int kb = 16;

        for (int jj = 0; jj < q; jj += jb) {
            int lastjb = std::min(jj + jb, q);
            for (int kk = 0; kk < m; kk += kb) {
                int lastKb = std::min(kk + kb, m);
                for (int i = 0; i < n; ++i) {
                    int j = jj;
                    int jEnd = lastjb - 16;
                    for (; j <= jEnd; j += 16) {
                        auto sum1 = avx::loadUnaligned(&result[j + i * q]);
                        auto sum2 = avx::loadUnaligned(&result[8 + j + i * q]);
                        for (int k = kk; k < lastKb; ++k) {
                            auto aValue = avx::setAllElements(a[k + i * m]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * q]);
                            sum1 += aValue * bVectr1;
                            auto bVectr2 = avx::loadUnaligned(&b[8 + j + k * q]);
                            sum2 += aValue * bVectr2;
                        }
                        avx::storeUnaligned(&result[j + i * q], sum1);
                        avx::storeUnaligned(&result[8 + j + i * q], sum2);
                    }
                    for (; j < lastjb; ++j) {
                        int sum = 0;
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
            }
        }
#endif
    }

    void avxMul6(int* result, const int* a, const int* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n*q; ++i) {
            result[i] = 0;
        }

        int jb = 256;
        int kb = 16;

        for (int jj = 0; jj < q; jj += jb) {
            int lastjb = std::min(jj + jb, q);
            for (int kk = 0; kk < m; kk += kb) {
                int lastKb = std::min(kk + kb, m);
                int i = 0;
                for (; i <= n - 2; i += 2) {
                    int j = jj;
                    int jEnd = lastjb - 16;
                    for (; j <= jEnd; j += 16) {
                        auto sum1i1 = avx::loadUnaligned(&result[j + i * q]);
                        auto sum1i2 = avx::loadUnaligned(&result[8 + j + i * q]);
                        auto sum2i1 = avx::loadUnaligned(&result[j + (i + 1) * q]);
                        auto sum2i2 = avx::loadUnaligned(&result[8 + j + (i + 1) * q]);
                        for (int k = kk; k < lastKb; ++k) {
                            auto aValue1 = avx::setAllElements(a[k + i * m]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * q]);
                            sum1i1 += aValue1 * bVectr1;
                            auto bVectr2 = avx::loadUnaligned(&b[8 + j + k * q]);
                            sum1i2 += aValue1 * bVectr2;
                            auto aValue2 = avx::setAllElements(a[k + (i + 1) * m]);
                            sum2i1 += aValue2 * bVectr1;
                            sum2i2 += aValue2 * bVectr2;
                        }
                        avx::storeUnaligned(&result[j + i * q], sum1i1);
                        avx::storeUnaligned(&result[8 + j + i * q], sum1i2);
                        avx::storeUnaligned(&result[j + (i + 1) * q], sum1i1);
                        avx::storeUnaligned(&result[8 + j + (i + 1) * q], sum1i2);
                    }
                    for (; j < lastjb; ++j) {
                        int sum = 0;
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
                if (i < n) {
                    int j = jj;
                    int jEnd = lastjb - 16;
                    for (; j <= jEnd; j += 16) {
                        auto sum1 = avx::loadUnaligned(&result[j + i * q]);
                        auto sum2 = avx::loadUnaligned(&result[8 + j + i * q]);
                        for (int k = kk; k < lastKb; ++k) {
                            auto aValue = avx::setAllElements(a[k + i * m]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * q]);
                            sum1 += aValue * bVectr1;
                            auto bVectr2 = avx::loadUnaligned(&b[8 + j + k * q]);
                            sum2 += aValue * bVectr2;
                        }
                        avx::storeUnaligned(&result[j + i * q], sum1);
                        avx::storeUnaligned(&result[8 + j + i * q], sum2);
                    }
                    for (; j < lastjb; ++j) {
                        int sum = 0;
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
            }
        }
#endif
    }

    template<typename T> void avxMul8(T* result, const T* a, const T* b, int n) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n*n; ++i) {
            result[i] = T{};
        }

        int ib = std::min(256, n);
        int jb = std::min(512, n);
        int kb = std::min(16, n);

        T* mat2 = avx::alignedArrayAlloc<T>(n*n);
        int m2idx = 0;
        for (int jj = 0; jj < n; jj += jb) {
            for (int kk = 0; kk < n; kk += kb) {
                for (int j = jj; j < jj + jb; j += 2 * avx::packedCount<T>()) {
                    for (int k = kk; k < kk + kb; k++) {
                        auto vecA_mat2 = avx::loadUnaligned(&b[j + k * n]);
                        auto vecB_mat2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * n]);
                        avx::storeUnaligned(&mat2[m2idx], vecA_mat2);
                        avx::storeUnaligned(&mat2[m2idx + avx::packedCount<T>()], vecB_mat2);
                        m2idx += 2 * avx::packedCount<T>();
                    }
                }
            }
        }

        for (int ii = 0; ii < n; ii += ib) {
            for (int jj = 0; jj < n; jj += jb) {
                for (int kk = 0; kk < n; kk += kb) {
                    for (int i = ii; i < ii + ib; i += 2) {
                        for (int j = jj; j < jj + jb; j += 2 * avx::packedCount<T>()) {
                            int m2idx = (j - jj) * kb + kk * jb + jj * n;
                            auto sumA_1 = avx::loadUnaligned(&result[i*n + j]);
                            auto sumB_1 = avx::loadUnaligned(&result[i*n + j + avx::packedCount<T>()]);
                            auto sumA_2 = avx::loadUnaligned(&result[(i + 1)*n + j]);
                            auto sumB_2 = avx::loadUnaligned(&result[(i + 1)*n + j + avx::packedCount<T>()]);
                            for (int k = kk; k < kk + kb; k++) {
                                auto bc_mat1_1 = avx::setAllElements(a[i*n + k]);
                                auto vecA_mat2 = avx::loadUnaligned(&mat2[m2idx]);
                                auto vecB_mat2 = avx::loadUnaligned(&mat2[m2idx + avx::packedCount<T>()]);
                                sumA_1 += bc_mat1_1 * vecA_mat2;
                                sumB_1 += bc_mat1_1 * vecB_mat2;
                                auto bc_mat1_2 = avx::setAllElements(a[(i + 1)*n + k]);
                                sumA_2 += bc_mat1_2 * vecA_mat2;
                                sumB_2 += bc_mat1_2 * vecB_mat2;
                                m2idx += 2 * avx::packedCount<T>();
                            }
                            avx::storeUnaligned(&result[i*n + j], sumA_1);
                            avx::storeUnaligned(&result[i*n + j + avx::packedCount<T>()], sumB_1);
                            avx::storeUnaligned(&result[(i + 1)*n + j], sumA_2);
                            avx::storeUnaligned(&result[(i + 1) * n + j + avx::packedCount<T>()], sumB_2);
                        }
                    }
                }
            }
        }
        avx::alignedArrayDealloc(mat2);
#endif
    }

    template<typename T> void mul7Task(T* result, const T* a, const T* b, int ii, int iEnd, int n, int m, int q) {
        constexpr int jb = 512;
        constexpr int kb = 16;

        for (int jj = 0; jj < q; jj += jb) {
            int jEnd = std::min(jj + jb, q);
            for (int kk = 0; kk < m; kk += kb) {
                int kEnd = std::min(kk + kb, m);
                int i = ii;
                for (; i <= iEnd - 2; i += 2) {
                    int j = jj;
                    for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                        auto sum1i1 = avx::loadUnaligned(&result[j + i * q]);
                        auto sum1i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * q]);
                        auto sum2i1 = avx::loadUnaligned(&result[j + (i + 1) * q]);
                        auto sum2i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * q]);
                        for (int k = kk; k < kEnd; ++k) {
                            auto aValue1 = avx::setAllElements(a[k + i * m]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * q]);
                            sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                            auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * q]);
                            sum1i2 = avx::fma(aValue1, bVectr2, sum1i2);
                            auto aValue2 = avx::setAllElements(a[k + (i + 1) * m]);
                            sum2i1 = avx::fma(aValue2, bVectr1, sum2i1);
                            sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                        }
                        avx::storeUnaligned(&result[j + i * q], sum1i1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * q], sum1i2);
                        avx::storeUnaligned(&result[j + (i + 1) * q], sum2i1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * q], sum2i2);
                    }
                    for (; j < jEnd; ++j) {
                        T& sum1 = result[j + i * q];
                        for (int k = kk; k < kEnd; ++k) {
                            sum1 += a[k + i * m] * b[j + k * q];
                        }
                        T& sum2 = result[j + (i + 1) * q];
                        for (int k = kk; k < kEnd; ++k) {
                            sum2 += a[k + (i + 1) * m] * b[j + k * q];
                        }
                    }
                }
                if (i < iEnd) {
                    int j = jj;
                    for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                        auto sum1 = avx::loadUnaligned(&result[j + i * q]);
                        auto sum2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * q]);
                        for (int k = kk; k < kEnd; ++k) {
                            auto aValue = avx::setAllElements(a[k + i * m]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * q]);
                            sum1 = avx::fma(aValue, bVectr1, sum1);
                            auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * q]);
                            sum2 = avx::fma(aValue, bVectr2, sum2);
                        }
                        avx::storeUnaligned(&result[j + i * q], sum1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * q], sum2);
                    }
                    for (; j < jEnd; ++j) {
                        T& sum = result[j + i * q];
                        for (int k = kk; k < kEnd; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                    }
                }
            }
        }
    }
    template<typename T> void parallelMul7(T* result, const T* a, const T* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n*q; ++i) {
            result[i] = T{};
        }

        constexpr int ib = 256;

        ThreadPool pool;
        for (int ii = 0; ii < n; ii += ib) {
            int iEnd = std::min(ii + ib, n);
            pool.addTask([=]() { mul7Task(result, a, b, ii, iEnd, n, m, q); });
        }
#endif
    }
    template<typename T> void mul7(T* result, const T* a, const T* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n*q; ++i) {
            result[i] = T{};
        }

        constexpr int ib = 256;

        for (int ii = 0; ii < n; ii += ib) {
            int iEnd = std::min(ii + ib, n);
            mul7Task(result, a, b, ii, iEnd, n, m, q);
        }
#endif
    }


    template<typename T> void mul7Task(T* result, const T* a, const T* b, int ii, int iEnd, int n, int m, int q, int effR, int effA, int effB) {
        constexpr int jb = 512;
        constexpr int kb = 16;

        for (int jj = 0; jj < q; jj += jb) {
            int jEnd = std::min(jj + jb, q);
            for (int kk = 0; kk < m; kk += kb) {
                int kEnd = std::min(kk + kb, m);
                int i = ii;
                for (; i <= iEnd - 2; i += 2) {
                    int j = jj;
                    for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                        auto sum1i1 = avx::loadUnaligned(&result[j + i * effR]);
                        auto sum1i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * effR]);
                        auto sum2i1 = avx::loadUnaligned(&result[j + (i + 1) * effR]);
                        auto sum2i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * effR]);
                        for (int k = kk; k < kEnd; ++k) {
                            auto aValue1 = avx::setAllElements(a[k + i * effA]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * effB]);
                            sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                            auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * effB]);
                            sum1i2 = avx::fma(aValue1, bVectr2, sum1i2);
                            auto aValue2 = avx::setAllElements(a[k + (i + 1) * effA]);
                            sum2i1 = avx::fma(aValue2, bVectr1, sum2i1);
                            sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                        }
                        avx::storeUnaligned(&result[j + i * effR], sum1i1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * effR], sum1i2);
                        avx::storeUnaligned(&result[j + (i + 1) * effR], sum2i1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * effR], sum2i2);
                    }
                    for (; j < jEnd; ++j) {
                        T& sum1 = result[j + i * effR];
                        for (int k = kk; k < kEnd; ++k) {
                            sum1 += a[k + i * effA] * b[j + k * effB];
                        }
                        T& sum2 = result[j + (i + 1) * effR];
                        for (int k = kk; k < kEnd; ++k) {
                            sum2 += a[k + (i + 1) * effA] * b[j + k * effB];
                        }
                    }
                }
                if (i < iEnd) {
                    int j = jj;
                    for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                        auto sum1 = avx::loadUnaligned(&result[j + i * effR]);
                        auto sum2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * effR]);
                        for (int k = kk; k < kEnd; ++k) {
                            auto aValue = avx::setAllElements(a[k + i * effA]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * effB]);
                            sum1 = avx::fma(aValue, bVectr1, sum1);
                            auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * effB]);
                            sum2 = avx::fma(aValue, bVectr2, sum2);
                        }
                        avx::storeUnaligned(&result[j + i * effR], sum1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * effR], sum2);
                    }
                    for (; j < jEnd; ++j) {
                        T& sum = result[j + i * effR];
                        for (int k = kk; k < kEnd; ++k) {
                            sum += a[k + i * effA] * b[j + k * effB];
                        }
                    }
                }
            }
        }
    }


    template<typename T> void parallelMul7(T* result, const T* a, const T* b, int n, int m, int q, int effR, int effA, int effB) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < q; ++j) {
                result[j + i * effR] = T{};
            }
        }

        constexpr int ib = 256;

        ThreadPool pool;
        for (int ii = 0; ii < n; ii += ib) {
            int iEnd = std::min(ii + ib, n);
            pool.addTask([=]() { mul7Task(result, a, b, ii, iEnd, n, m, q, effR, effA, effB); });
        }
#endif
    }
    template<typename T> void mul7(T* result, const T* a, const T* b, int n, int m, int q, int effR, int effA, int effB) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < q; ++j) {
                result[j + i * effR] = T{};
            }
        }

        constexpr int ib = 256;

        for (int ii = 0; ii < n; ii += ib) {
            int iEnd = std::min(ii + ib, n);
            mul7Task(result, a, b, ii, iEnd, n, m, q, effR, effA, effB);
        }
#endif
    }

    template<int n, int m, int q, int effR, int effA, int effB, typename T> void mul7(T* result, const T* a, const T* b) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < q; ++j) {
                result[j + i * effR] = T{};
            }
        }

        constexpr int ib = 256;
        constexpr int jb = 512;
        constexpr int kb = 16;

        for (int ii = 0; ii < n; ii += ib) {
            int iEnd = std::min(ii + ib, n);
            for (int jj = 0; jj < q; jj += jb) {
                int jEnd = std::min(jj + jb, q);
                for (int kk = 0; kk < m; kk += kb) {
                    int kEnd = std::min(kk + kb, m);
                    int i = ii;
                    for (; i <= iEnd - 2; i += 2) {
                        int j = jj;
                        for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                            auto sum1i1 = avx::loadUnaligned(&result[j + i * effR]);
                            auto sum1i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * effR]);
                            auto sum2i1 = avx::loadUnaligned(&result[j + (i + 1) * effR]);
                            auto sum2i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * effR]);
                            for (int k = kk; k < kEnd; ++k) {
                                auto aValue1 = avx::setAllElements(a[k + i * effA]);
                                auto bVectr1 = avx::loadUnaligned(&b[j + k * effB]);
                                sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                                auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * effB]);
                                sum1i2 = avx::fma(aValue1, bVectr2, sum1i2);
                                auto aValue2 = avx::setAllElements(a[k + (i + 1) * effA]);
                                sum2i1 = avx::fma(aValue2, bVectr1, sum2i1);
                                sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                            }
                            avx::storeUnaligned(&result[j + i * effR], sum1i1);
                            avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * effR], sum1i2);
                            avx::storeUnaligned(&result[j + (i + 1) * effR], sum2i1);
                            avx::storeUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * effR], sum2i2);
                        }
                        for (; j < jEnd; ++j) {
                            T& sum1 = result[j + i * effR];
                            for (int k = kk; k < kEnd; ++k) {
                                sum1 += a[k + i * effA] * b[j + k * effB];
                            }
                            T& sum2 = result[j + (i + 1) * effR];
                            for (int k = kk; k < kEnd; ++k) {
                                sum2 += a[k + (i + 1) * effA] * b[j + k * effB];
                            }
                        }
                    }
                    if (i < iEnd) {
                        int j = jj;
                        for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                            auto sum1 = avx::loadUnaligned(&result[j + i * effR]);
                            auto sum2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * effR]);
                            for (int k = kk; k < kEnd; ++k) {
                                auto aValue = avx::setAllElements(a[k + i * effA]);
                                auto bVectr1 = avx::loadUnaligned(&b[j + k * effB]);
                                sum1 = avx::fma(aValue, bVectr1, sum1);
                                auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * effB]);
                                sum2 = avx::fma(aValue, bVectr2, sum2);
                            }
                            avx::storeUnaligned(&result[j + i * effR], sum1);
                            avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * effR], sum2);
                        }
                        for (; j < jEnd; ++j) {
                            T& sum = result[j + i * effR];
                            for (int k = kk; k < kEnd; ++k) {
                                sum += a[k + i * effA] * b[j + k * effB];
                            }
                        }
                    }
                }
            }
        }
#endif
    }

    template<int n, int m, int q, int effR, int effA, int effB, typename T> void parallelMul7Task(T* result, const T* a, const T* b, int ii, int iEnd) {
        constexpr int jb = 512;
        constexpr int kb = 16;

        for (int jj = 0; jj < q; jj += jb) {
            int jEnd = std::min(jj + jb, q);
            for (int kk = 0; kk < m; kk += kb) {
                int kEnd = std::min(kk + kb, m);
                int i = ii;
                for (; i <= iEnd - 2; i += 2) {
                    int j = jj;
                    for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                        auto sum1i1 = avx::loadUnaligned(&result[j + i * effR]);
                        auto sum1i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * effR]);
                        auto sum2i1 = avx::loadUnaligned(&result[j + (i + 1) * effR]);
                        auto sum2i2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * effR]);
                        for (int k = kk; k < kEnd; ++k) {
                            auto aValue1 = avx::setAllElements(a[k + i * effA]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * effB]);
                            sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                            auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * effB]);
                            sum1i2 = avx::fma(aValue1, bVectr2, sum1i2);
                            auto aValue2 = avx::setAllElements(a[k + (i + 1) * effA]);
                            sum2i1 = avx::fma(aValue2, bVectr1, sum2i1);
                            sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                        }
                        avx::storeUnaligned(&result[j + i * effR], sum1i1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * effR], sum1i2);
                        avx::storeUnaligned(&result[j + (i + 1) * effR], sum2i1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + (i + 1) * effR], sum2i2);
                    }
                    for (; j < jEnd; ++j) {
                        T& sum1 = result[j + i * effR];
                        for (int k = kk; k < kEnd; ++k) {
                            sum1 += a[k + i * effA] * b[j + k * effB];
                        }
                        T& sum2 = result[j + (i + 1) * effR];
                        for (int k = kk; k < kEnd; ++k) {
                            sum2 += a[k + (i + 1) * effA] * b[j + k * effB];
                        }
                    }
                }
                if (i < iEnd) {
                    int j = jj;
                    for (; j <= jEnd - 2 * avx::packedCount<T>(); j += 2 * avx::packedCount<T>()) {
                        auto sum1 = avx::loadUnaligned(&result[j + i * effR]);
                        auto sum2 = avx::loadUnaligned(&result[avx::packedCount<T>() + j + i * effR]);
                        for (int k = kk; k < kEnd; ++k) {
                            auto aValue = avx::setAllElements(a[k + i * effA]);
                            auto bVectr1 = avx::loadUnaligned(&b[j + k * effB]);
                            sum1 = avx::fma(aValue, bVectr1, sum1);
                            auto bVectr2 = avx::loadUnaligned(&b[avx::packedCount<T>() + j + k * effB]);
                            sum2 = avx::fma(aValue, bVectr2, sum2);
                        }
                        avx::storeUnaligned(&result[j + i * effR], sum1);
                        avx::storeUnaligned(&result[avx::packedCount<T>() + j + i * effR], sum2);
                    }
                    for (; j < jEnd; ++j) {
                        T& sum = result[j + i * effR];
                        for (int k = kk; k < kEnd; ++k) {
                            sum += a[k + i * effA] * b[j + k * effB];
                        }
                    }
                }
            }
        }
    }

    template<int n, int m, int q, int effR, int effA, int effB, typename T> void parallelMul7(T* result, const T* a, const T* b) {
#ifdef AVX2_IS_AVAILABLE
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < q; ++j) {
                result[j + i * effR] = T{};
            }
        }

        constexpr int ib = 256;

        ThreadPool pool;
        for (int ii = 0; ii < n; ii += ib) {
            int iEnd = std::min(ii + ib, n);
            pool.addTask([=]() { parallelMul7Task<n, m, q, effR, effA, effB>(result, a, b, ii, iEnd); });
        }
#endif
    }

    

    template<int jb, int kb, typename T>
    void mul9Task(T* c, const T* a, const T* b, const T* newB, int ii, int iEnd, int n, int m, int p, int effC, int effA, int effB) {
        for (mtl::size_t kk = 0; kk < p; kk += kb) {
            mtl::size_t kEnd = std::min(kk + kb, p);
            for (mtl::size_t jj = 0; jj < m; jj += jb) {
                mtl::size_t jEnd = std::min(jj + jb, m);
                auto i = ii;
                for (; i <= iEnd - 2; i += 2) {
                    auto k = kk;
                    auto aPtr = a + i * effA;
                    auto cPtr = c + i * effC;
                    auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                    for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                        auto sum1i1 = avx::loadUnaligned(&cPtr[k]);
                        auto sum2i1 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                        auto sum1i2 = avx::loadUnaligned(&cPtr[k + effC]);
                        auto sum2i2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>() + effC]);
                        for (auto j = jj; j < jEnd; ++j) {
                            auto bVectr1 = avx::loadUnaligned(bPtr);
                            auto aValue1 = avx::setAllElements(aPtr[j]);
                            sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                            auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                            sum2i1 = avx::fma(aValue1, bVectr2, sum2i1);
                            auto aValue2 = avx::setAllElements(aPtr[j + effA]);
                            sum1i2 = avx::fma(aValue2, bVectr1, sum1i2);
                            sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                            bPtr += 2 * avx::packedCount<T>();
                        }
                        avx::storeUnaligned(&cPtr[k], sum1i1);
                        avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2i1);
                        avx::storeUnaligned(&cPtr[k + effC], sum1i2);
                        avx::storeUnaligned(&cPtr[k + avx::packedCount<T>() + effC], sum2i2);
                    }
                    for (; k < kEnd; ++k) {
                        auto sum1 = T{};
                        auto sum2 = T{};
                        for (auto j = jj; j < jEnd; ++j) {
                            sum1 += aPtr[j] * b[k + j * effB];
                            sum2 += aPtr[j + effA] * b[k + j * effB];
                        }
                        cPtr[k] += sum1;
                        cPtr[k + effC] += sum2;
                    }
                }
                if (i < iEnd) {
                    auto k = kk;
                    auto aPtr = a + i * effA;
                    auto cPtr = c + i * effC;
                    auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                    for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                        auto sum1 = avx::loadUnaligned(&cPtr[k]);
                        auto sum2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                        for (auto j = jj; j < jEnd; ++j) {
                            auto aValue = avx::setAllElements(aPtr[j]);
                            auto bVectr1 = avx::loadUnaligned(bPtr);
                            sum1 += aValue * bVectr1;
                            auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                            sum2 += aValue * bVectr2;
                            bPtr += 2 * avx::packedCount<T>();
                        }
                        avx::storeUnaligned(&cPtr[k], sum1);
                        avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2);
                    }
                    for (; k < kEnd; ++k) {
                        auto sum = T{};
                        for (auto j = jj; j < jEnd; ++j) {
                            sum += aPtr[j] * b[k + j * effB];
                        }
                        cPtr[k] += sum;
                    }
                }
            }
        }
    }

    template<int ib, int jb, int kb, typename T>
    T* getNewB(const T* b, int n, int m, int p, int effB) {
        auto newB = new T[m*p];
        auto newBPtr = newB;
        for (mtl::size_t kk = 0; kk < p; kk += kb) {
            mtl::size_t kEnd = std::min(kk + kb, p);
            for (mtl::size_t jj = 0; jj < m; jj += jb) {
                mtl::size_t jEnd = std::min(jj + jb, m);
                for (auto k = kk; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                    auto bPtr = b + k + jj * effB;
                    for (auto j = jj; j < jEnd; ++j) {
                        avx::storeUnaligned(newBPtr, avx::loadUnaligned(bPtr));
                        avx::storeUnaligned(newBPtr + avx::packedCount<T>(), avx::loadUnaligned(bPtr + avx::packedCount<T>()));
                        bPtr += effB;
                        newBPtr += 2 * avx::packedCount<T>();
                    }
                }
            }
        }
        return newB;
    }

    template<typename T>
    void mul9(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        static constexpr mtl::size_t ib = 128;
        static constexpr mtl::size_t jb = 16;
        static constexpr mtl::size_t kb = 256;

        for (mtl::size_t i = 0; i < n; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                c[j + i * effC] = T{};
            }
        }

        auto newB = getNewB<ib, jb, kb>(b, n, m, p, effB);

        for (mtl::size_t ii = 0; ii < n; ii += ib) {
            mtl::size_t iEnd = std::min(ii + ib, n);
            mul9Task<jb, kb>(c, a, b, newB, ii, iEnd, n, m, p, effC, effA, effB);
        }

        delete[] newB;
    }

    template<typename T>
    void parallelMul9(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        static constexpr mtl::size_t ib = 128;
        static constexpr mtl::size_t jb = 16;
        static constexpr mtl::size_t kb = 256;

        for (mtl::size_t i = 0; i < n; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                c[j + i * effC] = T{};
            }
        }

        auto newB = getNewB<ib, jb, kb>(b, n, m, p, effB);

        ThreadPool pool;
        for (mtl::size_t ii = 0; ii < n; ii += ib) {
            mtl::size_t iEnd = std::min(ii + ib, n);
            pool.addTask([=]() {mul9Task<jb, kb>(c, a, b, newB, ii, iEnd, n, m, p, effC, effA, effB); });
        }
        pool.completeTasksAndStop();

        delete[] newB;
    }



    template<int jb, int kb, int n, int m, int p, int effC, int effA, int effB, typename T>
    void mul9Task(T* c, const T* a, const T* b, const T* newB, int ii, int iEnd) {
        for (mtl::size_t kk = 0; kk < p; kk += kb) {
            mtl::size_t kEnd = std::min(kk + kb, p);
            for (mtl::size_t jj = 0; jj < m; jj += jb) {
                mtl::size_t jEnd = std::min(jj + jb, m);
                auto i = ii;
                for (; i <= iEnd - 2; i += 2) {
                    auto k = kk;
                    auto aPtr = a + i * effA;
                    auto cPtr = c + i * effC;
                    auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                    for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                        auto sum1i1 = avx::loadUnaligned(&cPtr[k]);
                        auto sum2i1 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                        auto sum1i2 = avx::loadUnaligned(&cPtr[k + effC]);
                        auto sum2i2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>() + effC]);
                        for (auto j = jj; j < jEnd; ++j) {
                            auto bVectr1 = avx::loadUnaligned(bPtr);
                            auto aValue1 = avx::setAllElements(aPtr[j]);
                            sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                            auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                            sum2i1 = avx::fma(aValue1, bVectr2, sum2i1);
                            auto aValue2 = avx::setAllElements(aPtr[j + effA]);
                            sum1i2 = avx::fma(aValue2, bVectr1, sum1i2);
                            sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                            bPtr += 2 * avx::packedCount<T>();
                        }
                        avx::storeUnaligned(&cPtr[k], sum1i1);
                        avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2i1);
                        avx::storeUnaligned(&cPtr[k + effC], sum1i2);
                        avx::storeUnaligned(&cPtr[k + avx::packedCount<T>() + effC], sum2i2);
                    }
                    for (; k < kEnd; ++k) {
                        auto sum1 = T{};
                        auto sum2 = T{};
                        for (auto j = jj; j < jEnd; ++j) {
                            sum1 += aPtr[j] * b[k + j * effB];
                            sum2 += aPtr[j + effA] * b[k + j * effB];
                        }
                        cPtr[k] += sum1;
                        cPtr[k + effC] += sum2;
                    }
                }
                if (i < iEnd) {
                    auto k = kk;
                    auto aPtr = a + i * effA;
                    auto cPtr = c + i * effC;
                    auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                    for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                        auto sum1 = avx::loadUnaligned(&cPtr[k]);
                        auto sum2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                        for (auto j = jj; j < jEnd; ++j) {
                            auto aValue = avx::setAllElements(aPtr[j]);
                            auto bVectr1 = avx::loadUnaligned(bPtr);
                            sum1 += aValue * bVectr1;
                            auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                            sum2 += aValue * bVectr2;
                            bPtr += 2 * avx::packedCount<T>();
                        }
                        avx::storeUnaligned(&cPtr[k], sum1);
                        avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2);
                    }
                    for (; k < kEnd; ++k) {
                        auto sum = T{};
                        for (auto j = jj; j < jEnd; ++j) {
                            sum += aPtr[j] * b[k + j * effB];
                        }
                        cPtr[k] += sum;
                    }
                }
            }
        }
    }
    template<int n, int m, int p, int effC, int effA, int effB, typename T>
    void mul9(T* c, const T* a, const T* b) {
        static constexpr mtl::size_t ib = 128;
        static constexpr mtl::size_t jb = 16;
        static constexpr mtl::size_t kb = 256;

        for (mtl::size_t i = 0; i < n; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                c[j + i * effC] = T{};
            }
        }

        auto newB = getNewB<ib, jb, kb>(b, n, m, p, effB);

        for (mtl::size_t ii = 0; ii < n; ii += ib) {
            mtl::size_t iEnd = std::min(ii + ib, n);
            mul9Task<jb, kb, n, m, p, effC, effA, effB>(c, a, b, newB, ii, iEnd);
        }

        delete[] newB;
    }

    template<int n, int m, int p, int effC, int effA, int effB, typename T>
    void parallelMul9(T* c, const T* a, const T* b) {
        static constexpr mtl::size_t ib = 128;
        static constexpr mtl::size_t jb = 16;
        static constexpr mtl::size_t kb = 256;

        for (mtl::size_t i = 0; i < n; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                c[j + i * effC] = T{};
            }
        }

        auto newB = getNewB<ib, jb, kb>(b, n, m, p, effB);

        ThreadPool pool;
        for (mtl::size_t ii = 0; ii < n; ii += ib) {
            mtl::size_t iEnd = std::min(ii + ib, n);
            pool.addTask([=]() {mul9Task<jb, kb, n, m, p, effC, effA, effB>(c, a, b, newB, ii, iEnd); });
        }
        pool.completeTasksAndStop();

        delete[] newB;
    }






    template<int jb, typename T>
    void mul10Task(T* c, const T* a, const T* b, const T* newB, int ii, int iEnd, int kk, int kEnd, int n, int m, int p, int effC, int effA, int effB) {        
        for (mtl::size_t jj = 0; jj < m; jj += jb) {
            mtl::size_t jEnd = std::min(jj + jb, m);
            auto i = ii;
            for (; i <= iEnd - 2; i += 2) {
                auto k = kk;
                auto aPtr = a + i * effA;
                auto cPtr = c + i * effC;
                auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                    auto sum1i1 = avx::loadUnaligned(&cPtr[k]);
                    auto sum2i1 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                    auto sum1i2 = avx::loadUnaligned(&cPtr[k + effC]);
                    auto sum2i2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>() + effC]);
                    for (auto j = jj; j < jEnd; ++j) {
                        auto bVectr1 = avx::loadUnaligned(bPtr);
                        auto aValue1 = avx::setAllElements(aPtr[j]);
                        sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                        auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                        sum2i1 = avx::fma(aValue1, bVectr2, sum2i1);
                        auto aValue2 = avx::setAllElements(aPtr[j + effA]);
                        sum1i2 = avx::fma(aValue2, bVectr1, sum1i2);
                        sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                        bPtr += 2 * avx::packedCount<T>();
                    }
                    avx::storeUnaligned(&cPtr[k], sum1i1);
                    avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2i1);
                    avx::storeUnaligned(&cPtr[k + effC], sum1i2);
                    avx::storeUnaligned(&cPtr[k + avx::packedCount<T>() + effC], sum2i2);
                }
                for (; k < kEnd; ++k) {
                    auto sum1 = T{};
                    auto sum2 = T{};
                    for (auto j = jj; j < jEnd; ++j) {
                        sum1 += aPtr[j] * b[k + j * effB];
                        sum2 += aPtr[j + effA] * b[k + j * effB];
                    }
                    cPtr[k] += sum1;
                    cPtr[k + effC] += sum2;
                }
            }
            if (i < iEnd) {
                auto k = kk;
                auto aPtr = a + i * effA;
                auto cPtr = c + i * effC;
                auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                    auto sum1 = avx::loadUnaligned(&cPtr[k]);
                    auto sum2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                    for (auto j = jj; j < jEnd; ++j) {
                        auto aValue = avx::setAllElements(aPtr[j]);
                        auto bVectr1 = avx::loadUnaligned(bPtr);
                        sum1 += aValue * bVectr1;
                        auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                        sum2 += aValue * bVectr2;
                        bPtr += 2 * avx::packedCount<T>();
                    }
                    avx::storeUnaligned(&cPtr[k], sum1);
                    avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2);
                }
                for (; k < kEnd; ++k) {
                    auto sum = T{};
                    for (auto j = jj; j < jEnd; ++j) {
                        sum += aPtr[j] * b[k + j * effB];
                    }
                    cPtr[k] += sum;
                }
            }
        }
    }
    template<typename T>
    void parallelMul10(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        static constexpr mtl::size_t ib = 128;
        static constexpr mtl::size_t jb = 16;
        static constexpr mtl::size_t kb = 256;

        for (mtl::size_t i = 0; i < n; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                c[j + i * effC] = T{};
            }
        }

        auto newB = getNewB<ib, jb, kb>(b, n, m, p, effB);

        ThreadPool pool;
        for (mtl::size_t ii = 0; ii < n; ii += ib) {
            mtl::size_t iEnd = std::min(ii + ib, n);
            for (mtl::size_t kk = 0; kk < p; kk += kb) {
                mtl::size_t kEnd = std::min(kk + kb, p);
                pool.addTask([=]() {mul10Task<jb>(c, a, b, newB, ii, iEnd, kk, kEnd, n, m, p, effC, effA, effB); });
            }
        }
        pool.completeTasksAndStop();

        delete[] newB;
    }



    template<int jb, int n, int m, int p, int effC, int effA, int effB, typename T>
    void mul10Task(T* c, const T* a, const T* b, const T* newB, int ii, int iEnd, int kk, int kEnd) {
        for (mtl::size_t jj = 0; jj < m; jj += jb) {
            mtl::size_t jEnd = std::min(jj + jb, m);
            auto i = ii;
            for (; i <= iEnd - 2; i += 2) {
                auto k = kk;
                auto aPtr = a + i * effA;
                auto cPtr = c + i * effC;
                auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                    auto sum1i1 = avx::loadUnaligned(&cPtr[k]);
                    auto sum2i1 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                    auto sum1i2 = avx::loadUnaligned(&cPtr[k + effC]);
                    auto sum2i2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>() + effC]);
                    for (auto j = jj; j < jEnd; ++j) {
                        auto bVectr1 = avx::loadUnaligned(bPtr);
                        auto aValue1 = avx::setAllElements(aPtr[j]);
                        sum1i1 = avx::fma(aValue1, bVectr1, sum1i1);
                        auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                        sum2i1 = avx::fma(aValue1, bVectr2, sum2i1);
                        auto aValue2 = avx::setAllElements(aPtr[j + effA]);
                        sum1i2 = avx::fma(aValue2, bVectr1, sum1i2);
                        sum2i2 = avx::fma(aValue2, bVectr2, sum2i2);
                        bPtr += 2 * avx::packedCount<T>();
                    }
                    avx::storeUnaligned(&cPtr[k], sum1i1);
                    avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2i1);
                    avx::storeUnaligned(&cPtr[k + effC], sum1i2);
                    avx::storeUnaligned(&cPtr[k + avx::packedCount<T>() + effC], sum2i2);
                }
                for (; k < kEnd; ++k) {
                    auto sum1 = T{};
                    auto sum2 = T{};
                    for (auto j = jj; j < jEnd; ++j) {
                        sum1 += aPtr[j] * b[k + j * effB];
                        sum2 += aPtr[j + effA] * b[k + j * effB];
                    }
                    cPtr[k] += sum1;
                    cPtr[k + effC] += sum2;
                }
            }
            if (i < iEnd) {
                auto k = kk;
                auto aPtr = a + i * effA;
                auto cPtr = c + i * effC;
                auto bPtr = newB + kk * m + jj * ((kEnd - kk) - (kEnd - kk) % (2 * avx::packedCount<T>()));
                for (; k <= kEnd - 2 * avx::packedCount<T>(); k += 2 * avx::packedCount<T>()) {
                    auto sum1 = avx::loadUnaligned(&cPtr[k]);
                    auto sum2 = avx::loadUnaligned(&cPtr[k + avx::packedCount<T>()]);
                    for (auto j = jj; j < jEnd; ++j) {
                        auto aValue = avx::setAllElements(aPtr[j]);
                        auto bVectr1 = avx::loadUnaligned(bPtr);
                        sum1 += aValue * bVectr1;
                        auto bVectr2 = avx::loadUnaligned(bPtr + avx::packedCount<T>());
                        sum2 += aValue * bVectr2;
                        bPtr += 2 * avx::packedCount<T>();
                    }
                    avx::storeUnaligned(&cPtr[k], sum1);
                    avx::storeUnaligned(&cPtr[k + avx::packedCount<T>()], sum2);
                }
                for (; k < kEnd; ++k) {
                    auto sum = T{};
                    for (auto j = jj; j < jEnd; ++j) {
                        sum += aPtr[j] * b[k + j * effB];
                    }
                    cPtr[k] += sum;
                }
            }
        }
    }
    template<int n, int m, int p, int effC, int effA, int effB, typename T>
    void parallelMul10(T* c, const T* a, const T* b) {
        static constexpr mtl::size_t ib = 128;
        static constexpr mtl::size_t jb = 16;
        static constexpr mtl::size_t kb = 256;

        for (mtl::size_t i = 0; i < n; ++i) {
            for (mtl::size_t j = 0; j < p; ++j) {
                c[j + i * effC] = T{};
            }
        }

        auto newB = getNewB<ib, jb, kb>(b, n, m, p, effB);

        ThreadPool pool;
        for (mtl::size_t ii = 0; ii < n; ii += ib) {
            mtl::size_t iEnd = std::min(ii + ib, n);
            for (mtl::size_t kk = 0; kk < p; kk += kb) {
                mtl::size_t kEnd = std::min(kk + kb, p);
                pool.addTask([=]() {mul10Task<jb, n, m, p, effC, effA, effB>(c, a, b, newB, ii, iEnd, kk, kEnd); });
            }
        }
        pool.completeTasksAndStop();

        delete[] newB;
    }
}

namespace avx {
    template<typename T> void parallelMul(T* c, const T* a, const T* b, int n, int m, int p) {
        detail::parallelMul9(c, a, b, n, m, p, p, m, p);
    }
    template<typename T> void parallelMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        detail::parallelMul9(c, a, b, n, m, p, effC, effA, effB);
    }
    template<int n, int m, int p, int effC, int effA, int effB, typename T> void parallelMul(T* c, const T* a, const T* b) {
        detail::parallelMul9<n, m, p, effC, effA, effB>(c, a, b);
    }
    template<typename T> void mul(T* c, const T* a, const T* b, int n, int m, int p) {
        detail::mul9(c, a, b, n, m, p, p, m, p);
    }
    template<typename T> void mul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        detail::mul9(c, a, b, n, m, p, effC, effA, effB);
    }
    template<int n, int m, int p, int effC, int effA, int effB, typename T> void mul(T* c, const T* a, const T* b) {
        detail::mul9<n, m, p, effC, effA, effB>(c, a, b);
    }

    template<typename M1, typename M2>
    auto parallelMul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        if constexpr (avx::IsAvailable) {
            if constexpr (M1::HasConstexprSizes() && M2::HasConstexprSizes()) {
                auto result = a.template createNew<M1::CRow(), M2::CCol()>();
                parallelMul<M1::CRow(), M1::CCol(), M2::CCol(), result.CEffCol(), M1::CEffCol(), M2::CEffCol()>(result.data(), a.data(), b.data());
                return result;
            } else {
                auto result = a.createNew(a.rowCount(), b.columnCount());
                parallelMul(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), result.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
                return result;
            }
        } else {
            static_assert(always_false_v<M1>, avx_StaticAssertMessage);
        }
    }

    template<typename M1, typename M2>
    auto mul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        if constexpr (avx::IsAvailable) {
            if constexpr (M1::HasConstexprSizes() && M2::HasConstexprSizes()) {
                auto result = a.template createNew<M1::CRow(), M2::CCol()>();
                mul<M1::CRow(), M1::CCol(), M2::CCol(), result.CEffCol(), M1::CEffCol(), M2::CEffCol()>(result.data(), a.data(), b.data());
                return result;
            } else {
                auto result = a.createNew(a.rowCount(), b.columnCount());
                mul(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), result.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
                return result;
            }
        } else {
            static_assert(always_false_v<M1>, avx_StaticAssertMessage);
        }
    }
}

#endif

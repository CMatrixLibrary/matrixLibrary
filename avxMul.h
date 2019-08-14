#ifndef AVX_MUL_H
#define AVX_MUL_H
#include "avxSimd.h"

void avxMul2(int* result, const int* a, const int* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
    for (int i = 0; i < n; ++i) {
        int j = 0;
        for (; j <= q - 8; j += 8) {
            auto sum = AVX256::zero<int>();
            for (int k = 0; k < m; ++k) {
                auto aValue = AVX256::setAllElements(a[k + i * m]);
                auto bVectr = AVX256::loadUnaligned(&b[j + k * q]);
                sum += aValue * bVectr;
            }
            AVX256::storeUnaligned(&result[j + i * q], sum);
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
            auto sum1 = AVX256::zero<int>();
            auto sum2 = AVX256::zero<int>();
            for (int k = 0; k < m; ++k) {
                auto aValue = AVX256::setAllElements(a[k + i * m]);
                auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                sum1 += aValue * bVectr1;
                auto bVectr2 = AVX256::loadUnaligned(&b[8 + j + k * q]);
                sum2 += aValue * bVectr2;
            }
            AVX256::storeUnaligned(&result[j + i * q], sum1);
            AVX256::storeUnaligned(&result[8 + j + i * q], sum2);
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
            auto sum1 = AVX256::zero<int>();
            auto sum2 = AVX256::zero<int>();
            auto sum3 = AVX256::zero<int>();
            auto sum4 = AVX256::zero<int>();
            for (int k = 0; k < m; ++k) {
                auto aValue = AVX256::setAllElements(a[k + i * m]);
                auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                sum1 += aValue * bVectr1;
                auto bVectr2 = AVX256::loadUnaligned(&b[8 + j + k * q]);
                sum2 += aValue * bVectr2;
                auto bVectr3 = AVX256::loadUnaligned(&b[16 + j + k * q]);
                sum3 += aValue * bVectr3;
                auto bVectr4 = AVX256::loadUnaligned(&b[24 + j + k * q]);
                sum4 += aValue * bVectr4;
            }
            AVX256::storeUnaligned(&result[j + i * q], sum1);
            AVX256::storeUnaligned(&result[8 + j + i * q], sum2);
            AVX256::storeUnaligned(&result[16 + j + i * q], sum3);
            AVX256::storeUnaligned(&result[24 + j + i * q], sum4);
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
                int jEnd = lastjb-16;
                for (; j <= jEnd; j += 16) {
                    auto sum1 = AVX256::loadUnaligned(&result[j + i * q]);
                    auto sum2 = AVX256::loadUnaligned(&result[8 + j + i * q]);
                    for (int k = kk; k < lastKb; ++k) {
                        auto aValue = AVX256::setAllElements(a[k + i * m]);
                        auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                        sum1 += aValue * bVectr1;
                        auto bVectr2 = AVX256::loadUnaligned(&b[8 + j + k * q]);
                        sum2 += aValue * bVectr2;
                    }
                    AVX256::storeUnaligned(&result[j + i * q], sum1);
                    AVX256::storeUnaligned(&result[8 + j + i * q], sum2);
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
            for (; i <= n-2; i += 2) {
                int j = jj;
                int jEnd = lastjb-16;
                for (; j <= jEnd; j += 16) {
                    auto sum1i1 = AVX256::loadUnaligned(&result[j + i * q]);
                    auto sum1i2 = AVX256::loadUnaligned(&result[8 + j + i * q]);
                    auto sum2i1 = AVX256::loadUnaligned(&result[j + (i+1) * q]);
                    auto sum2i2 = AVX256::loadUnaligned(&result[8 + j + (i+1) * q]);
                    for (int k = kk; k < lastKb; ++k) {
                        auto aValue1 = AVX256::setAllElements(a[k + i * m]);
                        auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                        sum1i1 += aValue1 * bVectr1;
                        auto bVectr2 = AVX256::loadUnaligned(&b[8 + j + k * q]);
                        sum1i2 += aValue1 * bVectr2;
                        auto aValue2 = AVX256::setAllElements(a[k + (i + 1) * m]);
                        sum2i1 += aValue2 * bVectr1;
                        sum2i2 += aValue2 * bVectr2;
                    }
                    AVX256::storeUnaligned(&result[j + i * q], sum1i1);
                    AVX256::storeUnaligned(&result[8 + j + i * q], sum1i2);
                    AVX256::storeUnaligned(&result[j + (i+1) * q], sum1i1);
                    AVX256::storeUnaligned(&result[8 + j + (i+1) * q], sum1i2);
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
                    auto sum1 = AVX256::loadUnaligned(&result[j + i * q]);
                    auto sum2 = AVX256::loadUnaligned(&result[8 + j + i * q]);
                    for (int k = kk; k < lastKb; ++k) {
                        auto aValue = AVX256::setAllElements(a[k + i * m]);
                        auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                        sum1 += aValue * bVectr1;
                        auto bVectr2 = AVX256::loadUnaligned(&b[8 + j + k * q]);
                        sum2 += aValue * bVectr2;
                    }
                    AVX256::storeUnaligned(&result[j + i * q], sum1);
                    AVX256::storeUnaligned(&result[8 + j + i * q], sum2);
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
template<typename T> void avxMul7(T* result, const T* a, const T* b, int n, int m, int q) {
#ifdef AVX2_IS_AVAILABLE
    for (int i = 0; i < n*q; ++i) {
        result[i] = T{};
    }

    int ib = 256;
    int jb = 512;
    int kb = 16;

    for (int ii = 0; ii < n; ii += ib) {
        int lastib = std::min(ii + ib, n);
        for (int jj = 0; jj < q; jj += jb) {
            int lastjb = std::min(jj + jb, q);
            for (int kk = 0; kk < m; kk += kb) {
                int lastKb = std::min(kk + kb, m);
                int i = ii;
                for (; i <= lastib - 2; i += 2) {
                    int j = jj;
                    int jEnd = lastjb - 2*AVX256::packedCount<T>();
                    for (; j <= jEnd; j += 2 * AVX256::packedCount<T>()) {
                        auto sum1i1 = AVX256::loadUnaligned(&result[j + i * q]);
                        auto sum1i2 = AVX256::loadUnaligned(&result[AVX256::packedCount<T>() + j + i * q]);
                        auto sum2i1 = AVX256::loadUnaligned(&result[j + (i + 1) * q]);
                        auto sum2i2 = AVX256::loadUnaligned(&result[AVX256::packedCount<T>() + j + (i + 1) * q]);
                        for (int k = kk; k < lastKb; ++k) {
                            auto aValue1 = AVX256::setAllElements(a[k + i * m]);
                            auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                            sum1i1 = AVX256::fma(aValue1, bVectr1, sum1i1);
                            auto bVectr2 = AVX256::loadUnaligned(&b[AVX256::packedCount<T>() + j + k * q]);
                            sum1i2 = AVX256::fma(aValue1, bVectr2, sum1i2);
                            auto aValue2 = AVX256::setAllElements(a[k + (i + 1) * m]);
                            sum2i1 = AVX256::fma(aValue2, bVectr1, sum2i1);
                            sum2i2 = AVX256::fma(aValue2, bVectr2, sum2i2);
                        }
                        AVX256::storeUnaligned(&result[j + i * q], sum1i1);
                        AVX256::storeUnaligned(&result[AVX256::packedCount<T>() + j + i * q], sum1i2);
                        AVX256::storeUnaligned(&result[j + (i + 1) * q], sum1i1);
                        AVX256::storeUnaligned(&result[AVX256::packedCount<T>() + j + (i + 1) * q], sum1i2);
                    }
                    for (; j < lastjb; ++j) {
                        T sum{};
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
                if (i < lastib) {
                    int j = jj;
                    int jEnd = lastjb - 2 * AVX256::packedCount<T>();
                    for (; j <= jEnd; j += 2 * AVX256::packedCount<T>()) {
                        auto sum1 = AVX256::loadUnaligned(&result[j + i * q]);
                        auto sum2 = AVX256::loadUnaligned(&result[AVX256::packedCount<T>() + j + i * q]);
                        for (int k = kk; k < lastKb; ++k) {
                            auto aValue = AVX256::setAllElements(a[k + i * m]);
                            auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                            sum1 = AVX256::fma(aValue, bVectr1, sum1);
                            auto bVectr2 = AVX256::loadUnaligned(&b[AVX256::packedCount<T>() + j + k * q]);
                            sum2 = AVX256::fma(aValue, bVectr2, sum2);
                        }
                        AVX256::storeUnaligned(&result[j + i * q], sum1);
                        AVX256::storeUnaligned(&result[AVX256::packedCount<T>() + j + i * q], sum2);
                    }
                    for (; j < lastjb; ++j) {
                        T sum{};
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
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

    T* mat2 = AVX256::alignedArrayAlloc<T>(n*n);
    int m2idx = 0;
    for (int jj = 0; jj < n; jj += jb) {
        for (int kk = 0; kk < n; kk += kb) {
            for (int j = jj; j < jj + jb; j += 2 * AVX256::packedCount<T>()) {
                for (int k = kk; k < kk + kb; k++) {
                    auto vecA_mat2 = AVX256::loadUnaligned(&b[j + k * n]);
                    auto vecB_mat2 = AVX256::loadUnaligned(&b[AVX256::packedCount<T>() + j + k * n]);
                    AVX256::storeUnaligned(&mat2[m2idx], vecA_mat2);
                    AVX256::storeUnaligned(&mat2[m2idx + AVX256::packedCount<T>()], vecB_mat2);
                    m2idx += 2 * AVX256::packedCount<T>();
                }
            }
        }
    }

    for (int ii = 0; ii < n; ii += ib) {
        for (int jj = 0; jj < n; jj += jb) {
            for (int kk = 0; kk < n; kk += kb) {
                for (int i = ii; i < ii + ib; i += 2) {
                    for (int j = jj; j < jj + jb; j += 2 * AVX256::packedCount<T>()) {
                        int m2idx = (j - jj) * kb + kk * jb + jj * n;
                        auto sumA_1 = AVX256::loadUnaligned(&result[i*n + j]);
                        auto sumB_1 = AVX256::loadUnaligned(&result[i*n + j + AVX256::packedCount<T>()]);
                        auto sumA_2 = AVX256::loadUnaligned(&result[(i + 1)*n + j]);
                        auto sumB_2 = AVX256::loadUnaligned(&result[(i + 1)*n + j + AVX256::packedCount<T>()]);
                        for (int k = kk; k < kk + kb; k++) {
                            auto bc_mat1_1 = AVX256::setAllElements(a[i*n + k]);
                            auto vecA_mat2 = AVX256::loadUnaligned(&mat2[m2idx]);
                            auto vecB_mat2 = AVX256::loadUnaligned(&mat2[m2idx + AVX256::packedCount<T>()]);
                            sumA_1 += bc_mat1_1 * vecA_mat2;
                            sumB_1 += bc_mat1_1 * vecB_mat2;
                            auto bc_mat1_2 = AVX256::setAllElements(a[(i+1)*n + k]);
                            sumA_2 += bc_mat1_2 * vecA_mat2;
                            sumB_2 += bc_mat1_2 * vecB_mat2;
                            m2idx += 2 * AVX256::packedCount<T>();
                        }
                        AVX256::storeUnaligned(&result[i*n + j], sumA_1);
                        AVX256::storeUnaligned(&result[i*n + j + AVX256::packedCount<T>()], sumB_1);
                        AVX256::storeUnaligned(&result[(i + 1)*n + j], sumA_2);
                        AVX256::storeUnaligned(&result[(i + 1) * n + j + AVX256::packedCount<T>()], sumB_2);
                    }
                }
            }
        }
    }
    AVX256::alignedArrayDealloc(mat2);
#endif
}
template<int n, int m, int q, typename T> void avxMul7(T* result, const T* a, const T* b) {
#ifdef AVX2_IS_AVAILABLE
    for (int i = 0; i < n*q; ++i) {
        result[i] = T{};
    }

    int ib = 256;
    int jb = 512;
    int kb = 16;

    for (int ii = 0; ii < n; ii += ib) {
        int lastib = std::min(ii + ib, n);
        for (int jj = 0; jj < q; jj += jb) {
            int lastjb = std::min(jj + jb, q);
            for (int kk = 0; kk < m; kk += kb) {
                int lastKb = std::min(kk + kb, m);
                int i = ii;
                for (; i <= lastib - 2; i += 2) {
                    int j = jj;
                    int jEnd = lastjb - 2*AVX256::packedCount<T>();
                    for (; j <= jEnd; j += 2 * AVX256::packedCount<T>()) {
                        auto sum1i1 = AVX256::loadUnaligned(&result[j + i * q]);
                        auto sum1i2 = AVX256::loadUnaligned(&result[AVX256::packedCount<T>() + j + i * q]);
                        auto sum2i1 = AVX256::loadUnaligned(&result[j + (i + 1) * q]);
                        auto sum2i2 = AVX256::loadUnaligned(&result[AVX256::packedCount<T>() + j + (i + 1) * q]);
                        for (int k = kk; k < lastKb; ++k) {
                            auto aValue1 = AVX256::setAllElements(a[k + i * m]);
                            auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                            sum1i1 = AVX256::fma(aValue1, bVectr1, sum1i1);
                            auto bVectr2 = AVX256::loadUnaligned(&b[AVX256::packedCount<T>() + j + k * q]);
                            sum1i2 = AVX256::fma(aValue1, bVectr2, sum1i2);
                            auto aValue2 = AVX256::setAllElements(a[k + (i + 1) * m]);
                            sum2i1 = AVX256::fma(aValue2, bVectr1, sum2i1);
                            sum2i2 = AVX256::fma(aValue2, bVectr2, sum2i2);
                        }
                        AVX256::storeUnaligned(&result[j + i * q], sum1i1);
                        AVX256::storeUnaligned(&result[AVX256::packedCount<T>() + j + i * q], sum1i2);
                        AVX256::storeUnaligned(&result[j + (i + 1) * q], sum1i1);
                        AVX256::storeUnaligned(&result[AVX256::packedCount<T>() + j + (i + 1) * q], sum1i2);
                    }
                    for (; j < lastjb; ++j) {
                        T sum{};
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
                if (i < lastib) {
                    int j = jj;
                    int jEnd = lastjb - 2 * AVX256::packedCount<T>();
                    for (; j <= jEnd; j += 2 * AVX256::packedCount<T>()) {
                        auto sum1 = AVX256::loadUnaligned(&result[j + i * q]);
                        auto sum2 = AVX256::loadUnaligned(&result[AVX256::packedCount<T>() + j + i * q]);
                        for (int k = kk; k < lastKb; ++k) {
                            auto aValue = AVX256::setAllElements(a[k + i * m]);
                            auto bVectr1 = AVX256::loadUnaligned(&b[j + k * q]);
                            sum1 = AVX256::fma(aValue, bVectr1, sum1);
                            auto bVectr2 = AVX256::loadUnaligned(&b[AVX256::packedCount<T>() + j + k * q]);
                            sum2 = AVX256::fma(aValue, bVectr2, sum2);
                        }
                        AVX256::storeUnaligned(&result[j + i * q], sum1);
                        AVX256::storeUnaligned(&result[AVX256::packedCount<T>() + j + i * q], sum2);
                    }
                    for (; j < lastjb; ++j) {
                        T sum{};
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
            }
        }
    }
#endif
}

#endif

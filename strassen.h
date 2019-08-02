#pragma once
#include "MatrixInterface.h"
#include "matrixOperators.h"
#include "MatrixExtendedFunctions.h"
#include "StackAllocator.h"
#include "avxSimd.h"

void avxMul2(int* result, const int* a, const int* b, int n, int m, int q) {
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
}
void avxMul3(int* result, const int* a, const int* b, int n, int m, int q) {
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
}
void avxMul4(int* result, const int* a, const int* b, int n, int m, int q) {
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
}
void avxMul5(int* result, const int* a, const int* b, int n, int m, int q) {
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
}

void avxMul6(int* result, const int* a, const int* b, int n, int m, int q) {
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
}
void avxMul7(int* result, const int* a, const int* b, int n, int m, int q) {
    for (int i = 0; i < n*q; ++i) {
        result[i] = 0;
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
                    int jEnd = lastjb - 16;
                    for (; j <= jEnd; j += 16) {
                        auto sum1i1 = AVX256::loadUnaligned(&result[j + i * q]);
                        auto sum1i2 = AVX256::loadUnaligned(&result[8 + j + i * q]);
                        auto sum2i1 = AVX256::loadUnaligned(&result[j + (i + 1) * q]);
                        auto sum2i2 = AVX256::loadUnaligned(&result[8 + j + (i + 1) * q]);
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
                        AVX256::storeUnaligned(&result[j + (i + 1) * q], sum1i1);
                        AVX256::storeUnaligned(&result[8 + j + (i + 1) * q], sum1i2);
                    }
                    for (; j < lastjb; ++j) {
                        int sum = 0;
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
                if (i < lastib) {
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
    }
}
void avxMul8(int* result, const int* a, const int* b, int n) {
    for (int i = 0; i < n*n; ++i) {
        result[i] = 0;
    }

    int ib = std::min(256, n);
    int jb = std::min(512, n);
    int kb = std::min(16, n);

    int* mat2 = AVX256::AlignedArrayAlloc<int>(n*n);
    int m2idx = 0;
    for (int jj = 0; jj < n; jj += jb) {
        for (int kk = 0; kk < n; kk += kb) {
            for (int j = jj; j < jj + jb; j += 16) {
                for (int k = kk; k < kk + kb; k++) {
                    auto vecA_mat2 = AVX256::loadUnaligned(&b[j + k * n]);
                    auto vecB_mat2 = AVX256::loadUnaligned(&b[8 + j + k * n]);
                    AVX256::storeUnaligned(&mat2[m2idx], vecA_mat2);
                    AVX256::storeUnaligned(&mat2[m2idx + 8], vecB_mat2);
                    m2idx += 16;
                }
            }
        }
    }

    for (int ii = 0; ii < n; ii += ib) {
        for (int jj = 0; jj < n; jj += jb) {
            for (int kk = 0; kk < n; kk += kb) {
                for (int i = ii; i < ii + ib; i += 2) {
                    for (int j = jj; j < jj + jb; j += 16) {
                        int m2idx = (j - jj) * kb + kk * jb + jj * n;
                        auto sumA_1 = AVX256::loadUnaligned(&result[i*n + j]);
                        auto sumB_1 = AVX256::loadUnaligned(&result[i*n + j + 8]);
                        auto sumA_2 = AVX256::loadUnaligned(&result[(i + 1)*n + j]);
                        auto sumB_2 = AVX256::loadUnaligned(&result[(i + 1)*n + j + 8]);
                        for (int k = kk; k < kk + kb; k++) {
                            auto bc_mat1_1 = AVX256::setAllElements(a[i*n + k]);
                            auto vecA_mat2 = AVX256::loadUnaligned(&mat2[m2idx]);
                            auto vecB_mat2 = AVX256::loadUnaligned(&mat2[m2idx + 8]);
                            sumA_1 += bc_mat1_1 * vecA_mat2;
                            sumB_1 += bc_mat1_1 * vecB_mat2;
                            auto bc_mat1_2 = AVX256::setAllElements(a[(i+1)*n + k]);
                            sumA_2 += bc_mat1_2 * vecA_mat2;
                            sumB_2 += bc_mat1_2 * vecB_mat2;
                            m2idx += 16;
                        }
                        AVX256::storeUnaligned(&result[i*n + j], sumA_1);
                        AVX256::storeUnaligned(&result[i*n + j + 8], sumB_1);
                        AVX256::storeUnaligned(&result[(i + 1)*n + j], sumA_2);
                        AVX256::storeUnaligned(&result[(i + 1) * n + j + 8], sumB_2);
                    }
                }
            }
        }
    }
    AVX256::AlignedArrayDealloc(mat2);
}
template<int n, int m, int q> void avxMul7(int* result, const int* a, const int* b) {
    for (int i = 0; i < n*q; ++i) {
        result[i] = 0;
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
                    int jEnd = lastjb - 16;
                    for (; j <= jEnd; j += 16) {
                        auto sum1i1 = AVX256::loadUnaligned(&result[j + i * q]);
                        auto sum1i2 = AVX256::loadUnaligned(&result[8 + j + i * q]);
                        auto sum2i1 = AVX256::loadUnaligned(&result[j + (i + 1) * q]);
                        auto sum2i2 = AVX256::loadUnaligned(&result[8 + j + (i + 1) * q]);
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
                        AVX256::storeUnaligned(&result[j + (i + 1) * q], sum1i1);
                        AVX256::storeUnaligned(&result[8 + j + (i + 1) * q], sum1i2);
                    }
                    for (; j < lastjb; ++j) {
                        int sum = 0;
                        for (int k = 0; k < m; ++k) {
                            sum += a[k + i * m] * b[j + k * q];
                        }
                        result[j + i * q] = sum;
                    }
                }
                if (i < lastib) {
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
    }
}

template<typename M1, typename M2> 
auto strassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if (steps <= 0) {
        return a * b;
    }

    auto dA = matrixDivide<2, 2>(a);
    auto dB = matrixDivide<2, 2>(b);

    auto m1 = strassen(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1], steps - 1);
    auto m2 = strassen(dA[1][0] + dA[1][1], dB[0][0], steps - 1);
    auto m3 = strassen(dA[0][0], dB[0][1] - dB[1][1], steps - 1);
    auto m4 = strassen(dA[1][1], dB[1][0] - dB[0][0], steps - 1);
    auto m5 = strassen(dA[0][0] + dA[0][1], dB[1][1], steps - 1);
    auto m6 = strassen(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1], steps - 1);
    auto m7 = strassen(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1], steps - 1);

    Matrix<typename M1::ValueType> c(a.rowCount(), b.columnCount());
    auto dC = matrixDivideView<2, 2>(c);
    dC[0][0].copy(m1 + m4 - m5 + m7);
    dC[0][1].copy(m3 + m5);
    dC[1][0].copy(m2 + m4);
    dC[1][1].copy(m1 + m3 - m2 + m6);

    return c;
}

template<int Steps, typename M1, typename M2> auto strassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, std::false_type) {
    return a * b;
}
template<int Steps, typename M1, typename M2> auto strassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, std::true_type unused = std::true_type()) {
    if (Steps <= 0) {
        return a * b;
    }

    auto dA = matrixDivide<2, 2>(a);
    auto dB = matrixDivide<2, 2>(b);

    auto m1 = strassen<Steps - 1>(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m2 = strassen<Steps - 1>(dA[1][0] + dA[1][1], dB[0][0], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m3 = strassen<Steps - 1>(dA[0][0], dB[0][1] - dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m4 = strassen<Steps - 1>(dA[1][1], dB[1][0] - dB[0][0], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m5 = strassen<Steps - 1>(dA[0][0] + dA[0][1], dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m6 = strassen<Steps - 1>(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m7 = strassen<Steps - 1>(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());

    Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> c;
    auto dC = matrixDivideView<2, 2>(c);
    dC[0][0].copy(m1 + m4 - m5 + m7);
    dC[0][1].copy(m3 + m5);
    dC[1][0].copy(m2 + m4);
    dC[1][1].copy(m1 + m3 - m2 + m6);

    return c;
}

template<typename M1, typename M2> 
auto strassenAvx(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if (steps <= 0) {
        Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
        avxMul7(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount());
        return result;
    }

    auto dA = matrixDivide<2, 2>(a);
    auto dB = matrixDivide<2, 2>(b);

    auto m1 = strassenAvx(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1], steps - 1);
    auto m2 = strassenAvx(dA[1][0] + dA[1][1], dB[0][0], steps - 1);
    auto m3 = strassenAvx(dA[0][0], dB[0][1] - dB[1][1], steps - 1);
    auto m4 = strassenAvx(dA[1][1], dB[1][0] - dB[0][0], steps - 1);
    auto m5 = strassenAvx(dA[0][0] + dA[0][1], dB[1][1], steps - 1);
    auto m6 = strassenAvx(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1], steps - 1);
    auto m7 = strassenAvx(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1], steps - 1);

    Matrix<typename M1::ValueType> c(a.rowCount(), b.columnCount());
    auto dC = matrixDivideView<2, 2>(c);
    dC[0][0].copy(m1 + m4 - m5 + m7);
    dC[0][1].copy(m3 + m5);
    dC[1][0].copy(m2 + m4);
    dC[1][1].copy(m1 + m3 - m2 + m6);

    return c;
}

template<int Steps, typename M1, typename M2> auto strassenAvx(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, std::false_type) {
    Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
    avxMul7<M1::CRow(), M1::CCol(), M2::CCol()>(result.data(), a.data(), b.data());
    return result;
}
template<int Steps, typename M1, typename M2> auto strassenAvx(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, std::true_type unused = std::true_type()) {
    if (Steps <= 0) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        avxMul7<M1::CRow(), M1::CCol(), M2::CCol()>(result.data(), a.data(), b.data());
        return result;
    }

    auto dA = matrixDivide<2, 2>(a);
    auto dB = matrixDivide<2, 2>(b);

    auto m1 = strassenAvx<Steps - 1>(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m2 = strassenAvx<Steps - 1>(dA[1][0] + dA[1][1], dB[0][0], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m3 = strassenAvx<Steps - 1>(dA[0][0], dB[0][1] - dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m4 = strassenAvx<Steps - 1>(dA[1][1], dB[1][0] - dB[0][0], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m5 = strassenAvx<Steps - 1>(dA[0][0] + dA[0][1], dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m6 = strassenAvx<Steps - 1>(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m7 = strassenAvx<Steps - 1>(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());

    Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> c;
    auto dC = matrixDivideView<2, 2>(c);
    dC[0][0].copy(m1 + m4 - m5 + m7);
    dC[0][1].copy(m3 + m5);
    dC[1][0].copy(m2 + m4);
    dC[1][1].copy(m1 + m3 - m2 + m6);

    return c;
}


void copy(int* dst, int* src, int n, int m, int effDst, int effSrc) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            dst[j + i * effDst] = src[j + i * effSrc];
        }
    }
}
template<int Rows, int Columns>
std::array<std::array<int*, Columns>, Rows> lowLevelStrassenDivide(int* a, int n, int m, StackAllocator& allocator) {
    std::array<std::array<int*, Columns>, Rows> result;
    auto rowSize = n / Rows;
    auto columnSize = m / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[y][x] = allocator.alloc(rowSize*columnSize);
            copy(result[y][x], a + x * columnSize + y * rowSize*m, rowSize, columnSize, columnSize, m);
        }
    }
    return result;
}
template<int Rows, int Columns>
std::array<std::array<int*, Columns>, Rows> lowLevelStrassenDivideView(int* a, int n, int m) {
    std::array<std::array<int*, Columns>, Rows> result;
    auto rowSize = n / Rows;
    auto columnSize = m / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[y][x] = a + x * columnSize + y * rowSize*m;
        }
    }
    return result;
}
int* lowLevelStrassenAdd(const int* a, const int* b, int n, int m, StackAllocator& allocator) {
    int* result = allocator.alloc(n*m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[j + i * m] = a[j + i * m] + b[j + i * m];
        }
    }
    return result;
}
int* lowLevelStrassenSub(const int* a, const int* b, int n, int m, StackAllocator& allocator) {
    int* result = allocator.alloc(n*m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[j + i * m] = a[j + i * m] - b[j + i * m];
        }
    }
    return result;
}
template<int n, int m> int* lowLevelStrassenAdd(const int* a, const int* b, StackAllocator& allocator) {
    int* result = allocator.alloc(n*m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[j + i * m] = a[j + i * m] + b[j + i * m];
        }
    }
    return result;
}
template<int n, int m> int* lowLevelStrassenSub(const int* a, const int* b, StackAllocator& allocator) {
    int* result = allocator.alloc(n*m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[j + i * m] = a[j + i * m] - b[j + i * m];
        }
    }
    return result;
}

void bestStrassenMul(int* result, const int* a, const int* b, int n, int m, int q) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < q; ++j) {
            int sum = 0;
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * m] * b[j + k * q];
            }
            result[j + i * q] = sum;
        }
    }
}
template<int n, int m, int q> void bestStrassenMul(int* result, const int* a, const int* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < q; ++j) {
            int sum = 0;
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * m] * b[j + k * q];
            }
            result[j + i * q] = sum;
        }
    }
}
void lowLevelStrassen(int* a, int* b, int n, int m, int q, int steps, int* c, StackAllocator& allocator) {
    if (steps <= 0) {
        bestStrassenMul(c, a, b, n, m, q);
        return;
    }

    auto dA = lowLevelStrassenDivide<2, 2>(a, n, m, allocator);
    auto dB = lowLevelStrassenDivide<2, 2>(b, m, q, allocator);

    int halfN = n / 2;
    int halfM = m / 2;
    int halfQ = q / 2;

    auto m1 = allocator.alloc(halfN * halfQ);
    auto m1_a = lowLevelStrassenAdd(dA[0][0], dA[1][1], halfN, halfM, allocator);
    auto m1_b = lowLevelStrassenAdd(dB[0][0], dB[1][1], halfM, halfQ, allocator);
    lowLevelStrassen(m1_a, m1_b, halfN, halfM, halfQ, steps - 1, m1, allocator);
    allocator.dealloc(m1_a);

    auto m2 = allocator.alloc(halfN * halfQ);
    auto m2_a = lowLevelStrassenAdd(dA[1][0], dA[1][1], halfN, halfM, allocator);
    lowLevelStrassen(m2_a, dB[0][0], halfN, halfM, halfQ, steps - 1, m2, allocator);
    allocator.dealloc(m2_a);

    auto m3 = allocator.alloc(halfN * halfQ);
    auto m3_b = lowLevelStrassenSub(dB[0][1], dB[1][1], halfM, halfQ, allocator);
    lowLevelStrassen(dA[0][0], m3_b, halfN, halfM, halfQ, steps - 1, m3, allocator);
    allocator.dealloc(m3_b);

    auto m4 = allocator.alloc(halfN * halfQ);
    auto m4_b = lowLevelStrassenSub(dB[1][0], dB[0][0], halfM, halfQ, allocator);
    lowLevelStrassen(dA[1][1], m4_b, halfN, halfM, halfQ, steps - 1, m4, allocator);
    allocator.dealloc(m4_b);

    auto m5 = allocator.alloc(halfN * halfQ);
    auto m5_a = lowLevelStrassenAdd(dA[0][0], dA[0][1], halfN, halfM, allocator);
    lowLevelStrassen(m5_a, dB[1][1], halfN, halfM, halfQ, steps - 1, m5, allocator);
    allocator.dealloc(m5_a);

    auto m6 = allocator.alloc(halfN * halfQ);
    auto m6_a = lowLevelStrassenSub(dA[1][0], dA[0][0], halfN, halfM, allocator);
    auto m6_b = lowLevelStrassenAdd(dB[0][0], dB[0][1], halfM, halfQ, allocator);
    lowLevelStrassen(m6_a, m6_b, halfN, halfM, halfQ, steps - 1, m6, allocator);
    allocator.dealloc(m6_a);

    auto m7 = allocator.alloc(halfN * halfQ);
    auto m7_a = lowLevelStrassenSub(dA[0][1], dA[1][1], halfN, halfM, allocator);
    auto m7_b = lowLevelStrassenAdd(dB[1][0], dB[1][1], halfM, halfQ, allocator);
    lowLevelStrassen(m7_a, m7_b, halfN, halfM, halfQ, steps - 1, m7, allocator);
    allocator.dealloc(m7_a);

    auto dC = lowLevelStrassenDivideView<2, 2>(c, n, q);
    for (int i = 0; i < halfN; ++i) {
        for (int j = 0; j < halfQ; ++j) {
            int dstIndex = j + i * q;
            int index = j + i * halfQ;
            dC[0][0][dstIndex] = m1[index] + m4[index] - m5[index] + m7[index];
            dC[0][1][dstIndex] = m3[index] + m5[index];
            dC[1][0][dstIndex] = m2[index] + m4[index];
            dC[1][1][dstIndex] = m1[index] + m3[index] - m2[index] + m6[index];
        }
    }

    allocator.dealloc(dA[0][0]);
}
void lowLevelStrassen(int* result, int* a, int* b, int n, int m, int q, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += 9*StackAllocator::Allign(eN*eM) + 9* StackAllocator::Allign(eM*eQ) + 7* StackAllocator::Allign(eN*eQ);
    }
    StackAllocator allocator(expected);
    lowLevelStrassen(a, b, n, m, q, steps, result, allocator);
}

template<int n, int m, int q> void lowLevelStrassen(int* a, int* b, int steps, int* c, StackAllocator& allocator) {
    if (steps <= 0) {
        bestStrassenMul<n, m, q>(c, a, b);
        return;
    }

    auto dA = lowLevelStrassenDivide<2, 2>(a, n, m, allocator);
    auto dB = lowLevelStrassenDivide<2, 2>(b, m, q, allocator);

    constexpr int halfN = n / 2;
    constexpr int halfM = m / 2;
    constexpr int halfQ = q / 2;

    auto m1 = allocator.alloc(halfN * halfQ);
    auto m1_a = lowLevelStrassenAdd<halfN, halfM>(dA[0][0], dA[1][1], allocator);
    auto m1_b = lowLevelStrassenAdd<halfM, halfQ>(dB[0][0], dB[1][1], allocator);
    lowLevelStrassen<halfN, halfM, halfQ>(m1_a, m1_b, steps - 1, m1, allocator);
    allocator.dealloc(m1_a);

    auto m2 = allocator.alloc(halfN * halfQ);
    auto m2_a = lowLevelStrassenAdd<halfN, halfM>(dA[1][0], dA[1][1], allocator);
    lowLevelStrassen<halfN, halfM, halfQ>(m2_a, dB[0][0], steps - 1, m2, allocator);
    allocator.dealloc(m2_a);

    auto m3 = allocator.alloc(halfN * halfQ);
    auto m3_b = lowLevelStrassenSub<halfM, halfQ>(dB[0][1], dB[1][1], allocator);
    lowLevelStrassen<halfN, halfM, halfQ>(dA[0][0], m3_b, steps - 1, m3, allocator);
    allocator.dealloc(m3_b);

    auto m4 = allocator.alloc(halfN * halfQ);
    auto m4_b = lowLevelStrassenSub<halfM, halfQ>(dB[1][0], dB[0][0], allocator);
    lowLevelStrassen<halfN, halfM, halfQ>(dA[1][1], m4_b, steps - 1, m4, allocator);
    allocator.dealloc(m4_b);

    auto m5 = allocator.alloc(halfN * halfQ);
    auto m5_a = lowLevelStrassenAdd<halfN, halfM>(dA[0][0], dA[0][1], allocator);
    lowLevelStrassen<halfN, halfM, halfQ>(m5_a, dB[1][1], steps - 1, m5, allocator);
    allocator.dealloc(m5_a);

    auto m6 = allocator.alloc(halfN * halfQ);
    auto m6_a = lowLevelStrassenSub<halfN, halfM>(dA[1][0], dA[0][0], allocator);
    auto m6_b = lowLevelStrassenAdd<halfM, halfQ>(dB[0][0], dB[0][1], allocator);
    lowLevelStrassen<halfN, halfM, halfQ>(m6_a, m6_b, steps - 1, m6, allocator);
    allocator.dealloc(m6_a);

    auto m7 = allocator.alloc(halfN * halfQ);
    auto m7_a = lowLevelStrassenSub<halfN, halfM>(dA[0][1], dA[1][1], allocator);
    auto m7_b = lowLevelStrassenAdd<halfM, halfQ>(dB[1][0], dB[1][1], allocator);
    lowLevelStrassen<halfN, halfM, halfQ>(m7_a, m7_b, steps - 1, m7, allocator);
    allocator.dealloc(m7_a);

    auto dC = lowLevelStrassenDivideView<2, 2>(c, n, q);
    for (int i = 0; i < halfN; ++i) {
        for (int j = 0; j < halfQ; ++j) {
            int dstIndex = j + i * q;
            int index = j + i * halfQ;
            dC[0][0][dstIndex] = m1[index] + m4[index] - m5[index] + m7[index];
            dC[0][1][dstIndex] = m3[index] + m5[index];
            dC[1][0][dstIndex] = m2[index] + m4[index];
            dC[1][1][dstIndex] = m1[index] + m3[index] - m2[index] + m6[index];
        }
    }

    allocator.dealloc(dA[0][0]);
}
template<int n, int m, int q> void lowLevelStrassen(int* result, int* a, int* b, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += 9*StackAllocator::Allign(eN*eM) + 9* StackAllocator::Allign(eM*eQ) + 7* StackAllocator::Allign(eN*eQ);
    }
    StackAllocator allocator(expected);
    lowLevelStrassen<n, m, q>(a, b, steps, result, allocator);
}

template<typename M1, typename M2>
auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        lowLevelStrassen<M1::CRow(), M1::CCol(), M2::CCol()>(result.data(), a.data(), b.data(), steps);
        return result;
    } else {
        auto result = a.createNew(a.rowCount(), b.columnCount());
        lowLevelStrassen(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), steps);
        return result;
    }
}

void lowLevelAvxStrassen(int* a, int* b, int n, int m, int q, int steps, int* c, StackAllocator& allocator) {
    if (steps <= 0) {
        avxMul7(c, a, b, n, m, q);
        return;
    }

    auto dA = lowLevelStrassenDivide<2, 2>(a, n, m, allocator);
    auto dB = lowLevelStrassenDivide<2, 2>(b, m, q, allocator);

    int halfN = n / 2;
    int halfM = m / 2;
    int halfQ = q / 2;

    auto m1 = allocator.alloc(halfN * halfQ);
    auto m1_a = lowLevelStrassenAdd(dA[0][0], dA[1][1], halfN, halfM, allocator);
    auto m1_b = lowLevelStrassenAdd(dB[0][0], dB[1][1], halfM, halfQ, allocator);
    lowLevelAvxStrassen(m1_a, m1_b, halfN, halfM, halfQ, steps - 1, m1, allocator);
    allocator.dealloc(m1_a);

    auto m2 = allocator.alloc(halfN * halfQ);
    auto m2_a = lowLevelStrassenAdd(dA[1][0], dA[1][1], halfN, halfM, allocator);
    lowLevelAvxStrassen(m2_a, dB[0][0], halfN, halfM, halfQ, steps - 1, m2, allocator);
    allocator.dealloc(m2_a);

    auto m3 = allocator.alloc(halfN * halfQ);
    auto m3_b = lowLevelStrassenSub(dB[0][1], dB[1][1], halfM, halfQ, allocator);
    lowLevelAvxStrassen(dA[0][0], m3_b, halfN, halfM, halfQ, steps - 1, m3, allocator);
    allocator.dealloc(m3_b);

    auto m4 = allocator.alloc(halfN * halfQ);
    auto m4_b = lowLevelStrassenSub(dB[1][0], dB[0][0], halfM, halfQ, allocator);
    lowLevelAvxStrassen(dA[1][1], m4_b, halfN, halfM, halfQ, steps - 1, m4, allocator);
    allocator.dealloc(m4_b);

    auto m5 = allocator.alloc(halfN * halfQ);
    auto m5_a = lowLevelStrassenAdd(dA[0][0], dA[0][1], halfN, halfM, allocator);
    lowLevelAvxStrassen(m5_a, dB[1][1], halfN, halfM, halfQ, steps - 1, m5, allocator);
    allocator.dealloc(m5_a);

    auto m6 = allocator.alloc(halfN * halfQ);
    auto m6_a = lowLevelStrassenSub(dA[1][0], dA[0][0], halfN, halfM, allocator);
    auto m6_b = lowLevelStrassenAdd(dB[0][0], dB[0][1], halfM, halfQ, allocator);
    lowLevelAvxStrassen(m6_a, m6_b, halfN, halfM, halfQ, steps - 1, m6, allocator);
    allocator.dealloc(m6_a);

    auto m7 = allocator.alloc(halfN * halfQ);
    auto m7_a = lowLevelStrassenSub(dA[0][1], dA[1][1], halfN, halfM, allocator);
    auto m7_b = lowLevelStrassenAdd(dB[1][0], dB[1][1], halfM, halfQ, allocator);
    lowLevelAvxStrassen(m7_a, m7_b, halfN, halfM, halfQ, steps - 1, m7, allocator);
    allocator.dealloc(m7_a);

    auto dC = lowLevelStrassenDivideView<2, 2>(c, n, q);
    for (int i = 0; i < halfN; ++i) {
        for (int j = 0; j < halfQ; ++j) {
            int dstIndex = j + i * q;
            int index = j + i * halfQ;
            dC[0][0][dstIndex] = m1[index] + m4[index] - m5[index] + m7[index];
            dC[0][1][dstIndex] = m3[index] + m5[index];
            dC[1][0][dstIndex] = m2[index] + m4[index];
            dC[1][1][dstIndex] = m1[index] + m3[index] - m2[index] + m6[index];
        }
    }

    allocator.dealloc(dA[0][0]);
}
void lowLevelAvxStrassen(int* result, int* a, int* b, int n, int m, int q, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += 9*StackAllocator::Allign(eN*eM) + 9* StackAllocator::Allign(eM*eQ) + 7* StackAllocator::Allign(eN*eQ);
    }
    StackAllocator allocator(expected);
    lowLevelAvxStrassen(a, b, n, m, q, steps, result, allocator);
}

template<int n, int m, int q> void lowLevelAvxStrassen(int* a, int* b, int steps, int* c, StackAllocator& allocator) {
    if (steps <= 0) {
        avxMul7<n, m, q>(c, a, b);
        return;
    }

    auto dA = lowLevelStrassenDivide<2, 2>(a, n, m, allocator);
    auto dB = lowLevelStrassenDivide<2, 2>(b, m, q, allocator);

    constexpr int halfN = n / 2;
    constexpr int halfM = m / 2;
    constexpr int halfQ = q / 2;

    auto m1 = allocator.alloc(halfN * halfQ);
    auto m1_a = lowLevelStrassenAdd<halfN, halfM>(dA[0][0], dA[1][1], allocator);
    auto m1_b = lowLevelStrassenAdd<halfM, halfQ>(dB[0][0], dB[1][1], allocator);
    lowLevelAvxStrassen<halfN, halfM, halfQ>(m1_a, m1_b, steps - 1, m1, allocator);
    allocator.dealloc(m1_a);

    auto m2 = allocator.alloc(halfN * halfQ);
    auto m2_a = lowLevelStrassenAdd<halfN, halfM>(dA[1][0], dA[1][1], allocator);
    lowLevelAvxStrassen<halfN, halfM, halfQ>(m2_a, dB[0][0], steps - 1, m2, allocator);
    allocator.dealloc(m2_a);

    auto m3 = allocator.alloc(halfN * halfQ);
    auto m3_b = lowLevelStrassenSub<halfM, halfQ>(dB[0][1], dB[1][1], allocator);
    lowLevelAvxStrassen<halfN, halfM, halfQ>(dA[0][0], m3_b, steps - 1, m3, allocator);
    allocator.dealloc(m3_b);

    auto m4 = allocator.alloc(halfN * halfQ);
    auto m4_b = lowLevelStrassenSub<halfM, halfQ>(dB[1][0], dB[0][0], allocator);
    lowLevelAvxStrassen<halfN, halfM, halfQ>(dA[1][1], m4_b, steps - 1, m4, allocator);
    allocator.dealloc(m4_b);

    auto m5 = allocator.alloc(halfN * halfQ);
    auto m5_a = lowLevelStrassenAdd<halfN, halfM>(dA[0][0], dA[0][1], allocator);
    lowLevelAvxStrassen<halfN, halfM, halfQ>(m5_a, dB[1][1], steps - 1, m5, allocator);
    allocator.dealloc(m5_a);

    auto m6 = allocator.alloc(halfN * halfQ);
    auto m6_a = lowLevelStrassenSub<halfN, halfM>(dA[1][0], dA[0][0], allocator);
    auto m6_b = lowLevelStrassenAdd<halfM, halfQ>(dB[0][0], dB[0][1], allocator);
    lowLevelAvxStrassen<halfN, halfM, halfQ>(m6_a, m6_b, steps - 1, m6, allocator);
    allocator.dealloc(m6_a);

    auto m7 = allocator.alloc(halfN * halfQ);
    auto m7_a = lowLevelStrassenSub<halfN, halfM>(dA[0][1], dA[1][1], allocator);
    auto m7_b = lowLevelStrassenAdd<halfM, halfQ>(dB[1][0], dB[1][1], allocator);
    lowLevelAvxStrassen<halfN, halfM, halfQ>(m7_a, m7_b, steps - 1, m7, allocator);
    allocator.dealloc(m7_a);

    auto dC = lowLevelStrassenDivideView<2, 2>(c, n, q);
    for (int i = 0; i < halfN; ++i) {
        for (int j = 0; j < halfQ; ++j) {
            int dstIndex = j + i * q;
            int index = j + i * halfQ;
            dC[0][0][dstIndex] = m1[index] + m4[index] - m5[index] + m7[index];
            dC[0][1][dstIndex] = m3[index] + m5[index];
            dC[1][0][dstIndex] = m2[index] + m4[index];
            dC[1][1][dstIndex] = m1[index] + m3[index] - m2[index] + m6[index];
        }
    }

    allocator.dealloc(dA[0][0]);
}
template<int n, int m, int q> void lowLevelAvxStrassen(int* result, int* a, int* b, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += 9*StackAllocator::Allign(eN*eM) + 9* StackAllocator::Allign(eM*eQ) + 7* StackAllocator::Allign(eN*eQ);
    }
    StackAllocator allocator(expected);
    lowLevelAvxStrassen<n, m, q>(a, b, steps, result, allocator);
}

template<typename M1, typename M2>
auto lowLevelAvxStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        lowLevelAvxStrassen<M1::CRow(), M1::CCol(), M2::CCol()>(result.data(), a.data(), b.data(), steps);
        return result;
    } else {
        auto result = a.createNew(a.rowCount(), b.columnCount());
        lowLevelAvxStrassen(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), steps);
        return result;
    }
}
#ifndef STRASSEN_H
#define STRASSEN_H
#include "MatrixInterface.h"
#include "matrixOperators.h"
#include "MatrixExtendedFunctions.h"
#include "StackAllocator.h"
#include "avxMul.h"
#include "MatrixView.h"
#include "MatrixExtendedFunctions.h"
#include "naiveBasicOperations.h"
#include "blasMul.h"
#include <utility>
#include <cmath>
#include <optional>


double log(double base, double value) {
    return ::log(value) / ::log(base);
}
int staticPaddingNewSize(int size, int base, int steps) {
    if (steps <= 0) return size;
    int divisor = pow(base, steps);
    auto remainder = size % divisor;
    if (remainder == 0) return size;
    else                return size + divisor - remainder;
}
std::optional<std::array<int, 3>> staticPaddingNewSizes(int n, int m, int q, int baseN, int baseM, int baseQ, int steps) {
    std::array<int, 3> paddedSizes = {
        staticPaddingNewSize(n, baseN, steps),
        staticPaddingNewSize(m, baseM, steps),
        staticPaddingNewSize(q, baseQ, steps)
    };
    if (paddedSizes[0] > n || paddedSizes[1] > m || paddedSizes[2] > q) {
        return paddedSizes;
    } else {
        return std::nullopt;
    }
}

enum class BaseMulType {
    Naive,
    Avx,
    Blas
};

template<typename ValueType> constexpr BaseMulType AutomaticBaseMulType = BaseMulType::Naive;
#if defined(BLAS_IS_AVAILABLE)
template<> constexpr BaseMulType AutomaticBaseMulType<float> = BaseMulType::Blas;
template<> constexpr BaseMulType AutomaticBaseMulType<double> = BaseMulType::Blas;
template<> constexpr BaseMulType AutomaticBaseMulType<std::complex<float>> = BaseMulType::Blas;
template<> constexpr BaseMulType AutomaticBaseMulType<std::complex<double>> = BaseMulType::Blas;
#elif defined(AVX2_IS_AVAILABLE)
template<> constexpr BaseMulType AutomaticBaseMulType<float> = BaseMulType::Avx;
template<> constexpr BaseMulType AutomaticBaseMulType<double> = BaseMulType::Avx;
#endif
#if defined(AVX2_IS_AVAILABLE)
template<> constexpr BaseMulType AutomaticBaseMulType<int32_t> = BaseMulType::Avx;
#endif

enum class ArithmeticOperation {
    Assign,
    Add,
    Sub,
    AddAssign,
    SubAssign
};

template<ArithmeticOperation op> struct FundamentalOperation {
    static const auto Value = op;
};
template<> struct FundamentalOperation<ArithmeticOperation::AddAssign> {
    static const auto Value = ArithmeticOperation::Add;
};
template<> struct FundamentalOperation<ArithmeticOperation::SubAssign> {
    static const auto Value = ArithmeticOperation::Sub;
};

template<ArithmeticOperation op> struct IsAssignOperation {static const auto Value = false;};
template<> struct IsAssignOperation<ArithmeticOperation::Assign> {static const auto Value = true;};
template<> struct IsAssignOperation<ArithmeticOperation::AddAssign> {static const auto Value = true;};
template<> struct IsAssignOperation<ArithmeticOperation::SubAssign> {static const auto Value = true;};

template<ArithmeticOperation op, typename T1, typename T2> inline auto calculate(T1&& a, T2&& b) {
    if constexpr (op == ArithmeticOperation::Assign) {
        return a = b;
    } else if constexpr (op == ArithmeticOperation::Add) {
        return a + b;
    } else if constexpr (op == ArithmeticOperation::Sub) {
        return a - b;
    } else if constexpr (op == ArithmeticOperation::AddAssign) {
        return a += b;
    } else if constexpr (op == ArithmeticOperation::SubAssign) {
        return a -= b;
    }
}
template<ArithmeticOperation op, typename T1, typename T2, typename T3> inline auto calculate(T1&& dst, T2&& a, T3&& b) {
    if constexpr (IsAssignOperation<op>::Value) {
        return calculate<op>(std::forward<T1>(dst), calculate<FundamentalOperation<op>::Value>(std::forward<T2>(a), std::forward<T3>(b)));
    } else {
        return dst = calculate<op>(std::forward<T2>(a), std::forward<T3>(b));
    }
}
template<ArithmeticOperation op, typename T> void operation(T* dst, const T* a, const T* b, int n, int m) {
    static_assert(!IsAssignOperation<op>::Value);
    for (int i = 0; i < n*m; ++i) {
        calculate<op>(dst[i], a[i], b[i]);
    }
}
template<ArithmeticOperation op, typename T> void operation(T* dst, const T* src, int n, int m, int effDst, int effSrc) {
    static_assert(IsAssignOperation<op>::Value);
    for (int i = 0; i < n; ++i) {
        auto dstR = &dst[i*effDst];
        auto srcR = &src[i*effSrc];
        for (int j = 0; j < m; ++j) {
            calculate<op>(dstR[j], srcR[j]);
        }
    }
}
template<ArithmeticOperation op, typename T> void operation(T* dst, const T* a, const T* b, int n, int m, int effDst, int effA, int effB) {
    static_assert(!IsAssignOperation<op>::Value);
    for (int i = 0; i < n; ++i) {
        auto dstR = &dst[i*effDst];
        auto aR = &a[i*effA];
        auto bR = &b[i*effB];
        for (int j = 0; j < m; ++j) {
            calculate<op>(dstR[j], aR[j], bR[j]);
        }
    }
}
template<int n, int m, ArithmeticOperation op, typename T> void operation(T* dst, const T* a, const T* b) {
    static_assert(!IsAssignOperation<op>::Value);
    for (int i = 0; i < n*m; ++i) {
        calculate<op>(dst[i], a[i], b[i]);
    }
}
template<int n, int m, int effDst, int effSrc, ArithmeticOperation op, typename T> void operation(T* dst, const T* src) {
    static_assert(IsAssignOperation<op>::Value);
    for (int i = 0; i < n; ++i) {
        auto dstR = &dst[i*effDst];
        auto srcR = &src[i*effSrc];
        for (int j = 0; j < m; ++j) {
            calculate<op>(dstR[j], srcR[j]);
        }
    }
}
template<int n, int m, int effDst, int effA, int effB, ArithmeticOperation op, typename T> void operation(T* dst, const T* a, const T* b) {
    static_assert(!IsAssignOperation<op>::Value);
    for (int i = 0; i < n; ++i) {
        auto dstR = &dst[i*effDst];
        auto aR = &a[i*effA];
        auto bR = &b[i*effB];
        for (int j = 0; j < m; ++j) {
            calculate<op>(dstR[j], aR[j], bR[j]);
        }
    }
}

template<ArithmeticOperation op, typename T> T* operation(const T* a, const T* b, int n, int m, StackAllocator<T>& allocator) {
    T* result = allocator.alloc(n*m);
    operation<op>(result, a, b, n, m);
    return result;
}
template<int n, int m, ArithmeticOperation op, typename T> T* operation(const T* a, const T* b, StackAllocator<T>& allocator) {
    T* result = allocator.alloc(n*m);
    operation<n, m, op>(result, a, b);
    return result;
}
template<ArithmeticOperation op, typename T> void operation(MatrixView<T> dst, MatrixConstView<T> src) {
    return operation<op>(dst.data(), src.data(), dst.rowCount(), dst.columnCount(), dst.effectiveColumnCount(), src.effectiveColumnCount());
}
template<ArithmeticOperation op, typename T> void operation(MatrixView<T> dst, MatrixConstView<T> a, MatrixConstView<T> b) {
    return operation<op>(dst.data(), a.data(), b.data(), dst.rowCount(), dst.columnCount(), dst.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
}

template<typename T> void naiveMul(T* result, const T* a, const T* b, int n, int m, int q) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < q; ++j) {
            T sum{};
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * m] * b[j + k * q];
            }
            result[j + i * q] = sum;
        }
    }
}
template<int n, int m, int q, typename T> void naiveMul(T* result, const T* a, const T* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < q; ++j) {
            T sum{};
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * m] * b[j + k * q];
            }
            result[j + i * q] = sum;
        }
    }
}

template<BaseMulType opType, typename T>
void mul(T* result, const T* a, const T* b, int n, int m, int q) {
    if constexpr (opType == BaseMulType::Naive) {
        naiveMul(result, a, b, n, m, q);
    } else if constexpr (opType == BaseMulType::Avx) {
        avxMul7(result, a, b, n, m, q);
    } else if constexpr (opType == BaseMulType::Blas) {
        for (int i = 0; i < n*q; ++i) result[i] = 0;
        blas::mul(result, a, b, n, m, q);
    }
}

template<int n, int m, int q, BaseMulType opType, typename T>
void mul(T* result, const T* a, const T* b) {
    if constexpr (opType == BaseMulType::Naive) {
        naiveMul<n, m, q>(result, a, b);
    } else if constexpr (opType == BaseMulType::Avx) {
        avxMul7<n, m, q>(result, a, b);
    } else if constexpr (opType == BaseMulType::Blas) {
        for (int i = 0; i < n*q; ++i) result[i] = 0;
        blas::mul(result, a, b, n, m, q);
    }
}

template<BaseMulType opType, typename M1, typename M2> 
auto mul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (opType == BaseMulType::Naive) return naiveMul(a, b);
    if constexpr (opType == BaseMulType::Avx) return avxMul(a, b);
    if constexpr (opType == BaseMulType::Blas) return blasMul(a, b);
}

template<BaseMulType opType, typename M1, typename M2>
auto strassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if (steps <= 0) {
        return mul<opType>(a, b);
    }

    auto dA = matrixDivide<2, 2>(a);
    auto dB = matrixDivide<2, 2>(b);

    auto m1 = strassen<opType>(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1], steps - 1);
    auto m2 = strassen<opType>(dA[1][0] + dA[1][1], dB[0][0], steps - 1);
    auto m3 = strassen<opType>(dA[0][0], dB[0][1] - dB[1][1], steps - 1);
    auto m4 = strassen<opType>(dA[1][1], dB[1][0] - dB[0][0], steps - 1);
    auto m5 = strassen<opType>(dA[0][0] + dA[0][1], dB[1][1], steps - 1);
    auto m6 = strassen<opType>(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1], steps - 1);
    auto m7 = strassen<opType>(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1], steps - 1);

    Matrix<typename M1::ValueType> c(a.rowCount(), b.columnCount());
    auto dC = matrixDivideView<2, 2>(c);
    dC[0][0].copy(m1 + m4 - m5 + m7);
    dC[0][1].copy(m3 + m5);
    dC[1][0].copy(m2 + m4);
    dC[1][1].copy(m1 + m3 - m2 + m6);

    return c;
}
template<BaseMulType opType, typename M1, typename M2>
auto strassenWithStaticPadding(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    auto paddedSizesOpt = staticPaddingNewSizes(a.rowCount(), a.columnCount(), b.columnCount(), 2, 2, 2, steps);
    if (paddedSizesOpt) {
        auto& paddedSizes = *paddedSizesOpt;
        auto newA = a.createNew(paddedSizes[0], paddedSizes[1]);
        auto newB = b.createNew(paddedSizes[1], paddedSizes[2]);
        newA.copy(a);
        newB.copy(b);
        auto result = strassen<opType>(newA, newB, steps);
        return shrink(result, a.rowCount(), b.columnCount());
    } else {
        return strassen<opType>(a, b, steps);
    }
}
template<typename M1, typename M2>
auto strassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return strassenWithStaticPadding<BaseMulType::Naive>(a, b, steps);
}
template<typename M1, typename M2>
auto strassenAvx(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if constexpr (AVX256::IsAvailable) return strassenWithStaticPadding<BaseMulType::Avx>(a, b, steps);
    else static_assert(AVX256::IsAvailable && always_false_v<M1>, AVX256_StaticAssertMessage);
}
template<typename M1, typename M2>
auto strassenAuto(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return strassenWithStaticPadding<AutomaticBaseMulType<typename M1::ValueType>>(a, b, steps);
}

template<BaseMulType opType, int Steps, typename M1, typename M2> auto strassenImpl(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, std::false_type) {
    return mul<opType>(a, b);
}
template<BaseMulType opType, int Steps, typename M1, typename M2> auto strassenImpl(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, std::true_type unused = std::true_type()) {
    if (Steps <= 0) {
        return strassenImpl<opType, Steps>(a, b, std::false_type{});
    }

    auto dA = matrixDivide<2, 2>(a);
    auto dB = matrixDivide<2, 2>(b);

    auto m1 = strassenImpl<opType, Steps - 1>(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m2 = strassenImpl<opType, Steps - 1>(dA[1][0] + dA[1][1], dB[0][0], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m3 = strassenImpl<opType, Steps - 1>(dA[0][0], dB[0][1] - dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m4 = strassenImpl<opType, Steps - 1>(dA[1][1], dB[1][0] - dB[0][0], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m5 = strassenImpl<opType, Steps - 1>(dA[0][0] + dA[0][1], dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m6 = strassenImpl<opType, Steps - 1>(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    auto m7 = strassenImpl<opType, Steps - 1>(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1], typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());

    Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> c;
    auto dC = matrixDivideView<2, 2>(c);
    dC[0][0].copy(m1 + m4 - m5 + m7);
    dC[0][1].copy(m3 + m5);
    dC[1][0].copy(m2 + m4);
    dC[1][1].copy(m1 + m3 - m2 + m6);

    return c;
}
template<int Steps, typename M1, typename M2> auto strassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return strassenImpl<BaseMulType::Naive, Steps>(a, b);
}
template<int Steps, typename M1, typename M2> auto strassenAvx(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (AVX256::IsAvailable) return strassenImpl<BaseMulType::Avx, Steps>(a, b);
    else static_assert(AVX256::IsAvailable && always_false_v<M1>, AVX256_StaticAssertMessage);
}
template<int Steps, typename M1, typename M2> auto strassenAuto(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return strassenImpl<AutomaticBaseMulType<typename M1::ValueType>, Steps>(a, b);
}

template<int Rows, int Columns, typename T>
std::array<std::array<T*, Columns>, Rows> strassenDivide(T* a, int n, int m, StackAllocator<T>& allocator) {
    std::array<std::array<T*, Columns>, Rows> result;
    auto rowSize = n / Rows;
    auto columnSize = m / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[y][x] = allocator.alloc(rowSize*columnSize);
            operation<ArithmeticOperation::Assign>(result[y][x], a + x * columnSize + y * rowSize*m, rowSize, columnSize, columnSize, m);
        }
    }
    return result;
}
template<int Rows, int Columns, typename T>
std::array<std::array<T*, Columns>, Rows> strassenDivideView(T* a, int n, int m, int effM) {
    std::array<std::array<T*, Columns>, Rows> result;
    auto rowSize = n / Rows;
    auto columnSize = m / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[y][x] = a + x * columnSize + y * rowSize*effM;
        }
    }
    return result;
}
template<int Rows, int Columns, typename T>
std::array<std::array<T*, Columns>, Rows> strassenDivideView(T* a, int n, int m) {
    return strassenDivideView<Rows, Columns>(a, n, m, m);
}


template<BaseMulType opType, typename T>
void lowLevelStrassen(T* a, T* b, int n, int m, int q, int steps, T* c, StackAllocator<T>& allocator) {
    if (steps <= 0) {
        mul<opType>(c, a, b, n, m, q);
        return;
    }

    auto dA = strassenDivide<2, 2>(a, n, m, allocator);
    auto dB = strassenDivide<2, 2>(b, m, q, allocator);

    int halfN = n / 2;
    int halfM = m / 2;
    int halfQ = q / 2;

    auto m1 = allocator.alloc(halfN * halfQ);
    auto m1_a = operation<ArithmeticOperation::Add>(dA[0][0], dA[1][1], halfN, halfM, allocator);
    auto m1_b = operation<ArithmeticOperation::Add>(dB[0][0], dB[1][1], halfM, halfQ, allocator);
    lowLevelStrassen<opType>(m1_a, m1_b, halfN, halfM, halfQ, steps - 1, m1, allocator);
    allocator.dealloc(m1_b, halfM*halfQ);
    allocator.dealloc(m1_a, halfN*halfM);

    auto m2 = allocator.alloc(halfN * halfQ);
    auto m2_a = operation<ArithmeticOperation::Add>(dA[1][0], dA[1][1], halfN, halfM, allocator);
    lowLevelStrassen<opType>(m2_a, dB[0][0], halfN, halfM, halfQ, steps - 1, m2, allocator);
    allocator.dealloc(m2_a, halfN*halfM);

    auto m3 = allocator.alloc(halfN * halfQ);
    auto m3_b = operation<ArithmeticOperation::Sub>(dB[0][1], dB[1][1], halfM, halfQ, allocator);
    lowLevelStrassen<opType>(dA[0][0], m3_b, halfN, halfM, halfQ, steps - 1, m3, allocator);
    allocator.dealloc(m3_b, halfM*halfQ);

    auto m4 = allocator.alloc(halfN * halfQ);
    auto m4_b = operation<ArithmeticOperation::Sub>(dB[1][0], dB[0][0], halfM, halfQ, allocator);
    lowLevelStrassen<opType>(dA[1][1], m4_b, halfN, halfM, halfQ, steps - 1, m4, allocator);
    allocator.dealloc(m4_b, halfM*halfQ);

    auto m5 = allocator.alloc(halfN * halfQ);
    auto m5_a = operation<ArithmeticOperation::Add>(dA[0][0], dA[0][1], halfN, halfM, allocator);
    lowLevelStrassen<opType>(m5_a, dB[1][1], halfN, halfM, halfQ, steps - 1, m5, allocator);
    allocator.dealloc(m5_a, halfN*halfM);

    auto m6 = allocator.alloc(halfN * halfQ);
    auto m6_a = operation<ArithmeticOperation::Sub>(dA[1][0], dA[0][0], halfN, halfM, allocator);
    auto m6_b = operation<ArithmeticOperation::Add>(dB[0][0], dB[0][1], halfM, halfQ, allocator);
    lowLevelStrassen<opType>(m6_a, m6_b, halfN, halfM, halfQ, steps - 1, m6, allocator);
    allocator.dealloc(m6_b, halfM*halfQ);
    allocator.dealloc(m6_a, halfN*halfM);

    auto m7 = allocator.alloc(halfN * halfQ);
    auto m7_a = operation<ArithmeticOperation::Sub>(dA[0][1], dA[1][1], halfN, halfM, allocator);
    auto m7_b = operation<ArithmeticOperation::Add>(dB[1][0], dB[1][1], halfM, halfQ, allocator);
    lowLevelStrassen<opType>(m7_a, m7_b, halfN, halfM, halfQ, steps - 1, m7, allocator);
    allocator.dealloc(m7_b, halfM*halfQ);
    allocator.dealloc(m7_a, halfN*halfM);

    auto dC = strassenDivideView<2, 2>(c, n, q);
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

    allocator.dealloc(m7, halfN*halfQ);
    allocator.dealloc(m6, halfN*halfQ);
    allocator.dealloc(m5, halfN*halfQ);
    allocator.dealloc(m4, halfN*halfQ);
    allocator.dealloc(m3, halfN*halfQ);
    allocator.dealloc(m2, halfN*halfQ);
    allocator.dealloc(m1, halfN*halfQ);
    allocator.dealloc(dB[1][1], halfM*halfQ);
    allocator.dealloc(dB[1][0], halfM*halfQ);
    allocator.dealloc(dB[0][1], halfM*halfQ);
    allocator.dealloc(dB[0][0], halfM*halfQ);
    allocator.dealloc(dA[1][1], halfN*halfM);
    allocator.dealloc(dA[1][0], halfN*halfM);
    allocator.dealloc(dA[0][1], halfN*halfM);
    allocator.dealloc(dA[0][0], halfN*halfM);
}

template<BaseMulType opType, typename T>
void lowLevelStrassenWithStaticPadding(T* result, T* a, T* b, int n, int m, int q, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    
    auto paddedSizesOpt = staticPaddingNewSizes(n, m, q, 2, 2, 2, steps);
    if (paddedSizesOpt) {
        auto& paddedSizes = *paddedSizesOpt;
        eN = paddedSizes[0];
        eM = paddedSizes[1];
        eQ = paddedSizes[2];
        expected += StackAllocator<T>::Allign(eN*eM);
        expected += StackAllocator<T>::Allign(eM*eQ);
        expected += StackAllocator<T>::Allign(eN*eQ);
    }

    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += 5*StackAllocator<T>::Allign(eN*eM) + 5*StackAllocator<T>::Allign(eM*eQ) + 7*StackAllocator<T>::Allign(eN*eQ);
    }

    StackAllocator<T> allocator(expected);

    if (paddedSizesOpt) {
        auto& paddedSizes = *paddedSizesOpt;
        auto newA = allocator.alloc(paddedSizes[0] * paddedSizes[1]);
        auto newB = allocator.alloc(paddedSizes[1] * paddedSizes[2]);
        auto newResult = allocator.alloc(paddedSizes[0] * paddedSizes[2]);
        for (int i = 0; i < paddedSizes[0] * paddedSizes[1]; ++i) {
            newA[i] = T{};
        }
        for (int i = 0; i < paddedSizes[1] * paddedSizes[2]; ++i) {
            newB[i] = T{};
        }
        operation<ArithmeticOperation::Assign>(newA, a, n, m, paddedSizes[1], m);
        operation<ArithmeticOperation::Assign>(newB, b, m, q, paddedSizes[2], q);
        lowLevelStrassen<opType>(newA, newB, paddedSizes[0], paddedSizes[1], paddedSizes[2], steps, newResult, allocator);
        operation<ArithmeticOperation::Assign>(result, newResult, n, q, q, paddedSizes[2]);
        allocator.dealloc(newResult, paddedSizes[0] * paddedSizes[2]);
        allocator.dealloc(newB, paddedSizes[1] * paddedSizes[2]);
        allocator.dealloc(newA, paddedSizes[0] * paddedSizes[1]);
    } else {
        lowLevelStrassen<opType>(a, b, n, m, q, steps, result, allocator);
    }
}

template<int n, int m, int q, BaseMulType opType, typename T>
void lowLevelStrassen(T* a, T* b, int steps, T* c, StackAllocator<T>& allocator) {
    if (steps <= 0) {
        mul<n, m, q, opType>(c, a, b);
        return;
    }

    auto dA = strassenDivide<2, 2>(a, n, m, allocator);
    auto dB = strassenDivide<2, 2>(b, m, q, allocator);

    constexpr int hn = n / 2;
    constexpr int hm = m / 2;
    constexpr int hq = q / 2;

    auto m1 = allocator.alloc(hn * hq);
    auto m1_a = operation<hn, hm, ArithmeticOperation::Add>(dA[0][0], dA[1][1], allocator);
    auto m1_b = operation<hm, hq, ArithmeticOperation::Add>(dB[0][0], dB[1][1], allocator);
    lowLevelStrassen<hn, hm, hq, opType>(m1_a, m1_b, steps - 1, m1, allocator);
    allocator.dealloc(m1_b, hm*hq);
    allocator.dealloc(m1_a, hn*hm);

    auto m2 = allocator.alloc(hn * hq);
    auto m2_a = operation<hn, hm, ArithmeticOperation::Add>(dA[1][0], dA[1][1], allocator);
    lowLevelStrassen<hn, hm, hq, opType>(m2_a, dB[0][0], steps - 1, m2, allocator);
    allocator.dealloc(m2_a, hn*hm);

    auto m3 = allocator.alloc(hn * hq);
    auto m3_b = operation<hm, hq, ArithmeticOperation::Sub>(dB[0][1], dB[1][1], allocator);
    lowLevelStrassen<hn, hm, hq, opType>(dA[0][0], m3_b, steps - 1, m3, allocator);
    allocator.dealloc(m3_b, hm*hq);

    auto m4 = allocator.alloc(hn * hq);
    auto m4_b = operation<hm, hq, ArithmeticOperation::Sub>(dB[1][0], dB[0][0], allocator);
    lowLevelStrassen<hn, hm, hq, opType>(dA[1][1], m4_b, steps - 1, m4, allocator);
    allocator.dealloc(m4_b, hm*hq);

    auto m5 = allocator.alloc(hn * hq);
    auto m5_a = operation<hn, hm, ArithmeticOperation::Add>(dA[0][0], dA[0][1], allocator);
    lowLevelStrassen<hn, hm, hq, opType>(m5_a, dB[1][1], steps - 1, m5, allocator);
    allocator.dealloc(m5_a, hn*hm);

    auto m6 = allocator.alloc(hn * hq);
    auto m6_a = operation<hn, hm, ArithmeticOperation::Sub>(dA[1][0], dA[0][0], allocator);
    auto m6_b = operation<hm, hq, ArithmeticOperation::Add>(dB[0][0], dB[0][1], allocator);
    lowLevelStrassen<hn, hm, hq, opType>(m6_a, m6_b, steps - 1, m6, allocator);
    allocator.dealloc(m6_b, hm*hq);
    allocator.dealloc(m6_a, hn*hm);

    auto m7 = allocator.alloc(hn * hq);
    auto m7_a = operation<hn, hm, ArithmeticOperation::Sub>(dA[0][1], dA[1][1], allocator);
    auto m7_b = operation<hm, hq, ArithmeticOperation::Add>(dB[1][0], dB[1][1], allocator);
    lowLevelStrassen<hn, hm, hq, opType>(m7_a, m7_b, steps - 1, m7, allocator);
    allocator.dealloc(m7_b, hm*hq);
    allocator.dealloc(m7_a, hn*hm);

    auto dC = strassenDivideView<2, 2>(c, n, q);
    for (int i = 0; i < hn; ++i) {
        for (int j = 0; j < hq; ++j) {
            int dstIndex = j + i * q;
            int index = j + i * hq;
            dC[0][0][dstIndex] = m1[index] + m4[index] - m5[index] + m7[index];
            dC[0][1][dstIndex] = m3[index] + m5[index];
            dC[1][0][dstIndex] = m2[index] + m4[index];
            dC[1][1][dstIndex] = m1[index] + m3[index] - m2[index] + m6[index];
        }
    }

    allocator.dealloc(m7, hn*hq);
    allocator.dealloc(m6, hn*hq);
    allocator.dealloc(m5, hn*hq);
    allocator.dealloc(m4, hn*hq);
    allocator.dealloc(m3, hn*hq);
    allocator.dealloc(m2, hn*hq);
    allocator.dealloc(m1, hn*hq);
    allocator.dealloc(dB[1][1], hm*hq);
    allocator.dealloc(dB[1][0], hm*hq);
    allocator.dealloc(dB[0][1], hm*hq);
    allocator.dealloc(dB[0][0], hm*hq);
    allocator.dealloc(dA[1][1], hn*hm);
    allocator.dealloc(dA[1][0], hn*hm);
    allocator.dealloc(dA[0][1], hn*hm);
    allocator.dealloc(dA[0][0], hn*hm);
}
template<int n, int m, int q, BaseMulType opType, typename T>
void lowLevelStrassen(T* result, T* a, T* b, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += 5*StackAllocator<T>::Allign(eN*eM) + 5*StackAllocator<T>::Allign(eM*eQ) + 7*StackAllocator<T>::Allign(eN*eQ);
    }
    StackAllocator<T> allocator(expected);
    lowLevelStrassen<n, m, q, opType>(a, b, steps, result, allocator);
}

template<BaseMulType opType, typename M1, typename M2>
auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        lowLevelStrassen<M1::CRow(), M1::CCol(), M2::CCol(), opType>(result.data(), a.data(), b.data(), steps);
        return result;
    } else {
        auto result = a.createNew(a.rowCount(), b.columnCount());
        lowLevelStrassenWithStaticPadding<opType>(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), steps);
        return result;
    }
}
template<typename M1, typename M2>
auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return lowLevelStrassen<BaseMulType::Naive>(a, b, steps);
}
template<typename M1, typename M2>
auto lowLevelAvxStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if constexpr (AVX256::IsAvailable) return lowLevelStrassen<BaseMulType::Avx>(a, b, steps);
    else static_assert(AVX256::IsAvailable && always_false_v<M1>, AVX256_StaticAssertMessage);
}
template<typename M1, typename M2>
auto lowLevelAutoStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return lowLevelStrassen<AutomaticBaseMulType<typename M1::ValueType>>(a, b, steps);
}

template<BaseMulType opType, typename T>
void minSpaceStrassenMul(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, StackAllocator<T>& allocator) {
    T* cx = c;
    T* ax = a;
    T* bx = b;
    if (effA != m) {
        ax = allocator.alloc(n*m);
        operation<ArithmeticOperation::Assign>(ax, a, n, m, m, effA);
    }
    if (effB != q) {
        bx = allocator.alloc(m*q);
        operation<ArithmeticOperation::Assign>(bx, b, m, q, q, effB);
    }
    if (effC != q) {
        cx = allocator.alloc(n*q);
    }

    mul<opType>(cx, ax, bx, n, m, q);

    if (effC != q) {
        operation<ArithmeticOperation::Assign>(c, cx, n, q, effC, q);
        allocator.dealloc(cx, n*q);
    }
    if (effB != q) allocator.dealloc(bx, m*q);
    if (effA != m) allocator.dealloc(ax, n*m);
}

template<BaseMulType opType, typename T>
void minSpaceStrassen(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
    if (steps <= 0) {
        minSpaceStrassenMul<opType>(c, a, b, n, m, q, effC, effA, effB, allocator);
        return;
    }

    auto dA = strassenDivideView<2, 2>(a, n, m, effA);
    auto dB = strassenDivideView<2, 2>(b, m, q, effB);
    auto dC = strassenDivideView<2, 2>(c, n, q, effC);
    
    int hn = n / 2;
    int hm = m / 2;
    int hq = q / 2;

    auto tempA = allocator.alloc(hn * hm);
    auto tempB = allocator.alloc(hm * hq);
    auto tempC = allocator.alloc(hn * hq);

    operation<ArithmeticOperation::Add>(tempA, dA[0][0], dA[1][1], hn, hm, hm, effA, effA);
    operation<ArithmeticOperation::Add>(tempB, dB[0][0], dB[1][1], hm, hq, hq, effB, effB);
    minSpaceStrassen<opType>(dC[0][0], tempA, tempB, hn, hm, hq, effC, hm, hq, steps - 1, allocator);
    operation<ArithmeticOperation::Assign>(dC[1][1], dC[0][0], hn, hq, effC, effC);

    operation<ArithmeticOperation::Add>(tempA, dA[1][0], dA[1][1], hn, hm, hm, effA, effA);
    minSpaceStrassen<opType>(dC[1][0], tempA, dB[0][0], hn, hm, hq, effC, hm, effB, steps - 1, allocator);
    operation<ArithmeticOperation::SubAssign>(dC[1][1], dC[1][0], hn, hq, effC, effC);

    operation<ArithmeticOperation::Sub>(tempB, dB[0][1], dB[1][1], hm, hq, hq, effB, effB);
    minSpaceStrassen<opType>(dC[0][1], dA[0][0], tempB, hn, hm, hq, effC, effA, hq, steps - 1, allocator);
    operation<ArithmeticOperation::AddAssign>(dC[1][1], dC[0][1], hn, hq, effC, effC);

    operation<ArithmeticOperation::Sub>(tempB, dB[1][0], dB[0][0], hm, hq, hq, effB, effB);
    minSpaceStrassen<opType>(tempC, dA[1][1], tempB, hn, hm, hq, hq, effA, hq, steps - 1, allocator);
    operation<ArithmeticOperation::AddAssign>(dC[0][0], tempC, hn, hq, effC, hq);
    operation<ArithmeticOperation::AddAssign>(dC[1][0], tempC, hn, hq, effC, hq);

    operation<ArithmeticOperation::Add>(tempA, dA[0][0], dA[0][1], hn, hm, hm, effA, effA);
    minSpaceStrassen<opType>(tempC, tempA, dB[1][1], hn, hm, hq, hq, hm, effB, steps - 1, allocator);
    operation<ArithmeticOperation::SubAssign>(dC[0][0], tempC, hn, hq, effC, hq);
    operation<ArithmeticOperation::AddAssign>(dC[0][1], tempC, hn, hq, effC, hq);

    operation<ArithmeticOperation::Sub>(tempA, dA[1][0], dA[0][0], hn, hm, hm, effA, effA);
    operation<ArithmeticOperation::Add>(tempB, dB[0][0], dB[0][1], hm, hq, hq, effB, effB);
    minSpaceStrassen<opType>(tempC, tempA, tempB, hn, hm, hq, hq, hm, hq, steps - 1, allocator);
    operation<ArithmeticOperation::AddAssign>(dC[1][1], tempC, hn, hq, effC, hq);

    operation<ArithmeticOperation::Sub>(tempA, dA[0][1], dA[1][1], hn, hm, hm, effA, effA);
    operation<ArithmeticOperation::Add>(tempB, dB[1][0], dB[1][1], hm, hq, hq, effB, effB);
    minSpaceStrassen<opType>(tempC, tempA, tempB, hn, hm, hq, hq, hm, hq, steps - 1, allocator);
    operation<ArithmeticOperation::AddAssign>(dC[0][0], tempC, hn, hq, effC, hq);

    allocator.dealloc(tempC, hn*hq);
    allocator.dealloc(tempB, hm*hq);
    allocator.dealloc(tempA, hn*hm);
}


template<BaseMulType opType, typename T> void minSpaceStrassenWithStaticPadding(T* c, T* a, T* b, int n, int m, int q, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    
    auto paddedSizesOpt = staticPaddingNewSizes(n, m, q, 2, 2, 2, steps);
    if (paddedSizesOpt) {
        auto& paddedSizes = *paddedSizesOpt;
        eN = paddedSizes[0];
        eM = paddedSizes[1];
        eQ = paddedSizes[2];
        expected += StackAllocator<T>::Allign(eN*eM);
        expected += StackAllocator<T>::Allign(eM*eQ);
        expected += StackAllocator<T>::Allign(eN*eQ);
    }

    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += StackAllocator<T>::Allign(eN*eM) + StackAllocator<T>::Allign(eM*eQ) + StackAllocator<T>::Allign(eN*eQ);
    }
    if (steps >= 1) {
        expected += std::max(StackAllocator<T>::Allign(eN*eM), StackAllocator<T>::Allign(eM*eQ)) + StackAllocator<T>::Allign(eN*eQ);
    }

    StackAllocator<T> allocator(expected);

    if (paddedSizesOpt) {
        auto& paddedSizes = *paddedSizesOpt;
        auto newA = allocator.alloc(paddedSizes[0] * paddedSizes[1]);
        auto newB = allocator.alloc(paddedSizes[1] * paddedSizes[2]);
        auto newC = allocator.alloc(paddedSizes[0] * paddedSizes[2]);
        for (int i = 0; i < paddedSizes[0] * paddedSizes[1]; ++i) {
            newA[i] = T{};
        }
        for (int i = 0; i < paddedSizes[1] * paddedSizes[2]; ++i) {
            newB[i] = T{};
        }
        operation<ArithmeticOperation::Assign>(newA, a, n, m, paddedSizes[1], m);
        operation<ArithmeticOperation::Assign>(newB, b, m, q, paddedSizes[2], q);
        minSpaceStrassen<opType>(newC, newA, newB, paddedSizes[0], paddedSizes[1], paddedSizes[2], paddedSizes[2], paddedSizes[1], paddedSizes[2], steps, allocator);
        operation<ArithmeticOperation::Assign>(c, newC, n, q, q, paddedSizes[2]);
        allocator.dealloc(newC, paddedSizes[0] * paddedSizes[2]);
        allocator.dealloc(newB, paddedSizes[1] * paddedSizes[2]);
        allocator.dealloc(newA, paddedSizes[0] * paddedSizes[1]);
    } else {
        minSpaceStrassen<opType>(c, a, b, n, m, q, q, m, q, steps, allocator);
    }
}

template<int n, int m, int q, int effC, int effA, int effB, BaseMulType opType, typename T>
void minSpaceStrassenMul(T* c, T* a, T* b, StackAllocator<T>& allocator) {
    T* cx = c;
    T* ax = a;
    T* bx = b;
    if (effA != m) {
        ax = allocator.alloc(n*m);
        operation<n, m, m, effA, ArithmeticOperation::Assign>(ax, a);
    }
    if (effB != q) {
        bx = allocator.alloc(m*q);
        operation<m, q, q, effB, ArithmeticOperation::Assign>(bx, b);
    }
    if (effC != q) {
        cx = allocator.alloc(n*q);
    }

    mul<opType>(cx, ax, bx, n, m, q);

    if (effC != q) {
        operation<n, q, effC, q, ArithmeticOperation::Assign>(c, cx);
        allocator.dealloc(cx, n*q);
    }
    if (effB != q) allocator.dealloc(bx, m*q);
    if (effA != m) allocator.dealloc(ax, n*m);
}


template<int n, int m, int q, int effC, int effA, int effB, int Steps, BaseMulType opType, typename T>
void minSpaceStrassen(T* c, T* a, T* b, StackAllocator<T>& allocator, std::false_type) {
    minSpaceStrassenMul<n, m, q, effC, effA, effB, opType>(c, a, b, allocator);
}
template<int n, int m, int q, int effC, int effA, int effB, int Steps, BaseMulType opType, typename T>
void minSpaceStrassen(T* c, T* a, T* b, StackAllocator<T>& allocator, std::true_type unused = std::true_type()) {
    if constexpr (Steps <= 0) {
        minSpaceStrassen<n, m, q, effC, effA, effB, Steps, opType>(c, a, b, allocator, std::false_type{});
        return;
    }

    auto dA = strassenDivideView<2, 2>(a, n, m, effA);
    auto dB = strassenDivideView<2, 2>(b, m, q, effB);
    auto dC = strassenDivideView<2, 2>(c, n, q, effC);

    constexpr int hn = n / 2;
    constexpr int hm = m / 2;
    constexpr int hq = q / 2;

    auto tempA = allocator.alloc(hn * hm);
    auto tempB = allocator.alloc(hm * hq);
    auto tempC = allocator.alloc(hn * hq);

    operation<hn, hm, hm, effA, effA, ArithmeticOperation::Add>(tempA, dA[0][0], dA[1][1]);
    operation<hm, hq, hq, effB, effB, ArithmeticOperation::Add>(tempB, dB[0][0], dB[1][1]);
    minSpaceStrassen<hn, hm, hq, effC, hm, hq, Steps - 1, opType>(dC[0][0], tempA, tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operation<hn, hq, effC, effC, ArithmeticOperation::SubAssign>(dC[1][1], dC[0][0]);

    operation<hn, hm, hm, effA, effA, ArithmeticOperation::Add>(tempA, dA[1][0], dA[1][1]);
    minSpaceStrassen<hn, hm, hq, effC, hm, effB, Steps - 1, opType>(dC[1][0], tempA, dB[0][0], allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operation<hn, hq, effC, effC, ArithmeticOperation::SubAssign>(dC[1][1], dC[1][0]);

    operation<hm, hq, hq, effB, effB, ArithmeticOperation::Sub>(tempB, dB[0][1], dB[1][1]);
    minSpaceStrassen<hn, hm, hq, effC, effA, hq, Steps - 1, opType>(dC[0][1], dA[0][0], tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operation<hn, hq, effC, effC, ArithmeticOperation::AddAssign>(dC[1][1], dC[0][1]);

    operation<hm, hq, hq, effB, effB, ArithmeticOperation::Sub>(tempB, dB[1][0], dB[0][0]);
    minSpaceStrassen<hn, hm, hq, hq, effA, hq, Steps - 1, opType>(tempC, dA[1][1], tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operation<hn, hq, effC, hq, ArithmeticOperation::AddAssign>(dC[0][0], tempC);
    operation<hn, hq, effC, hq, ArithmeticOperation::AddAssign>(dC[1][0], tempC);

    operation<hn, hm, hm, effA, effA, ArithmeticOperation::Add>(tempA, dA[0][0], dA[0][1]);
    minSpaceStrassen<hn, hm, hq, hq, hm, effB, Steps - 1, opType>(tempC, tempA, dB[1][1], allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operation<hn, hq, effC, hq, ArithmeticOperation::SubAssign>(dC[0][0], tempC);
    operation<hn, hq, effC, hq, ArithmeticOperation::AddAssign>(dC[0][1], tempC);

    operation<hn, hm, hm, effA, effA, ArithmeticOperation::Sub>(tempA, dA[1][0], dA[0][0]);
    operation<hm, hq, hq, effB, effB, ArithmeticOperation::Add>(tempB, dB[0][0], dB[0][1]);
    minSpaceStrassen< hn, hm, hq, hq, hm, hq, Steps - 1, opType>(tempC, tempA, tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operation<hn, hq, effC, hq, ArithmeticOperation::AddAssign>(dC[1][1], tempC);

    operation<hn, hm, hm, effA, effA, ArithmeticOperation::Sub>(tempA, dA[0][1], dA[1][1]);
    operation<hm, hq, hq, effB, effB, ArithmeticOperation::Add>(tempB, dB[1][0], dB[1][1]);
    minSpaceStrassen<hn, hm, hq, hq, hm, hq, Steps - 1, opType>(tempC, tempA, tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operation<hn, hq, effC, hq, ArithmeticOperation::AddAssign>(dC[0][0], tempC);

    allocator.dealloc(tempC, hn*hq);
    allocator.dealloc(tempB, hm*hq);
    allocator.dealloc(tempA, hn*hm);
}


template<int n, int m, int q, int Steps, BaseMulType opType, typename T> void minSpaceStrassen(T* c, T* a, T* b) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    for (int i = 0; i < Steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += StackAllocator<T>::Allign(eN*eM) + StackAllocator<T>::Allign(eM*eQ) + StackAllocator<T>::Allign(eN*eQ);
    }
    if (Steps >= 1) {
        expected += std::max(StackAllocator<T>::Allign(eN*eM), StackAllocator<T>::Allign(eM*eQ)) + StackAllocator<T>::Allign(eN*eQ);
    }
    StackAllocator<T> allocator(expected);
    minSpaceStrassen<n, m, q, q, m, q, Steps, opType>(c, a, b, allocator);
}

template<BaseMulType opType, typename M1, typename M2>
auto minSpaceStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    auto result = a.createNew(a.rowCount(), b.columnCount());
    minSpaceStrassenWithStaticPadding<opType>(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), steps);
    return result;
}
template<int Steps, BaseMulType opType, typename M1, typename M2>
auto minSpaceStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        minSpaceStrassen<M1::CRow(), M1::CCol(), M2::CCol(), Steps, opType>(result.data(), a.data(), b.data());
        return result;
    } else {
        return minSpaceStrassen<opType>(a, b, Steps);
    }
}
template<typename M1, typename M2>
auto minSpaceStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return minSpaceStrassen<BaseMulType::Naive>(a, b, steps);
}
template<typename M1, typename M2>
auto minSpaceAvxStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if constexpr (AVX256::IsAvailable) return minSpaceStrassen<BaseMulType::Avx>(a, b, steps);
    else static_assert(AVX256::IsAvailable && always_false_v<M1>, AVX256_StaticAssertMessage);
}
template<typename M1, typename M2>
auto minSpaceAutoStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return minSpaceStrassen<AutomaticBaseMulType<typename M1::ValueType>>(a, b, steps);
}
template<int Steps, typename M1, typename M2>
auto minSpaceStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return minSpaceStrassen<Steps, BaseMulType::Naive>(a, b);
}
template<int Steps, typename M1, typename M2>
auto minSpaceAvxStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (AVX256::IsAvailable) return minSpaceStrassen<Steps, BaseMulType::Avx>(a, b);
    else static_assert(AVX256::IsAvailable && always_false_v<M1>, AVX256_StaticAssertMessage);
}
template<int Steps, typename M1, typename M2>
auto minSpaceAutoStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return minSpaceStrassen<Steps, AutomaticBaseMulType<typename M1::ValueType>>(a, b);
}

#endif

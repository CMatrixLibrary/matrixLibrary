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
#include "ThreadPool.h"
#include "genericArithmeticOperations.h"
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
    ParallelAvx,
    Blas,
    NaiveEff,
    AvxEff
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

template<ArithmeticOperation::OpType... ops, typename T, typename... Ts> T* operationWithAlloc(StackAllocator<T>& allocator, int n, int m, const T* arg, Ts&&... args) {
    auto result = allocator.alloc(n*m);
    operation<ops...>(n, m, result, arg, args...);
    return result;
}

template<int n, int m, ArithmeticOperation::OpType... ops, typename T, typename... Ts> T* operationWithAlloc(StackAllocator<T>& allocator, const T* arg, Ts&&... args) {
    auto result = allocator.alloc(n*m);
    operation<n, m, ops...>(result, arg, args...);
    return result;
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
        avx::mul(result, a, b, n, m, q);
    } else if constexpr (opType == BaseMulType::ParallelAvx) {
        avx::parallelMul(result, a, b, n, m, q);
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
        avx::mul<n, m, q, q, m, q>(result, a, b);
    } else if constexpr (opType == BaseMulType::ParallelAvx) {
        avx::parallelMul<n, m, q, q, m, q>(result, a, b);
    } else if constexpr (opType == BaseMulType::Blas) {
        for (int i = 0; i < n*q; ++i) result[i] = 0;
        blas::mul(result, a, b, n, m, q);
    }
}

template<BaseMulType opType, typename M1, typename M2> 
auto mul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (opType == BaseMulType::Naive) return naiveMul(a, b);
    if constexpr (opType == BaseMulType::Avx) return avx::mul(a, b);
    if constexpr (opType == BaseMulType::ParallelAvx) return avx::parallelMul(a, b);
    if constexpr (opType == BaseMulType::Blas) return blas::mul(a, b);
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
auto strassenPool(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if (steps <= 0) {
        return mul<opType>(a, b);
    }

    auto dA = matrixDivide<2, 2>(a);
    auto dB = matrixDivide<2, 2>(b);

    Matrix<typename M1::ValueType> m1, m2, m3, m4, m5, m6, m7;
    ThreadPool p;
    p.addTask([&]() {m1 = strassen<opType>(dA[0][0] + dA[1][1], dB[0][0] + dB[1][1], steps - 1); });
    p.addTask([&]() {m2 = strassen<opType>(dA[1][0] + dA[1][1], dB[0][0], steps - 1); });
    p.addTask([&]() {m3 = strassen<opType>(dA[0][0], dB[0][1] - dB[1][1], steps - 1); });
    p.addTask([&]() {m4 = strassen<opType>(dA[1][1], dB[1][0] - dB[0][0], steps - 1); });
    p.addTask([&]() {m5 = strassen<opType>(dA[0][0] + dA[0][1], dB[1][1], steps - 1); });
    p.addTask([&]() {m6 = strassen<opType>(dA[1][0] - dA[0][0], dB[0][0] + dB[0][1], steps - 1); });
    p.addTask([&]() {m7 = strassen<opType>(dA[0][1] - dA[1][1], dB[1][0] + dB[1][1], steps - 1); });
    p.completeTasksAndStop();
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
    if constexpr (avx::IsAvailable) return strassenWithStaticPadding<BaseMulType::Avx>(a, b, steps);
    else static_assert(avx::IsAvailable && always_false_v<M1>, avx_StaticAssertMessage);
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
    if constexpr (avx::IsAvailable) return strassenImpl<BaseMulType::Avx, Steps>(a, b);
    else static_assert(avx::IsAvailable && always_false_v<M1>, avx_StaticAssertMessage);
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
            operationEff<ArithmeticOperation::Assign>(rowSize, columnSize, columnSize, m, result[y][x], a + x * columnSize + y * rowSize*m);
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
    using namespace ArithmeticOperation;
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
    auto m1_a = operationWithAlloc<Add>(allocator, halfN, halfM, dA[0][0], dA[1][1]);
    auto m1_b = operationWithAlloc<Add>(allocator, halfM, halfQ, dB[0][0], dB[1][1]);
    lowLevelStrassen<opType>(m1_a, m1_b, halfN, halfM, halfQ, steps - 1, m1, allocator);
    allocator.dealloc(m1_b, halfM*halfQ);
    allocator.dealloc(m1_a, halfN*halfM);

    auto m2 = allocator.alloc(halfN * halfQ);
    auto m2_a = operationWithAlloc<Add>(allocator, halfN, halfM, dA[1][0], dA[1][1]);
    lowLevelStrassen<opType>(m2_a, dB[0][0], halfN, halfM, halfQ, steps - 1, m2, allocator);
    allocator.dealloc(m2_a, halfN*halfM);

    auto m3 = allocator.alloc(halfN * halfQ);
    auto m3_b = operationWithAlloc<Sub>(allocator, halfM, halfQ, dB[0][1], dB[1][1]);
    lowLevelStrassen<opType>(dA[0][0], m3_b, halfN, halfM, halfQ, steps - 1, m3, allocator);
    allocator.dealloc(m3_b, halfM*halfQ);

    auto m4 = allocator.alloc(halfN * halfQ);
    auto m4_b = operationWithAlloc<Sub>(allocator, halfM, halfQ, dB[1][0], dB[0][0]);
    lowLevelStrassen<opType>(dA[1][1], m4_b, halfN, halfM, halfQ, steps - 1, m4, allocator);
    allocator.dealloc(m4_b, halfM*halfQ);

    auto m5 = allocator.alloc(halfN * halfQ);
    auto m5_a = operationWithAlloc<Add>(allocator, halfN, halfM, dA[0][0], dA[0][1]);
    lowLevelStrassen<opType>(m5_a, dB[1][1], halfN, halfM, halfQ, steps - 1, m5, allocator);
    allocator.dealloc(m5_a, halfN*halfM);

    auto m6 = allocator.alloc(halfN * halfQ);
    auto m6_a = operationWithAlloc<Sub>(allocator, halfN, halfM, dA[1][0], dA[0][0]);
    auto m6_b = operationWithAlloc<Add>(allocator, halfM, halfQ, dB[0][0], dB[0][1]);
    lowLevelStrassen<opType>(m6_a, m6_b, halfN, halfM, halfQ, steps - 1, m6, allocator);
    allocator.dealloc(m6_b, halfM*halfQ);
    allocator.dealloc(m6_a, halfN*halfM);

    auto m7 = allocator.alloc(halfN * halfQ);
    auto m7_a = operationWithAlloc<Sub>(allocator, halfN, halfM, dA[0][1], dA[1][1]);
    auto m7_b = operationWithAlloc<Add>(allocator, halfM, halfQ, dB[1][0], dB[1][1]);
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
        operationEff<ArithmeticOperation::Assign>(n, m, paddedSizes[1], m, newA, a);
        operationEff<ArithmeticOperation::Assign>(m, q, paddedSizes[2], q, newB, b);
        lowLevelStrassen<opType>(newA, newB, paddedSizes[0], paddedSizes[1], paddedSizes[2], steps, newResult, allocator);
        operationEff<ArithmeticOperation::Assign>(n, q, q, paddedSizes[2], result, newResult);
        allocator.dealloc(newResult, paddedSizes[0] * paddedSizes[2]);
        allocator.dealloc(newB, paddedSizes[1] * paddedSizes[2]);
        allocator.dealloc(newA, paddedSizes[0] * paddedSizes[1]);
    } else {
        lowLevelStrassen<opType>(a, b, n, m, q, steps, result, allocator);
    }
}

template<BaseMulType opType, int Steps, int n, int m, int q, typename T>
void lowLevelStrassen(T* a, T* b, T* c, StackAllocator<T>& allocator, std::false_type) {
    mul<n, m, q, opType>(c, a, b);
}
template<BaseMulType opType, int Steps, int n, int m, int q, typename T>
void lowLevelStrassen(T* a, T* b, T* c, StackAllocator<T>& allocator, std::true_type unused = std::true_type{}) {
    using namespace ArithmeticOperation;
    if constexpr (Steps <= 0) {
        lowLevelStrassen<opType, Steps, n, m, q>(a, b, c, allocator, std::false_type{});
        return;
    }

    auto dA = strassenDivide<2, 2>(a, n, m, allocator);
    auto dB = strassenDivide<2, 2>(b, m, q, allocator);

    constexpr int hn = n / 2;
    constexpr int hm = m / 2;
    constexpr int hq = q / 2;

    auto m1 = allocator.alloc(hn * hq);
    auto m1_a = operationWithAlloc<hn, hm, Add>(allocator, dA[0][0], dA[1][1]);
    auto m1_b = operationWithAlloc<hm, hq, Add>(allocator, dB[0][0], dB[1][1]);
    lowLevelStrassen<opType, Steps - 1, hn, hm, hq>(m1_a, m1_b, m1, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    allocator.dealloc(m1_b, hm*hq);
    allocator.dealloc(m1_a, hn*hm);

    auto m2 = allocator.alloc(hn * hq);
    auto m2_a = operationWithAlloc<hn, hm, Add>(allocator, dA[1][0], dA[1][1]);
    lowLevelStrassen<opType, Steps - 1, hn, hm, hq>(m2_a, dB[0][0], m2, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    allocator.dealloc(m2_a, hn*hm);

    auto m3 = allocator.alloc(hn * hq);
    auto m3_b = operationWithAlloc<hm, hq, Sub>(allocator, dB[0][1], dB[1][1]);
    lowLevelStrassen<opType, Steps - 1, hn, hm, hq>(dA[0][0], m3_b, m3, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    allocator.dealloc(m3_b, hm*hq);

    auto m4 = allocator.alloc(hn * hq);
    auto m4_b = operationWithAlloc<hm, hq, Sub>(allocator, dB[1][0], dB[0][0]);
    lowLevelStrassen<opType, Steps - 1, hn, hm, hq>(dA[1][1], m4_b, m4, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    allocator.dealloc(m4_b, hm*hq);

    auto m5 = allocator.alloc(hn * hq);
    auto m5_a = operationWithAlloc<hn, hm, Add>(allocator, dA[0][0], dA[0][1]);
    lowLevelStrassen<opType, Steps - 1, hn, hm, hq>(m5_a, dB[1][1], m5, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    allocator.dealloc(m5_a, hn*hm);

    auto m6 = allocator.alloc(hn * hq);
    auto m6_a = operationWithAlloc<hn, hm, Sub>(allocator, dA[1][0], dA[0][0]);
    auto m6_b = operationWithAlloc<hm, hq, Add>(allocator, dB[0][0], dB[0][1]);
    lowLevelStrassen<opType, Steps - 1, hn, hm, hq>(m6_a, m6_b, m6, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    allocator.dealloc(m6_b, hm*hq);
    allocator.dealloc(m6_a, hn*hm);

    auto m7 = allocator.alloc(hn * hq);
    auto m7_a = operationWithAlloc<hn, hm, Sub>(allocator, dA[0][1], dA[1][1]);
    auto m7_b = operationWithAlloc<hm, hq, Add>(allocator, dB[1][0], dB[1][1]);
    lowLevelStrassen<opType, Steps - 1, hn, hm, hq>(m7_a, m7_b, m7, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
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
template<BaseMulType opType, int Steps, int n, int m, int q, typename T>
void lowLevelStrassen(T* result, T* a, T* b) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    for (int i = 0; i < Steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += 5*StackAllocator<T>::Allign(eN*eM) + 5*StackAllocator<T>::Allign(eM*eQ) + 7*StackAllocator<T>::Allign(eN*eQ);
    }
    StackAllocator<T> allocator(expected);
    lowLevelStrassen<opType, Steps, n, m, q>(a, b, result, allocator);
}

template<BaseMulType opType, typename M1, typename M2>
auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    auto result = a.createNew(a.rowCount(), b.columnCount());
    lowLevelStrassenWithStaticPadding<opType>(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), steps);
    return result;
}
template<BaseMulType opType, int Steps, typename M1, typename M2>
auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        lowLevelStrassen<opType, Steps, M1::CRow(), M1::CCol(), M2::CCol()>(result.data(), a.data(), b.data());
        return result;
    }
    else {
        auto result = a.createNew(a.rowCount(), b.columnCount());
        lowLevelStrassenWithStaticPadding<opType>(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), Steps);
        return result;
    }
}
template<typename M1, typename M2>
auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return lowLevelStrassen<BaseMulType::Naive>(a, b, steps);
}
template<typename M1, typename M2>
auto lowLevelAvxStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if constexpr (avx::IsAvailable) return lowLevelStrassen<BaseMulType::Avx>(a, b, steps);
    else static_assert(avx::IsAvailable && always_false_v<M1>, avx_StaticAssertMessage);
}
template<typename M1, typename M2>
auto lowLevelAutoStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    return lowLevelStrassen<AutomaticBaseMulType<typename M1::ValueType>>(a, b, steps);
}
template<int Steps, typename M1, typename M2>
auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return lowLevelStrassen<BaseMulType::Naive, Steps>(a, b);
}
template<int Steps, typename M1, typename M2>
auto lowLevelAvxStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (avx::IsAvailable) return lowLevelStrassen<BaseMulType::Avx, Steps>(a, b);
    else static_assert(avx::IsAvailable && always_false_v<M1>, avx_StaticAssertMessage);
}
template<int Steps, typename M1, typename M2>
auto lowLevelAutoStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return lowLevelStrassen<AutomaticBaseMulType<typename M1::ValueType>, Steps>(a, b);
}

template<BaseMulType opType, typename T>
void minSpaceStrassenMul(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, StackAllocator<T>& allocator) {
    if constexpr (opType == BaseMulType::NaiveEff) {
        naiveMul(c, a, b, n, m, q, effC, effA, effB);
        return;
    }
    if constexpr (opType == BaseMulType::AvxEff) {
        avxParallelMul(c, a, b, n, m, q, effC, effA, effB);
        return;
    }
    T* cx = c;
    T* ax = a;
    T* bx = b;
    if (effA != m) {
        ax = allocator.alloc(n*m);
        operationEff<ArithmeticOperation::Assign>(n, m, m, effA, ax, a);
    }
    if (effB != q) {
        bx = allocator.alloc(m*q);
        operationEff<ArithmeticOperation::Assign>(m, q, q, effB, bx, b);
    }
    if (effC != q) {
        cx = allocator.alloc(n*q);
    }

    mul<opType>(cx, ax, bx, n, m, q);

    if (effC != q) {
        operationEff<ArithmeticOperation::Assign>(n, q, effC, q, c, cx);
        allocator.dealloc(cx, n*q);
    }
    if (effB != q) allocator.dealloc(bx, m*q);
    if (effA != m) allocator.dealloc(ax, n*m);
}

template<BaseMulType opType, typename T>
void minSpaceStrassen(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
    using namespace ArithmeticOperation;
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

    operationEff<Add>(hn, hm, hm, effA, tempA, dA[0][0], dA[1][1]);
    operationEff<Add>(hm, hq, hq, effB, tempB, dB[0][0], dB[1][1]);
    minSpaceStrassen<opType>(dC[0][0], tempA, tempB, hn, hm, hq, effC, hm, hq, steps - 1, allocator);
    operationEff<Assign>(hn, hq, effC, effC, dC[1][1], dC[0][0]);

    operationEff<Add>(hn, hm, hm, effA, tempA, dA[1][0], dA[1][1]);
    minSpaceStrassen<opType>(dC[1][0], tempA, dB[0][0], hn, hm, hq, effC, hm, effB, steps - 1, allocator);
    operationEff<SubAssign>(hn, hq, effC, effC, dC[1][1], dC[1][0]);

    operationEff<Sub>(hm, hq, hq, effB, tempB, dB[0][1], dB[1][1]);
    minSpaceStrassen<opType>(dC[0][1], dA[0][0], tempB, hn, hm, hq, effC, effA, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, effC, dC[1][1], dC[0][1]);

    operationEff<Sub>(hm, hq, hq, effB, tempB, dB[1][0], dB[0][0]);
    minSpaceStrassen<opType>(tempC, dA[1][1], tempB, hn, hm, hq, hq, effA, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[0][0], tempC);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[1][0], tempC);

    operationEff<Add>(hn, hm, hm, effA, tempA, dA[0][0], dA[0][1]);
    minSpaceStrassen<opType>(tempC, tempA, dB[1][1], hn, hm, hq, hq, hm, effB, steps - 1, allocator);
    operationEff<SubAssign>(hn, hq, effC, hq, dC[0][0], tempC);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[0][1], tempC);

    operationEff<Sub>(hn, hm, hm, effA, tempA, dA[1][0], dA[0][0]);
    operationEff<Add>(hm, hq, hq, effB, tempB, dB[0][0], dB[0][1]);
    minSpaceStrassen<opType>(tempC, tempA, tempB, hn, hm, hq, hq, hm, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[1][1], tempC);

    operationEff<Sub>(hn, hm, hm, effA, tempA, dA[0][1], dA[1][1]);
    operationEff<Add>(hm, hq, hq, effB, tempB, dB[1][0], dB[1][1]);
    minSpaceStrassen<opType>(tempC, tempA, tempB, hn, hm, hq, hq, hm, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[0][0], tempC);

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
    if (steps >= 1 && opType != BaseMulType::NaiveEff && opType != BaseMulType::AvxEff) {
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
        operationEff<ArithmeticOperation::Assign>(n, m, paddedSizes[1], m, newA, a);
        operationEff<ArithmeticOperation::Assign>(m, q, paddedSizes[2], q, newB, b);
        minSpaceStrassen<opType>(newC, newA, newB, paddedSizes[0], paddedSizes[1], paddedSizes[2], paddedSizes[2], paddedSizes[1], paddedSizes[2], steps, allocator);
        operationEff<ArithmeticOperation::Assign>(n, q, q, paddedSizes[2], c, newC);
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
        operationEff<n, m, m, effA, ArithmeticOperation::Assign>(ax, a);
    }
    if (effB != q) {
        bx = allocator.alloc(m*q);
        operationEff<m, q, q, effB, ArithmeticOperation::Assign>(bx, b);
    }
    if (effC != q) {
        cx = allocator.alloc(n*q);
    }

    mul<opType>(cx, ax, bx, n, m, q);

    if (effC != q) {
        operationEff<n, q, effC, q, ArithmeticOperation::Assign>(c, cx);
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
    using namespace ArithmeticOperation;
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

    operationEff<hn, hm, hm, effA, Add>(tempA, dA[0][0], dA[1][1]);
    operationEff<hm, hq, hq, effB, Add>(tempB, dB[0][0], dB[1][1]);
    minSpaceStrassen<hn, hm, hq, effC, hm, hq, Steps - 1, opType>(dC[0][0], tempA, tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operationEff<hn, hq, effC, effC, Assign>(dC[1][1], dC[0][0]);

    operationEff<hn, hm, hm, effA, Add>(tempA, dA[1][0], dA[1][1]);
    minSpaceStrassen<hn, hm, hq, effC, hm, effB, Steps - 1, opType>(dC[1][0], tempA, dB[0][0], allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operationEff<hn, hq, effC, effC, SubAssign>(dC[1][1], dC[1][0]);

    operationEff<hm, hq, hq, effB, Sub>(tempB, dB[0][1], dB[1][1]);
    minSpaceStrassen<hn, hm, hq, effC, effA, hq, Steps - 1, opType>(dC[0][1], dA[0][0], tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operationEff<hn, hq, effC, effC, AddAssign>(dC[1][1], dC[0][1]);

    operationEff<hm, hq, hq, effB, Sub>(tempB, dB[1][0], dB[0][0]);
    minSpaceStrassen<hn, hm, hq, hq, effA, hq, Steps - 1, opType>(tempC, dA[1][1], tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operationEff<hn, hq, effC, hq, AddAssign>(dC[0][0], tempC);
    operationEff<hn, hq, effC, hq, AddAssign>(dC[1][0], tempC);

    operationEff<hn, hm, hm, effA, Add>(tempA, dA[0][0], dA[0][1]);
    minSpaceStrassen<hn, hm, hq, hq, hm, effB, Steps - 1, opType>(tempC, tempA, dB[1][1], allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operationEff<hn, hq, effC, hq, SubAssign>(dC[0][0], tempC);
    operationEff<hn, hq, effC, hq, AddAssign>(dC[0][1], tempC);

    operationEff<hn, hm, hm, effA, Sub>(tempA, dA[1][0], dA[0][0]);
    operationEff<hm, hq, hq, effB, Add>(tempB, dB[0][0], dB[0][1]);
    minSpaceStrassen< hn, hm, hq, hq, hm, hq, Steps - 1, opType>(tempC, tempA, tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operationEff<hn, hq, effC, hq, AddAssign>(dC[1][1], tempC);

    operationEff<hn, hm, hm, effA, Sub>(tempA, dA[0][1], dA[1][1]);
    operationEff<hm, hq, hq, effB, Add>(tempB, dB[1][0], dB[1][1]);
    minSpaceStrassen<hn, hm, hq, hq, hm, hq, Steps - 1, opType>(tempC, tempA, tempB, allocator, typename std::conditional<(Steps > 1), std::true_type, std::false_type>::type());
    operationEff<hn, hq, effC, hq, AddAssign>(dC[0][0], tempC);

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
    if constexpr (avx::IsAvailable) return minSpaceStrassen<BaseMulType::Avx>(a, b, steps);
    else static_assert(avx::IsAvailable && always_false_v<M1>, avx_StaticAssertMessage);
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
    if constexpr (avx::IsAvailable) return minSpaceStrassen<Steps, BaseMulType::Avx>(a, b);
    else static_assert(avx::IsAvailable && always_false_v<M1>, avx_StaticAssertMessage);
}
template<int Steps, typename M1, typename M2>
auto minSpaceAutoStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    return minSpaceStrassen<Steps, AutomaticBaseMulType<typename M1::ValueType>>(a, b);
}

#endif





template<typename T> void naiveMul(T* result, const T* a, const T* b, int n, int m, int q, int effR, int effA, int effB) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < q; ++j) {
            T sum{};
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * effA] * b[j + k * effB];
            }
            result[j + i * effR] = sum;
        }
    }
}
template<typename T> void peelingMulAdd(T* result, const T* a, const T* b, int n, int m, int q, int effR, int effA, int effB) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < q; ++j) {
            T& sum = result[j + i * effR];
            for (int k = 0; k < m; ++k) {
                sum += a[k + i * effA] * b[j + k * effB];
            }
        }
    }
}

template<BaseMulType opType, typename T>
void minSpaceStrassenWithDynamicPeeling2(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator);

template<BaseMulType opType, typename T>
void minSpaceStrassenWithDynamicPeeling(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
    using namespace ArithmeticOperation;
    auto dA = strassenDivideView<2, 2>(a, n, m, effA);
    auto dB = strassenDivideView<2, 2>(b, m, q, effB);
    auto dC = strassenDivideView<2, 2>(c, n, q, effC);
    
    int hn = n / 2;
    int hm = m / 2;
    int hq = q / 2;

    auto tempA = allocator.alloc(hn * hm);
    auto tempB = allocator.alloc(hm * hq);
    auto tempC = allocator.alloc(hn * hq);

    operationEff<Add>(hn, hm, hm, effA, tempA, dA[0][0], dA[1][1]);
    operationEff<Add>(hm, hq, hq, effB, tempB, dB[0][0], dB[1][1]);
    minSpaceStrassenWithDynamicPeeling2<opType>(dC[0][0], tempA, tempB, hn, hm, hq, effC, hm, hq, steps - 1, allocator);
    operationEff<Assign>(hn, hq, effC, effC, dC[1][1], dC[0][0]);

    operationEff<Add>(hn, hm, hm, effA, tempA, dA[1][0], dA[1][1]);
    minSpaceStrassenWithDynamicPeeling2<opType>(dC[1][0], tempA, dB[0][0], hn, hm, hq, effC, hm, effB, steps - 1, allocator);
    operationEff<SubAssign>(hn, hq, effC, effC, dC[1][1], dC[1][0]);

    operationEff<Sub>(hm, hq, hq, effB, tempB, dB[0][1], dB[1][1]);
    minSpaceStrassenWithDynamicPeeling2<opType>(dC[0][1], dA[0][0], tempB, hn, hm, hq, effC, effA, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, effC, dC[1][1], dC[0][1]);

    operationEff<Sub>(hm, hq, hq, effB, tempB, dB[1][0], dB[0][0]);
    minSpaceStrassenWithDynamicPeeling2<opType>(tempC, dA[1][1], tempB, hn, hm, hq, hq, effA, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[0][0], tempC);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[1][0], tempC);

    operationEff<Add>(hn, hm, hm, effA, tempA, dA[0][0], dA[0][1]);
    minSpaceStrassenWithDynamicPeeling2<opType>(tempC, tempA, dB[1][1], hn, hm, hq, hq, hm, effB, steps - 1, allocator);
    operationEff<SubAssign>(hn, hq, effC, hq, dC[0][0], tempC);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[0][1], tempC);

    operationEff<Sub>(hn, hm, hm, effA, tempA, dA[1][0], dA[0][0]);
    operationEff<Add>(hm, hq, hq, effB, tempB, dB[0][0], dB[0][1]);
    minSpaceStrassenWithDynamicPeeling2<opType>(tempC, tempA, tempB, hn, hm, hq, hq, hm, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[1][1], tempC);

    operationEff<Sub>(hn, hm, hm, effA, tempA, dA[0][1], dA[1][1]);
    operationEff<Add>(hm, hq, hq, effB, tempB, dB[1][0], dB[1][1]);
    minSpaceStrassenWithDynamicPeeling2<opType>(tempC, tempA, tempB, hn, hm, hq, hq, hm, hq, steps - 1, allocator);
    operationEff<AddAssign>(hn, hq, effC, hq, dC[0][0], tempC);

    allocator.dealloc(tempC, hn*hq);
    allocator.dealloc(tempB, hm*hq);
    allocator.dealloc(tempA, hn*hm);
}

template<BaseMulType opType, typename T>
void minSpaceStrassenWithDynamicPeeling2(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
    if (steps <= 0) {
        minSpaceStrassenMul<opType>(c, a, b, n, m, q, effC, effA, effB, allocator);
        return;
    }
    int nPeeled = n;
    int mPeeled = m;
    int qPeeled = q;
    if (n & 1) nPeeled -= 1;
    if (m & 1) mPeeled -= 1;
    if (q & 1) qPeeled -= 1;
    minSpaceStrassenWithDynamicPeeling<opType>(c, a, b, nPeeled, mPeeled, qPeeled, effC, effA, effB, steps, allocator);

    // C11 += A12 * B21
    peelingMulAdd(c, a+mPeeled, b+effB*mPeeled, nPeeled, m&1, qPeeled, effC, effA, effB);

    // C12 = A11 * B12 + A12 * B22
    naiveMul(c + qPeeled, a, b + qPeeled, nPeeled, m, q&1, effC, effA, effB);

    // C21 = A21 * B11 + A22 * B21
    naiveMul(c + effC*nPeeled, a+effA*nPeeled, b, n&1, m, qPeeled, effC, effA, effB);

    // C22 = A21 * B12 + A22 * B22
    naiveMul(c+ effC*nPeeled+qPeeled, a+effA*nPeeled, b+qPeeled, n&1, m, q&1, effC, effA, effB);
}

template<BaseMulType opType, typename T> void minSpaceStrassenWithDynamicPeeling(T* c, T* a, T* b, int n, int m, int q, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    
    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += StackAllocator<T>::Allign(eN*eM) + StackAllocator<T>::Allign(eM*eQ) + StackAllocator<T>::Allign(eN*eQ);
    }
    if (steps >= 1 && opType != BaseMulType::NaiveEff && opType != BaseMulType::AvxEff) {
        expected += std::max(StackAllocator<T>::Allign(eN*eM), StackAllocator<T>::Allign(eM*eQ)) + StackAllocator<T>::Allign(eN*eQ);
    }
    StackAllocator<T> allocator(expected);

    minSpaceStrassenWithDynamicPeeling2<opType>(c, a, b, n, m, q, q, m, q, steps, allocator);
}

template<BaseMulType opType, typename M1, typename M2>
auto minSpaceStrassenWithDynamicPeeling(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    auto result = a.createNew(a.rowCount(), b.columnCount());
    minSpaceStrassenWithDynamicPeeling<opType>(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), steps);
    return result;
}














template<BaseMulType opType, typename T> void LowLevelParallelStrassenMinSpace(T* c, T* a, T* b, int n, int m, int q, int effA, int effB, int steps) {
    int expected = 0;
    int eN = n;
    int eM = m;
    int eQ = q;
    
    for (int i = 0; i < steps; ++i) {
        eN /= 2;
        eM /= 2;
        eQ /= 2;
        expected += StackAllocator<T>::Allign(eN*eM) + StackAllocator<T>::Allign(eM*eQ) + StackAllocator<T>::Allign(eN*eQ);
    }
    if (opType != BaseMulType::NaiveEff && opType != BaseMulType::AvxEff) {
        expected += std::max(StackAllocator<T>::Allign(eN*eM), StackAllocator<T>::Allign(eM*eQ)) + StackAllocator<T>::Allign(eN*eQ);
    }
    StackAllocator<T> allocator(expected);

    minSpaceStrassenWithDynamicPeeling2<opType>(c, a, b, n, m, q, q, effA, effB, steps, allocator);
}

template<BaseMulType opType, typename T>
void LowLevelParallelStrassen(T* c, T* a, T* b, int n, int m, int q, int effC, int effA, int effB, int steps) {
    using namespace ArithmeticOperation;
    
    int hn = n / 2;
    int hm = m / 2;
    int hq = q / 2;

    StackAllocator<T> allocator(5 * StackAllocator<T>::Allign(hn*hm) + 5 * StackAllocator<T>::Allign(hm*hq) + 7 * StackAllocator<T>::Allign(hn*hq));

    auto dA = strassenDivideView<2, 2>(a, n, m, effA);
    auto dB = strassenDivideView<2, 2>(b, m, q, effB);

    auto m1 = allocator.alloc(hn * hq);
    auto m2 = allocator.alloc(hn * hq);
    auto m3 = allocator.alloc(hn * hq);
    auto m4 = allocator.alloc(hn * hq);
    auto m5 = allocator.alloc(hn * hq);
    auto m6 = allocator.alloc(hn * hq);
    auto m7 = allocator.alloc(hn * hq);
    auto m1_a = allocator.alloc(hn*hm);
    auto m1_b = allocator.alloc(hm*hq);
    auto m2_a = allocator.alloc(hn*hm);
    auto m3_b = allocator.alloc(hm*hq);
    auto m4_b = allocator.alloc(hm*hq);
    auto m5_a = allocator.alloc(hn*hm);
    auto m6_a = allocator.alloc(hn*hm);
    auto m6_b = allocator.alloc(hm*hq);
    auto m7_a = allocator.alloc(hn*hm);
    auto m7_b = allocator.alloc(hm*hq);
    
    ThreadPool pool;

    pool.addTask([=]() {
        operationEff<Add>(hn, hm, hm, effA, m1_a, dA[0][0], dA[1][1]);
        operationEff<Add>(hm, hq, hq, effB, m1_b, dB[0][0], dB[1][1]);
        LowLevelParallelStrassenMinSpace<opType>(m1, m1_a, m1_b, hn, hm, hq, hm, hq, steps - 1);
    });
    pool.addTask([=]() {
        operationEff<Add>(hn, hm, hm, effA, m2_a, dA[1][0], dA[1][1]);
        LowLevelParallelStrassenMinSpace<opType>(m2, m2_a, dB[0][0], hn, hm, hq, hm, effB, steps - 1);
    });
    pool.addTask([=]() {
        operationEff<Sub>(hm, hq, hq, effB, m3_b, dB[0][1], dB[1][1]);
        LowLevelParallelStrassenMinSpace<opType>(m3, dA[0][0], m3_b, hn, hm, hq, effA, hq, steps - 1);
    });
    pool.addTask([=]() {
        operationEff<Sub>(hm, hq, hq, effB, m4_b, dB[1][0], dB[0][0]);
        LowLevelParallelStrassenMinSpace<opType>(m4, dA[1][1], m4_b, hn, hm, hq, effA, hq, steps - 1);
    });
    pool.addTask([=]() {
        operationEff<Add>(hn, hm, hm, effA, m5_a, dA[0][0], dA[0][1]);
        LowLevelParallelStrassenMinSpace<opType>(m5, m5_a, dB[1][1], hn, hm, hq, hm, effB, steps - 1);
    });
    pool.addTask([=]() {
        operationEff<Sub>(hn, hm, hm, effA, m6_a, dA[1][0], dA[0][0]);
        operationEff<Add>(hm, hq, hq, effB, m6_b, dB[0][0], dB[0][1]);
        LowLevelParallelStrassenMinSpace<opType>(m6, m6_a, m6_b, hn, hm, hq, hm, hq, steps - 1);
    });
    pool.addTask([=]() {
        operationEff<Sub>(hn, hm, hm, effA, m7_a, dA[0][1], dA[1][1]);
        operationEff<Add>(hm, hq, hq, effB, m7_b, dB[1][0], dB[1][1]);
        LowLevelParallelStrassenMinSpace<opType>(m7, m7_a, m7_b, hn, hm, hq, hm, hq, steps - 1);
    });

    pool.completeTasksAndStop();

    auto dC = strassenDivideView<2, 2>(c, n, q, effC);
    for (int i = 0; i < hn; ++i) {
        for (int j = 0; j < hq; ++j) {
            int dstIndex = j + i * effC;
            int index = j + i * hq;
            dC[0][0][dstIndex] = m1[index] + m4[index] - m5[index] + m7[index];
            dC[0][1][dstIndex] = m3[index] + m5[index];
            dC[1][0][dstIndex] = m2[index] + m4[index];
            dC[1][1][dstIndex] = m1[index] + m3[index] - m2[index] + m6[index];
        }
    }

    allocator.dealloc(m7_b, hm*hq);
    allocator.dealloc(m7_a, hn*hm);
    allocator.dealloc(m6_b, hm*hq);
    allocator.dealloc(m6_a, hn*hm);
    allocator.dealloc(m5_a, hn*hm);
    allocator.dealloc(m4_b, hm*hq);
    allocator.dealloc(m3_b, hm*hq);
    allocator.dealloc(m2_a, hn*hm);
    allocator.dealloc(m1_b, hm*hq);
    allocator.dealloc(m1_a, hn*hm);
    allocator.dealloc(m7, hn*hq);
    allocator.dealloc(m6, hn*hq);
    allocator.dealloc(m5, hn*hq);
    allocator.dealloc(m4, hn*hq);
    allocator.dealloc(m3, hn*hq);
    allocator.dealloc(m2, hn*hq);
    allocator.dealloc(m1, hn*hq);
}

template<BaseMulType opType, typename T>
void LowLevelParallelStrassenX(T* c, T* a, T* b, int n, int m, int q, int steps) {
    int nPeeled = n;
    int mPeeled = m;
    int qPeeled = q;
    if (n & 1) nPeeled -= 1;
    if (m & 1) mPeeled -= 1;
    if (q & 1) qPeeled -= 1;

    LowLevelParallelStrassen<opType>(c, a, b, nPeeled, mPeeled, qPeeled, q, m, q, steps);

    peelingMulAdd(c, a + mPeeled, b + q * mPeeled, nPeeled, m & 1, qPeeled, q, m, q);
    naiveMul(c + qPeeled, a, b + qPeeled, nPeeled, m, q & 1, q, m, q);
    naiveMul(c + q * nPeeled, a + m * nPeeled, b, n & 1, m, qPeeled, q, m, q);
    naiveMul(c + q * nPeeled + qPeeled, a + m * nPeeled, b + qPeeled, n & 1, m, q & 1, q, m, q);
}



template<BaseMulType opType, typename M1, typename M2>
auto LowLevelParallelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
    if (steps <= 0) return mul<opType>(a, b);
    auto result = a.createNew(a.rowCount(), b.columnCount());
    LowLevelParallelStrassenX<opType>(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), steps);
    return result;
}






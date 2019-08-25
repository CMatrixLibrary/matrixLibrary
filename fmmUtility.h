#ifndef FMM_UTILITY_H
#define FMM_UTILITY_H

#include "MatrixInterface.h"
#include "Matrix.h"
#include "StackAllocator.h"
#include "genericArithmeticOperations.h"
#include "avxMul.h"
#include "blasMul.h"
#include "naiveMul.h"
#include "parallelMul.h"
#include "blockMul.h"
#include "parallelBlockMul.h"
#include <algorithm>
#include <optional>
#include <tuple>
#include <cmath>

namespace fmm {

    // enums used to encode running method.
    // For example MinSpace algorithm with Naive base multiplication, static padding
    // and normal base multiplication size scheme would be encoded as 0x01020104
    enum Algorithm {
        AlgorithmAutomatic      = 0x00,
        HighLevel               = 0x01,
        LowLevel                = 0x02,
        MinSpace                = 0x04,
        HighLevelParallel       = 0x08,
        LowLevelParallel        = 0x10,
        AlgorithmMASK           = 0xff
    };
    enum BaseMulType {
        BaseMulTypeAutomatic    = 0x00'00,
        Naive                   = 0x01'00,
        Parallel                = 0x02'00,
        Block                   = 0x04'00,
        ParallelBlock           = 0x08'00,
        Avx                     = 0x10'00,
        ParallelAvx             = 0x20'00,
        Blas                    = 0x40'00,
        BaseMulTypeMASK         = 0xff'00
    };
    enum ResizeStrategy {
        DynamicPeeling          = 0x01'00'00,
        StaticPadding           = 0x02'00'00,
        ResizeStrategyMASK      = 0xff'00'00,
        ResizeStrategyAutomatic = DynamicPeeling
    };
    enum BaseMulSize {
        Normal                  = 0x01'00'00'00,
        Effective               = 0x02'00'00'00,
        BaseMulSizeMASK         = 0xff'00'00'00,
        BaseMulSizeAutomatic    = Effective
    };
    
    template<int Method, Algorithm alg> 
    constexpr int getNewWithAlgorithm = (Method & (BaseMulTypeMASK | ResizeStrategyMASK | BaseMulSizeMASK)) | alg;

    template<int Method, int Flag> constexpr bool contains = (Method & Flag) == Flag;

}

// utility for unknown sizes at compile-time
namespace fmm::detail {

    // operations using calculate<>
    template<ArithmeticOperation::OpType... ops, typename T, typename... Ts> T* operationEffWithAlloc(StackAllocator<T>& allocator, int n, int m, int effM, const T* arg, Ts&&... args) {
        auto result = allocator.alloc(n*m);
        operationEff<ops...>(n, m, m, effM, result, arg, args...);
        return result;
    }

    // divide
    template<int Rows, int Columns, typename T>
    std::array<std::array<T*, Columns>, Rows> divide(T* a, int n, int m, int effM, StackAllocator<T>& allocator) {
        std::array<std::array<T*, Columns>, Rows> result;
        auto rowSize = n / Rows;
        auto columnSize = m / Columns;
        for (int y = 0; y < Rows; ++y) {
            for (int x = 0; x < Columns; ++x) {
                result[y][x] = allocator.alloc(rowSize*columnSize);
                operationEff<ArithmeticOperation::Assign>(rowSize, columnSize, columnSize, effM, result[y][x], a + x * columnSize + y * rowSize*effM);
            }
        }
        return result;
    }
    template<int Rows, int Columns, typename T> auto divide(T* a, int n, int m, StackAllocator<T>& allocator) {
        return divide<Rows, Columns>(a, n, m, m, allocator);
    }

    template<int Rows, int Columns, typename T>
    std::array<std::array<T*, Columns>, Rows> divideView(T* a, int n, int m, int effM) {
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
    template<int Rows, int Columns, typename T> auto divideView(T* a, int n, int m) {
        return divideView<Rows, Columns>(a, n, m, m);
    }

    // static padding
    int staticPaddingNewSize(int size, int base, int steps) {
        if (steps <= 0) return size;
        int divisor = pow(base, steps);
        auto remainder = size % divisor;
        if (remainder == 0) return size;
        else                return size + divisor - remainder;
    }

    std::optional<std::array<int, 3>> staticPaddingNewSizes(int n, int m, int p, int baseN, int baseM, int baseP, int steps) {
        std::array<int, 3> paddedSizes = {
            staticPaddingNewSize(n, baseN, steps),
            staticPaddingNewSize(m, baseM, steps),
            staticPaddingNewSize(p, baseP, steps)
        };
        if (paddedSizes[0] > n || paddedSizes[1] > m || paddedSizes[2] > p) {
            return paddedSizes;
        } else {
            return std::nullopt;
        }
    }
    template<typename T> std::tuple<T*, T*, T*> allocateAndCopy(
        T* c, T* a, T* b, int n, int m, int p, int effA, int effB,
        int newN, int newM, int newP, StackAllocator<T>& allocator)
    {
        T* newA = allocator.alloc(newN*newM);
        T* newB = allocator.alloc(newM*newP);
        T* newC = allocator.alloc(newN*newP);
        for (int i = 0; i < newN * newM; ++i) {
            newA[i] = T{};
        }
        for (int i = 0; i < newM * newP; ++i) {
            newB[i] = T{};
        }
        operationEff<ArithmeticOperation::Assign>(n, m, newM, effA, newA, a);
        operationEff<ArithmeticOperation::Assign>(m, p, newP, effB, newB, b);
        return std::tuple(newC, newA, newB);
    }
    template<typename T> void copyAndDeallocate(
        T* newC, T* newA, T* newB, int n, int p,
        int newN, int newM, int newP, T* c, int effC, StackAllocator<T>& allocator
    ) {
        operationEff<ArithmeticOperation::Assign>(n, p, effC, newP, c, newC);
        allocator.dealloc(newC, newN*newP);
        allocator.dealloc(newB, newM*newP);
        allocator.dealloc(newA, newN*newM);
    }

    // dynamic peeling
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

    template<int BaseN, int BaseM, int BaseP> 
    std::tuple<int, int, int> peelingSizes(int n, int m, int p) {
        return { n - n % BaseN, m - m % BaseM, p - p % BaseP };
    }

    template<int BaseN, int BaseM, int BaseP, typename T> 
    void peelingMul(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB) {
        auto[nPeeled, mPeeled, pPeeled] = peelingSizes<BaseN, BaseM, BaseP>(n, m, p);
        peelingMulAdd(c, a + mPeeled, b + effB * mPeeled, nPeeled, m % BaseM, pPeeled, effC, effA, effB);
        naiveMul(c + pPeeled, a, b + pPeeled, nPeeled, m, p % BaseP, effC, effA, effB);
        naiveMul(c + effC * nPeeled, a + effA * nPeeled, b, n % BaseN, m, pPeeled, effC, effA, effB);
        naiveMul(c + effC * nPeeled + pPeeled, a + effA * nPeeled, b + pPeeled, n % BaseN, m, p % BaseP, effC, effA, effB);
    }

    // calculating space needed
    template<int BaseN, int BaseM, int BaseP, typename T> 
    int additionalRunningSpaceSize(int steps, int n, int m, int p, int nm_mul, int mp_mul, int np_mul) {
        int sum = 0;
        for (int i = 0; i < steps; ++i) {
            n /= BaseN;
            m /= BaseM;
            p /= BaseP;
            sum += nm_mul * StackAllocator<T>::Allign(n*m)
                + mp_mul * StackAllocator<T>::Allign(m*p)
                + np_mul * StackAllocator<T>::Allign(n*p);
        }
        return sum;
    }

    template<int BaseN> int sizeAfterLastStep(int n, int steps) {
        return n / pow(BaseN, steps);
    }

    template<typename T> int additionalStaticPaddingSize(int nPaddded, int mPadded, int pPadded) {
        return StackAllocator<T>::Allign(nPaddded * mPadded)
            + StackAllocator<T>::Allign(mPadded * pPadded)
            + StackAllocator<T>::Allign(nPaddded * pPadded);
    }
    
    // base multiplication
    template<int Method, typename T> 
    void baseMul(T* c, const T* a, const T* b, int n, int m, int p) {
        if constexpr (contains<Method, BaseMulType::Naive>)         naiveMul(c, a, b, n, m, p);
        if constexpr (contains<Method, BaseMulType::Parallel>)      parallelMul(c, a, b, n, m, p);
        if constexpr (contains<Method, BaseMulType::Block>)         blockMul(c, a, b, n, m, p);
        if constexpr (contains<Method, BaseMulType::ParallelBlock>) parallelBlockMul(c, a, b, n, m, p);
        if constexpr (contains<Method, BaseMulType::Avx>)           avx::mul(c, a, b, n, m, p);
        if constexpr (contains<Method, BaseMulType::ParallelAvx>)   avx::parallelMul(c, a, b, n, m, p);
        if constexpr (contains<Method, BaseMulType::Blas>)          blas::mul(c, a, b, n, m, p);
    }

    template<int Method, typename T> 
    void baseMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        if constexpr (contains<Method, BaseMulType::Naive>)         naiveMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (contains<Method, BaseMulType::Parallel>)      parallelMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (contains<Method, BaseMulType::Block>)         blockMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (contains<Method, BaseMulType::ParallelBlock>) parallelBlockMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (contains<Method, BaseMulType::Avx>)           avx::mul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (contains<Method, BaseMulType::ParallelAvx>)   avx::parallelMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (contains<Method, BaseMulType::Blas>)          blas::mul(c, a, b, n, m, p, effC, effA, effB);
    }

    template<int Method, typename T>
    void baseMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB, StackAllocator<T>& allocator) {
        if (effC != p || effA != m || effB != p) {
            if constexpr (contains<Method, BaseMulSize::Effective>) {
                baseMul<Method>(c, a, b, n, m, p, effC, effA, effB);
            } else if constexpr (contains<Method, BaseMulSize::Normal>) {
                T* cx = c;
                T* ax = const_cast<T*>(a);
                T* bx = const_cast<T*>(b);
                if (effA != m) {
                    ax = allocator.alloc(n*m);
                    operationEff<ArithmeticOperation::Assign>(n, m, m, effA, ax, a);
                }
                if (effB != p) {
                    bx = allocator.alloc(m*p);
                    operationEff<ArithmeticOperation::Assign>(m, p, p, effB, bx, b);
                }
                if (effC != p) cx = allocator.alloc(n*p);

                baseMul<Method>(cx, ax, bx, n, m, p);

                if (effC != p) {
                    operationEff<ArithmeticOperation::Assign>(n, p, effC, p, c, cx);
                    allocator.dealloc(cx, n*p);
                }
                if (effB != p) allocator.dealloc(bx, m*p);
                if (effA != m) allocator.dealloc(ax, n*m);
            }
        } else {
            baseMul<Method>(c, a, b, n, m, p);
        }
    }

    // additional utility
    template<int BaseN, int BaseM, int BaseP> std::tuple<int, int, int> divideSizes(int n, int m, int p) {
        return { n / BaseN, m / BaseM, p / BaseP };
    }


    template<int Method, int BaseN, int BaseM, int BaseP, typename FunctionImpl, typename T>
    void nextStep(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
        if (steps <= 0) {
            baseMul<Method>(c, a, b, n, m, p, effC, effA, effB, allocator);
            return;
        }
        if constexpr (contains<Method, ResizeStrategy::DynamicPeeling>) {
            auto[nPeeled, mPeeled, pPeeled] = peelingSizes<BaseN, BaseM, BaseP>(n, m, p);
            FunctionImpl::template Run<Method>(c, a, b, nPeeled, mPeeled, pPeeled, effC, effA, effB, steps, allocator);
            peelingMul<BaseN, BaseM, BaseP>(c, a, b, n, m, p, effC, effA, effB);
        } else {
            FunctionImpl::template Run<Method>(c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }

    
    // Min-Space
    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, typename FunctionImpl, typename T>
    void minSpaceRun(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps) {
        auto paddingSizesOpt = staticPaddingNewSizes(n, m, p, BaseN, BaseM, BaseP, steps);
        int runningSpace = 0;
        int nPadded = n;
        int mPadded = m;
        int pPadded = p;
        if (contains<Method, ResizeStrategy::StaticPadding> && paddingSizesOpt) {
            auto& paddingSizes = *paddingSizesOpt;
            nPadded = paddingSizes[0];
            mPadded = paddingSizes[1];
            pPadded = paddingSizes[2];
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(steps, nPadded, mPadded, pPadded, 1, 1, 1);

        if (contains<Method, BaseMulSize::Normal> && !(steps == 0 && effC==p && effA==m && effB==p)) {
            int nn = sizeAfterLastStep<BaseN>(nPadded, steps);
            int mm = sizeAfterLastStep<BaseM>(mPadded, steps);
            int pp = sizeAfterLastStep<BaseP>(pPadded, steps);
            if (steps == 0 && effC == p) {
                runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp));
            } else {
                runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp)) + StackAllocator<T>::Allign(nn*pp);
            }
        }

        StackAllocator<T> allocator(runningSpace);
        if (contains<Method, ResizeStrategy::StaticPadding> && paddingSizesOpt) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP, FunctionImpl>(newC, newA, newB, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP, FunctionImpl>(c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }


    // Low-level
    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, typename FunctionImpl, typename T>
    void lowLevelRun(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps) {
        auto paddingSizesOpt = staticPaddingNewSizes(n, m, p, BaseN, BaseM, BaseP, steps);
        int runningSpace = 0;
        int nPadded = n;
        int mPadded = m;
        int pPadded = p;
        if (contains<Method, ResizeStrategy::StaticPadding> && paddingSizesOpt) {
            auto& paddingSizes = *paddingSizesOpt;
            nPadded = paddingSizes[0];
            mPadded = paddingSizes[1];
            pPadded = paddingSizes[2];
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(steps, nPadded, mPadded, pPadded, 1, 1, MulCount);

        if (contains<Method, BaseMulSize::Normal> && steps == 0 && (effC!=p || effA!=m || effB!=p)) {
            int nn = sizeAfterLastStep<BaseN>(nPadded, steps);
            int mm = sizeAfterLastStep<BaseM>(mPadded, steps);
            int pp = sizeAfterLastStep<BaseP>(pPadded, steps);
            runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp)) + StackAllocator<T>::Allign(nn*pp);
        }

        StackAllocator<T> allocator(runningSpace);
        if (contains<Method, ResizeStrategy::StaticPadding> && paddingSizesOpt) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP, FunctionImpl>(newC, newA, newB, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP, FunctionImpl>(c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }


    // Parallel Low-Level
    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, typename FunctionImpl, typename T>
    void lowLevelParallelRun(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps) {
        auto paddingSizesOpt = staticPaddingNewSizes(n, m, p, BaseN, BaseM, BaseP, steps);
        int runningSpace = 0;
        int nPadded = n;
        int mPadded = m;
        int pPadded = p;
        if (contains<Method, ResizeStrategy::StaticPadding> && paddingSizesOpt) {
            auto& paddingSizes = *paddingSizesOpt;
            nPadded = paddingSizes[0];
            mPadded = paddingSizes[1];
            pPadded = paddingSizes[2];
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(std::min(1, steps), nPadded, mPadded, pPadded, 1 + BaseN*BaseM, 1 + BaseM*BaseP, MulCount);

        StackAllocator<T> allocator(runningSpace);
        if (contains<Method, ResizeStrategy::StaticPadding> && paddingSizesOpt) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP, FunctionImpl>(newC, newA, newB, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP, FunctionImpl>(c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }

    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, typename FunctionImpl, typename M1, typename M2>
    auto runAlgorithm(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        auto c = a.createNew(a.rowCount(), b.columnCount());

        if constexpr (contains<Method, Algorithm::LowLevel>) {
            lowLevelRun<Method, BaseN, BaseM, BaseP, MulCount, FunctionImpl>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        else if constexpr (contains<Method, Algorithm::MinSpace>) {
            minSpaceRun<Method, BaseN, BaseM, BaseP, MulCount, FunctionImpl>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        else if constexpr (contains<Method, Algorithm::LowLevelParallel>) {
            lowLevelParallelRun<Method, BaseN, BaseM, BaseP, MulCount, FunctionImpl>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        else { // Automatic
            // Here the algorithm should be chosen based on type and size of matrices
            lowLevelParallelRun<Method, BaseN, BaseM, BaseP, MulCount, FunctionImpl>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        
        return c;
    }
}

// utility for sizes known at compile time
namespace fmm::detail {

    // operations using calculate<>
    template<int n, int m, int effM, ArithmeticOperation::OpType... ops, typename T, typename... Ts> T* operationEffWithAlloc(StackAllocator<T>& allocator, const T* arg, Ts&&... args) {
        auto result = allocator.alloc(n*m);
        operationEff<n, m, m, effM, ops...>(result, arg, args...);
        return result;
    }

    // static padding
    template<int Base, int Exp> struct ConstexprPow {
        static constexpr int Run() {
            return Base * ConstexprPow<Base, Exp - 1>::Run();
        }
    };
    template<int Base> struct ConstexprPow<Base, 0> {
        static constexpr int Run() {
            return 1;
        }
    };
    template<int Size, int Base, int Steps> constexpr int staticPaddingNewSize() {
        if constexpr (Steps <= 0) return Size;
        else return ((Size % ConstexprPow<Base, Steps>::Run()) == 0) ?
            Size
            :
            Size + ConstexprPow<Base, Steps>::Run() - (Size % ConstexprPow<Base, Steps>::Run());
    }

    // base multiplication
    template<int Method, int n, int m, int p, int effC, int effA, int effB, typename T>
    void baseMul(T* c, const T* a, const T* b) {
        if constexpr (contains<Method, BaseMulType::Naive>)         naiveMul<n, m, p, effC, effA, effB>(c, a, b);
        if constexpr (contains<Method, BaseMulType::Parallel>)      parallelMul<n, m, p, effC, effA, effB>(c, a, b);
        if constexpr (contains<Method, BaseMulType::Block>)         blockMul<n, m, p, effC, effA, effB>(c, a, b);
        if constexpr (contains<Method, BaseMulType::ParallelBlock>) parallelBlockMul<n, m, p, effC, effA, effB>(c, a, b);
        if constexpr (contains<Method, BaseMulType::Avx>)           avx::mul<n, m, p, effC, effA, effB>(c, a, b);
        if constexpr (contains<Method, BaseMulType::ParallelAvx>)   avx::parallelMul<n, m, p, effC, effA, effB>(c, a, b);
        if constexpr (contains<Method, BaseMulType::Blas>)          blas::mul<n, m, p, effC, effA, effB>(c, a, b);
    }

    template<int Method, int n, int m, int p, int effC, int effA, int effB, typename T>
    void baseMul(T* c, const T* a, const T* b, StackAllocator<T>& allocator) {
        if constexpr (effC != p || effA != m || effB != p) {
            if constexpr (contains<Method, BaseMulSize::Effective>) {
                baseMul<Method, n, m, p, effC, effA, effB>(c, a, b);
            } else if constexpr (contains<Method, BaseMulSize::Normal>) {
                T* cx = c;
                T* ax = const_cast<T*>(a);
                T* bx = const_cast<T*>(b);
                if (effA != m) {
                    ax = allocator.alloc(n*m);
                    operationEff<ArithmeticOperation::Assign>(n, m, m, effA, ax, a);
                }
                if (effB != p) {
                    bx = allocator.alloc(m*p);
                    operationEff<ArithmeticOperation::Assign>(m, p, p, effB, bx, b);
                }
                if (effC != p) cx = allocator.alloc(n*p);

                baseMul<Method, n, m, p, p, m, p>(cx, ax, bx);

                if (effC != p) {
                    operationEff<ArithmeticOperation::Assign>(n, p, effC, p, c, cx);
                    allocator.dealloc(cx, n*p);
                }
                if (effB != p) allocator.dealloc(bx, m*p);
                if (effA != m) allocator.dealloc(ax, n*m);
            }
        } else {
            baseMul<Method, n, m, p, p, m, p>(c, a, b);
        }
    }


    template<int Method, int BaseN, int BaseM, int BaseP, int n, int m, int p, int effC, int effA, int effB, int steps, typename FunctionImpl, typename T> 
    struct NextStepImpl {
        void static Run(T* c, T* a, T* b, StackAllocator<T>& allocator) {
            if constexpr (contains<Method, ResizeStrategy::DynamicPeeling>) {
                constexpr auto nPeeled = n - n % BaseN;
                constexpr auto mPeeled = m - m % BaseM;
                constexpr auto pPeeled = p - p % BaseP;
                FunctionImpl::template Run<Method, nPeeled, mPeeled, pPeeled, effC, effA, effB, steps>(c, a, b, allocator);
                peelingMul<BaseN, BaseM, BaseP>(c, a, b, n, m, p, effC, effA, effB);
            } else {
                FunctionImpl::template Run<Method, n, m, p, effC, effA, effB, steps>(c, a, b, allocator);
            }
        }
    };
    template<int Method, int BaseN, int BaseM, int BaseP, int n, int m, int p, int effC, int effA, int effB, typename FunctionImpl, typename T> 
    struct NextStepImpl<Method, BaseN, BaseM, BaseP, n, m, p, effC, effA, effB, 0, FunctionImpl, T> {
        void static Run(T* c, T* a, T* b, StackAllocator<T>& allocator) {
            baseMul<Method, n, m, p, effC, effA, effB>(c, a, b, allocator);
        }
    };
    template<int Method, int BaseN, int BaseM, int BaseP, int n, int m, int p, int effC, int effA, int effB, int steps, typename FunctionImpl, typename T>
    void nextStep(T* c, T* a, T* b, StackAllocator<T>& allocator) {
        NextStepImpl< Method, BaseN, BaseM, BaseP, n, m, p, effC, effA, effB, steps, FunctionImpl, T>::Run(c, a, b, allocator);
    }

    
    // Min-Space
    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, int n, int m, int p, int effC, int effA, int effB, int steps, typename FunctionImpl, typename T>
    void minSpaceRun(T* c, T* a, T* b) {
        constexpr int nPadded = staticPaddingNewSize<n, BaseN, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        constexpr int mPadded = staticPaddingNewSize<m, BaseM, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        constexpr int pPadded = staticPaddingNewSize<p, BaseP, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        int runningSpace = 0;
        if (n != nPadded || m != mPadded || p != pPadded) {
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(steps, nPadded, mPadded, pPadded, 1, 1, 1);

        if (contains<Method, BaseMulSize::Normal> && !(steps == 0 && effC==p && effA==m && effB==p)) {
            int nn = sizeAfterLastStep<BaseN>(nPadded, steps);
            int mm = sizeAfterLastStep<BaseM>(mPadded, steps);
            int pp = sizeAfterLastStep<BaseP>(pPadded, steps);
            if (steps == 0 && effC == p) {
                runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp));
            } else {
                runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp)) + StackAllocator<T>::Allign(nn*pp);
            }
        }

        StackAllocator<T> allocator(runningSpace);
        if (n != nPadded || m != mPadded || p != pPadded) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, FunctionImpl>(newC, newA, newB, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP, n, m, p, effC, effA, effB, steps, FunctionImpl>(c, a, b, allocator);
        }
    }

    // Low-Level
    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, int n, int m, int p, int effC, int effA, int effB, int steps, typename FunctionImpl, typename T>
    void lowLevelRun(T* c, T* a, T* b) {
        constexpr int nPadded = staticPaddingNewSize<n, BaseN, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        constexpr int mPadded = staticPaddingNewSize<m, BaseM, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        constexpr int pPadded = staticPaddingNewSize<p, BaseP, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        int runningSpace = 0;
        if (n != nPadded || m != mPadded || p != pPadded) {
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(steps, nPadded, mPadded, pPadded, 1, 1, MulCount);

        if (contains<Method, BaseMulSize::Normal> && steps == 0 && (effC != p || effA != m || effB != p)) {
            int nn = sizeAfterLastStep<BaseN>(nPadded, steps);
            int mm = sizeAfterLastStep<BaseM>(mPadded, steps);
            int pp = sizeAfterLastStep<BaseP>(pPadded, steps);
            if (steps == 0 && effC == p) {
                runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp));
            } else {
                runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp)) + StackAllocator<T>::Allign(nn*pp);
            }
        }

        StackAllocator<T> allocator(runningSpace);
        if (n != nPadded || m != mPadded || p != pPadded) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, FunctionImpl>(newC, newA, newB, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP, n, m, p, effC, effA, effB, steps, FunctionImpl>(c, a, b, allocator);
        }
    }

    // Parallel Low-Level
    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, int n, int m, int p, int effC, int effA, int effB, int steps, typename FunctionImpl, typename T>
    void lowLevelParallelRun(T* c, T* a, T* b) {
        constexpr int nPadded = staticPaddingNewSize<n, BaseN, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        constexpr int mPadded = staticPaddingNewSize<m, BaseM, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        constexpr int pPadded = staticPaddingNewSize<p, BaseP, contains<Method, ResizeStrategy::StaticPadding> ? steps : 0>();
        int runningSpace = 0;
        if (n != nPadded || m != mPadded || p != pPadded) {
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(steps, nPadded, mPadded, pPadded, 1 + BaseN * BaseM, 1 + BaseM * BaseP, MulCount);

        StackAllocator<T> allocator(runningSpace);
        if (n != nPadded || m != mPadded || p != pPadded) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, FunctionImpl>(newC, newA, newB, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP, n, m, p, effC, effA, effB, steps, FunctionImpl>(c, a, b, allocator);
        }
    }

    template<int Method, int BaseN, int BaseM, int BaseP, int MulCount, int n, int m, int p, int effC, int effA, int effB, int steps, typename FunctionImpl, typename M1, typename M2>
    auto runAlgorithm(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        auto c = a.template createNew<n, p>();

        if constexpr (contains<Method, Algorithm::LowLevel>) {
            lowLevelRun<Method, BaseN, BaseM, BaseP, MulCount, n, m, p, effC, effA, effB, steps, FunctionImpl>(
                c.data(), a.data(), b.data()
            );
        }
        else if constexpr (contains<Method, Algorithm::MinSpace>) {
            minSpaceRun<Method, BaseN, BaseM, BaseP, MulCount, n, m, p, effC, effA, effB, steps, FunctionImpl>(
                c.data(), a.data(), b.data()
            );
        }
        else if constexpr (contains<Method, Algorithm::LowLevelParallel>) {
            lowLevelParallelRun<Method, BaseN, BaseM, BaseP, MulCount, n, m, p, effC, effA, effB, steps, FunctionImpl>(
                c.data(), a.data(), b.data()
            );
        }
        else { // Automatic
            // Here the algorithm should be chosen based on type and size of matrices
            lowLevelParallelRun<Method, BaseN, BaseM, BaseP, MulCount, n, m, p, effC, effA, effB, steps, FunctionImpl>(
                c.data(), a.data(), b.data()
            );
        }
        
        return c;
    }
}

#endif

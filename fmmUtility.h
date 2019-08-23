#ifndef FMM_UTILITY_H
#define FMM_UTILITY_H

#include "Matrix.h"
#include "StackAllocator.h"
#include "genericArithmeticOperations.h"
#include "avxMul.h"
#include "blasMul.h"
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
        BlockParallel           = 0x08'00,
        Avx                     = 0x10'00,
        AvxParallel             = 0x20'00,
        Blas                    = 0x40'00,
        BaseMulTypeMASK         = 0xff'00
    };
    enum ResizeStrategy {
        ResizeStrategyAutomatic = 0x00'00'00,
        DynamicPeeling          = 0x01'00'00,
        StaticPadding           = 0x02'00'00,
        ResizeStrategyMASK      = 0xff'00'00
    };
    enum BaseMulSize {
        BaseMulSizeAutomatic    = 0x00'00'00'00,
        Normal                  = 0x01'00'00'00,
        Effective               = 0x02'00'00'00,
        BaseMulSizeMASK         = 0xff'00'00'00
    };
    
    template<int Method, Algorithm alg> 
    constexpr int getNewWithAlgorithm = (Method & (BaseMulTypeMASK | ResizeStrategyMASK | BaseMulSizeMASK)) | alg;
}

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
        if constexpr (Method & BaseMulType::Naive)         naiveMul(c, a, b, n, m, p);
        if constexpr (Method & BaseMulType::Parallel)      parallelMul(c, a, b, n, m, p);
        if constexpr (Method & BaseMulType::Block)         blockMul(c, a, b, n, m, p);
        if constexpr (Method & BaseMulType::BlockParallel) blockParallelMul(c, a, b, n, m, p);
        if constexpr (Method & BaseMulType::Avx)           avx::mul(c, a, b, n, m, p);
        if constexpr (Method & BaseMulType::AvxParallel)   avx::parallelMul(c, a, b, n, m, p);
        if constexpr (Method & BaseMulType::Blas)          blas::mul(c, a, b, n, m, p);
    }

    template<int Method, typename T> 
    void baseMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {
        if constexpr (Method & BaseMulType::Naive)         naiveMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (Method & BaseMulType::Parallel)      parallelMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (Method & BaseMulType::Block)         blockMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (Method & BaseMulType::BlockParallel) blockParallelMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (Method & BaseMulType::Avx)           avx::mul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (Method & BaseMulType::AvxParallel)   avx::parallelMul(c, a, b, n, m, p, effC, effA, effB);
        if constexpr (Method & BaseMulType::Blas)          blas::mul(c, a, b, n, m, p, effC, effA, effB);
    }

    template<int Method, typename T>
    void baseMul(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB, StackAllocator<T>& allocator) {
        if (effC != p || effA != m || effB != p) {
            if constexpr (Method & BaseMulSize::Effective) {
                baseMul<Method>(c, a, b, n, m, p, effC, effA, effB);
            } else if constexpr (Method & BaseMulSize::Normal) {
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


    template<int Method, int BaseN, int BaseM, int BaseP, typename T, typename Function>
    void nextStep(Function function, T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
        if (steps <= 0) {
            baseMul<Method>(c, a, b, n, m, p, effC, effA, effB, allocator);
            return;
        }
        if constexpr (Method & ResizeStrategy::DynamicPeeling) {
            auto[nPeeled, mPeeled, pPeeled] = peelingSizes<BaseN, BaseM, BaseP>(n, m, p);
            function(c, a, b, nPeeled, mPeeled, pPeeled, effC, effA, effB, steps, allocator);
            peelingMul<BaseN, BaseM, BaseP>(c, a, b, n, m, p, effC, effA, effB);
        } else {
            function(c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }


    // Min-Space
    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void minSpaceRecursive(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator);
    
    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void minSpaceRun(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps) {
        auto paddingSizesOpt = staticPaddingNewSizes(n, m, p, BaseN, BaseM, BaseP, steps);
        int runningSpace = 0;
        int nPadded = n;
        int mPadded = m;
        int pPadded = p;
        if ((Method & ResizeStrategy::StaticPadding) && paddingSizesOpt) {
            auto& paddingSizes = *paddingSizesOpt;
            nPadded = paddingSizes[0];
            mPadded = paddingSizes[1];
            pPadded = paddingSizes[2];
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(steps, nPadded, mPadded, pPadded, 1, 1, 1);

        if (Method & BaseMulSize::Normal && !(steps == 0 && effC==p && effA==m && effB==p)) {
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
        if ((Method & ResizeStrategy::StaticPadding) && paddingSizesOpt) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, newC, newA, newB, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }


    // Low-level
    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void lowLevelRecursive(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator);

    template<int Method, int BaseN, int BaseM, int BaseP, typename T> 
    void lowLevelRun(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps) {
        auto paddingSizesOpt = staticPaddingNewSizes(n, m, p, BaseN, BaseM, BaseP, steps);
        int runningSpace = 0;
        int nPadded = n;
        int mPadded = m;
        int pPadded = p;
        if ((Method & ResizeStrategy::StaticPadding) && paddingSizesOpt) {
            auto& paddingSizes = *paddingSizesOpt;
            nPadded = paddingSizes[0];
            mPadded = paddingSizes[1];
            pPadded = paddingSizes[2];
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(steps, nPadded, mPadded, pPadded, 5, 5, 7);

        if (Method & BaseMulSize::Normal && steps == 0 && (effC!=p || effA!=m || effB!=p)) {
            int nn = sizeAfterLastStep<BaseN>(nPadded, steps);
            int mm = sizeAfterLastStep<BaseM>(mPadded, steps);
            int pp = sizeAfterLastStep<BaseP>(pPadded, steps);
            runningSpace += std::max(StackAllocator<T>::Allign(nn*mm), StackAllocator<T>::Allign(mm*pp)) + StackAllocator<T>::Allign(nn*pp);
        }

        StackAllocator<T> allocator(runningSpace);
        if ((Method & ResizeStrategy::StaticPadding) && paddingSizesOpt) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, newC, newA, newB, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }


    // Parallel Low-Level
    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void lowLevelParallelRecursive(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator);

    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void lowLevelParallelRun(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps) {
        //constexpr int UpdatedMethod = 
        
        auto paddingSizesOpt = staticPaddingNewSizes(n, m, p, BaseN, BaseM, BaseP, steps);
        int runningSpace = 0;
        int nPadded = n;
        int mPadded = m;
        int pPadded = p;
        if ((Method & ResizeStrategy::StaticPadding) && paddingSizesOpt) {
            auto& paddingSizes = *paddingSizesOpt;
            nPadded = paddingSizes[0];
            mPadded = paddingSizes[1];
            pPadded = paddingSizes[2];
            runningSpace += additionalStaticPaddingSize<T>(nPadded, mPadded, pPadded);
        }
        runningSpace += additionalRunningSpaceSize<BaseN, BaseM, BaseP, T>(std::min(1, steps), nPadded, mPadded, pPadded, 5, 5, 7);

        StackAllocator<T> allocator(runningSpace);
        if ((Method & ResizeStrategy::StaticPadding) && paddingSizesOpt) {
            auto[newC, newA, newB] = allocateAndCopy(c, a, b, n, m, p, effA, effB, nPadded, mPadded, pPadded, allocator);
            nextStep<Method, BaseN, BaseM, BaseP>(lowLevelParallelRecursive<Method, BaseN, BaseM, BaseP, T>, newC, newA, newB, nPadded, mPadded, pPadded, pPadded, mPadded, pPadded, steps, allocator);
            copyAndDeallocate(newC, newA, newB, n, p, nPadded, mPadded, pPadded, c, effC, allocator);
        } else {
            nextStep<Method, BaseN, BaseM, BaseP>(lowLevelParallelRecursive<Method, BaseN, BaseM, BaseP, T>, c, a, b, n, m, p, effC, effA, effB, steps, allocator);
        }
    }

    template<int Method, int BaseN, int BaseM, int BaseP, typename M1, typename M2>
    auto runAlgorithm(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        auto c = a.createNew(a.rowCount(), b.columnCount());

        if constexpr (Method & Algorithm::LowLevel) {
            lowLevelRun<Method, BaseN, BaseM, BaseP>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        else if constexpr (Method & Algorithm::MinSpace) {
            minSpaceRun<Method, BaseN, BaseM, BaseP>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        else if constexpr (Method & Algorithm::LowLevelParallel) {
            lowLevelParallelRun<Method, BaseN, BaseM, BaseP>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        else { // Automatic
            // Here the algorithm should be chosen based on type and size of matrices
            lowLevelParallelRun<Method, BaseN, BaseM, BaseP>(
                c.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(),
                c.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount(), steps
            );
        }
        
        return c;
    }
}

// Min-Space <2,2,2:7> (Strassen)
namespace fmm::detail {
    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void minSpaceRecursive(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
        using namespace ArithmeticOperation;

        auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

        auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
        auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
        auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

        auto tempA = allocator.alloc(dn * dm);
        auto tempB = allocator.alloc(dm * dp);
        auto tempC = allocator.alloc(dn * dp);

        operationEff<Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][1]);
        operationEff<Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, dC[0][0], tempA, tempB, dn, dm, dp, effC, dm, dp, steps - 1, allocator);
        operationEff<Assign>(dn, dp, effC, effC, dC[1][1], dC[0][0]);

        operationEff<Add>(dn, dm, dm, effA, tempA, dA[1][0], dA[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, dC[1][0], tempA, dB[0][0], dn, dm, dp, effC, dm, effB, steps - 1, allocator);
        operationEff<SubAssign>(dn, dp, effC, effC, dC[1][1], dC[1][0]);

        operationEff<Sub>(dm, dp, dp, effB, tempB, dB[0][1], dB[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, dC[0][1], dA[0][0], tempB, dn, dm, dp, effC, effA, dp, steps - 1, allocator);
        operationEff<AddAssign>(dn, dp, effC, effC, dC[1][1], dC[0][1]);

        operationEff<Sub>(dm, dp, dp, effB, tempB, dB[1][0], dB[0][0]);
        nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, tempC, dA[1][1], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
        operationEff<AddAssign>(dn, dp, effC, dp, dC[0][0], tempC);
        operationEff<AddAssign>(dn, dp, effC, dp, dC[1][0], tempC);

        operationEff<Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, tempC, tempA, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
        operationEff<SubAssign>(dn, dp, effC, dp, dC[0][0], tempC);
        operationEff<AddAssign>(dn, dp, effC, dp, dC[0][1], tempC);

        operationEff<Sub>(dn, dm, dm, effA, tempA, dA[1][0], dA[0][0]);
        operationEff<Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        operationEff<AddAssign>(dn, dp, effC, dp, dC[1][1], tempC);

        operationEff<Sub>(dn, dm, dm, effA, tempA, dA[0][1], dA[1][1]);
        operationEff<Add>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(minSpaceRecursive<Method, BaseN, BaseM, BaseP, T>, tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        operationEff<AddAssign>(dn, dp, effC, dp, dC[0][0], tempC);

        allocator.dealloc(tempC, dn*dp);
        allocator.dealloc(tempB, dm*dp);
        allocator.dealloc(tempA, dn*dm);
    }
}

// Low-Level <2,2,2:7> (Strassen)
namespace fmm::detail {
    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void lowLevelRecursive(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
        using namespace ArithmeticOperation;

        auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

        auto dA = divide<BaseN, BaseM>(a, n, m, effA, allocator);
        auto dB = divide<BaseM, BaseP>(b, m, p, effB, allocator);
        auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

        auto m1 = allocator.alloc(dn * dp);
        auto m1_a = operationEffWithAlloc<Add>(allocator, dn, dm, dm, dA[0][0], dA[1][1]);
        auto m1_b = operationEffWithAlloc<Add>(allocator, dm, dp, dp, dB[0][0], dB[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, m1, m1_a, m1_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        allocator.dealloc(m1_b, dm*dp);
        allocator.dealloc(m1_a, dn*dm);

        auto m2 = allocator.alloc(dn * dp);
        auto m2_a = operationEffWithAlloc<Add>(allocator, dn, dm, dm, dA[1][0], dA[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, m2, m2_a, dB[0][0], dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        allocator.dealloc(m2_a, dn*dm);

        auto m3 = allocator.alloc(dn * dp);
        auto m3_b = operationEffWithAlloc<Sub>(allocator, dm, dp, dp, dB[0][1], dB[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, m3, dA[0][0], m3_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        allocator.dealloc(m3_b, dm*dp);

        auto m4 = allocator.alloc(dn * dp);
        auto m4_b = operationEffWithAlloc<Sub>(allocator, dm, dp, dp, dB[1][0], dB[0][0]);
        nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, m4, dA[1][1], m4_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        allocator.dealloc(m4_b, dm*dp);

        auto m5 = allocator.alloc(dn * dp);
        auto m5_a = operationEffWithAlloc<Add>(allocator, dn, dm, dm, dA[0][0], dA[0][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, m5, m5_a, dB[1][1], dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        allocator.dealloc(m5_a, dn*dm);

        auto m6 = allocator.alloc(dn * dp);
        auto m6_a = operationEffWithAlloc<Sub>(allocator, dn, dm, dm, dA[1][0], dA[0][0]);
        auto m6_b = operationEffWithAlloc<Add>(allocator, dm, dp, dp, dB[0][0], dB[0][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, m6, m6_a, m6_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        allocator.dealloc(m6_b, dm*dp);
        allocator.dealloc(m6_a, dn*dm);

        auto m7 = allocator.alloc(dn * dp);
        auto m7_a = operationEffWithAlloc<Sub>(allocator, dn, dm, dm, dA[0][1], dA[1][1]);
        auto m7_b = operationEffWithAlloc<Add>(allocator, dm, dp, dp, dB[1][0], dB[1][1]);
        nextStep<Method, BaseN, BaseM, BaseP>(lowLevelRecursive<Method, BaseN, BaseM, BaseP, T>, m7, m7_a, m7_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
        allocator.dealloc(m7_b, dm*dp);
        allocator.dealloc(m7_a, dn*dm);

        operationEff<Assign, Add, Sub, Add>(dn, dp, effC, dp, dC[0][0], m1, m4, m5, m7);
        operationEff<Assign, Add>(dn, dp, effC, dp, dC[0][1], m3, m5);
        operationEff<Assign, Add>(dn, dp, effC, dp, dC[1][0], m2, m4);
        operationEff<Assign, Add, Sub, Add>(dn, dp, effC, dp, dC[1][1], m1, m3, m2, m6);

        allocator.dealloc(m7, dn*dp);
        allocator.dealloc(m6, dn*dp);
        allocator.dealloc(m5, dn*dp);
        allocator.dealloc(m4, dn*dp);
        allocator.dealloc(m3, dn*dp);
        allocator.dealloc(m2, dn*dp);
        allocator.dealloc(m1, dn*dp);
        allocator.dealloc(dB[1][1], dm*dp);
        allocator.dealloc(dB[1][0], dm*dp);
        allocator.dealloc(dB[0][1], dm*dp);
        allocator.dealloc(dB[0][0], dm*dp);
        allocator.dealloc(dA[1][1], dn*dm);
        allocator.dealloc(dA[1][0], dn*dm);
        allocator.dealloc(dA[0][1], dn*dm);
        allocator.dealloc(dA[0][0], dn*dm);
    }
}


// Parallel Low-Level <2,2,2:7> (Strassen)
namespace fmm::detail {
    template<int Method, int BaseN, int BaseM, int BaseP, typename T>
    void lowLevelParallelRecursive(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
        using namespace ArithmeticOperation;

        auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

        auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
        auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
        auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

        auto m1 = allocator.alloc(dn * dp);
        auto m2 = allocator.alloc(dn * dp);
        auto m3 = allocator.alloc(dn * dp);
        auto m4 = allocator.alloc(dn * dp);
        auto m5 = allocator.alloc(dn * dp);
        auto m6 = allocator.alloc(dn * dp);
        auto m7 = allocator.alloc(dn * dp);
        auto m1_a = allocator.alloc(dn*dm);
        auto m1_b = allocator.alloc(dm*dp);
        auto m2_a = allocator.alloc(dn*dm);
        auto m3_b = allocator.alloc(dm*dp);
        auto m4_b = allocator.alloc(dm*dp);
        auto m5_a = allocator.alloc(dn*dm);
        auto m6_a = allocator.alloc(dn*dm);
        auto m6_b = allocator.alloc(dm*dp);
        auto m7_a = allocator.alloc(dn*dm);
        auto m7_b = allocator.alloc(dm*dp);

        ThreadPool pool;

        pool.addTask([=]() {
            operationEff<Add>(dn, dm, dm, effA, m1_a, dA[0][0], dA[1][1]);
            operationEff<Add>(dm, dp, dp, effB, m1_b, dB[0][0], dB[1][1]);
            minSpaceRun<Method, BaseN, BaseM, BaseP>(m1, m1_a, m1_b, dn, dm, dp, dp, dm, dp, steps - 1);
        });
        pool.addTask([=]() {
            operationEff<Add>(dn, dm, dm, effA, m2_a, dA[1][0], dA[1][1]);
            minSpaceRun<Method, BaseN, BaseM, BaseP>(m2, m2_a, dB[0][0], dn, dm, dp, dp, dm, effB, steps - 1);
        });
        pool.addTask([=]() {
            operationEff<Sub>(dm, dp, dp, effB, m3_b, dB[0][1], dB[1][1]);
            minSpaceRun<Method, BaseN, BaseM, BaseP>(m3, dA[0][0], m3_b, dn, dm, dp, dp, effA, dp, steps - 1);
        });
        pool.addTask([=]() {
            operationEff<Sub>(dm, dp, dp, effB, m4_b, dB[1][0], dB[0][0]);
            minSpaceRun<Method, BaseN, BaseM, BaseP>(m4, dA[1][1], m4_b, dn, dm, dp, dp, effA, dp, steps - 1);
        });
        pool.addTask([=]() {
            operationEff<Add>(dn, dm, dm, effA, m5_a, dA[0][0], dA[0][1]);
            minSpaceRun<Method, BaseN, BaseM, BaseP>(m5, m5_a, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1);
        });
        pool.addTask([=]() {
            operationEff<Sub>(dn, dm, dm, effA, m6_a, dA[1][0], dA[0][0]);
            operationEff<Add>(dm, dp, dp, effB, m6_b, dB[0][0], dB[0][1]);
            minSpaceRun<Method, BaseN, BaseM, BaseP>(m6, m6_a, m6_b, dn, dm, dp, dp, dm, dp, steps - 1);
        });
        pool.addTask([=]() {
            operationEff<Sub>(dn, dm, dm, effA, m7_a, dA[0][1], dA[1][1]);
            operationEff<Add>(dm, dp, dp, effB, m7_b, dB[1][0], dB[1][1]);
            minSpaceRun<Method, BaseN, BaseM, BaseP>(m7, m7_a, m7_b, dn, dm, dp, dp, dm, dp, steps - 1);
        });

        pool.completeTasksAndStop();

        operationEff<Assign, Add, Sub, Add>(dn, dp, effC, dp, dC[0][0], m1, m4, m5, m7);
        operationEff<Assign, Add>(dn, dp, effC, dp, dC[0][1], m3, m5);
        operationEff<Assign, Add>(dn, dp, effC, dp, dC[1][0], m2, m4);
        operationEff<Assign, Add, Sub, Add>(dn, dp, effC, dp, dC[1][1], m1, m3, m2, m6);

        allocator.dealloc(m7_b, dm*dp);
        allocator.dealloc(m7_a, dn*dm);
        allocator.dealloc(m6_b, dm*dp);
        allocator.dealloc(m6_a, dn*dm);
        allocator.dealloc(m5_a, dn*dm);
        allocator.dealloc(m4_b, dm*dp);
        allocator.dealloc(m3_b, dm*dp);
        allocator.dealloc(m2_a, dn*dm);
        allocator.dealloc(m1_b, dm*dp);
        allocator.dealloc(m1_a, dn*dm);
        allocator.dealloc(m7, dn*dp);
        allocator.dealloc(m6, dn*dp);
        allocator.dealloc(m5, dn*dp);
        allocator.dealloc(m4, dn*dp);
        allocator.dealloc(m3, dn*dp);
        allocator.dealloc(m2, dn*dp);
        allocator.dealloc(m1, dn*dp);
    }
}


namespace fmm {
    template<int Method, typename M1, typename M2>
    auto minSpaceStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 2, 2, 2>(a, b, steps);
    }

    template<int Method, typename M1, typename M2>
    auto lowLevelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevel>, 2, 2, 2>(a, b, steps);
    }

    template<int Method, typename M1, typename M2>
    auto lowLevelParallelStrassen(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevelParallel>, 2, 2, 2>(a, b, steps);
    }
}

#endif

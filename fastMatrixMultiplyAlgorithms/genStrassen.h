// this file was generated using fastMatrixMultiplyAlgorithms/generator
#ifndef GEN_FMM_2_2_2_7_H
#define GEN_FMM_2_2_2_7_H
#include "../fmmUtility.h"
#include "../ThreadPool.h"

namespace fmm::detail {
    struct GenStrassenRecursiveMinSpace {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 2;
            constexpr int BaseP = 2;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto tempA = allocator.alloc(dn * dm);
            auto tempB = allocator.alloc(dm * dp);
            auto tempC = allocator.alloc(dn * dp);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][1]);
            operationEff<Assign, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<Assign, Assign>(dn, dp, dp, effC, tempC, dC[0][0], dC[1][1]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveMinSpace>(tempC, tempA, dB[0][0], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<Assign, SubAssign>(dn, dp, dp, effC, tempC, dC[1][0], dC[1][1]);

            operationEff<Assign, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveMinSpace>(tempC, dA[0][0], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArg<Assign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[1][1]);

            operationEff<Assign, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[1][0]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveMinSpace>(tempC, dA[1][1], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[1][0]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveMinSpace>(tempC, tempA, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<SubAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][1]);

            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][0]);
            operationEff<Assign, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[1][1]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][1], dA[1][1]);
            operationEff<Assign, Add, Add>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[0][0]);

            allocator.dealloc(tempC, dn*dp);
            allocator.dealloc(tempB, dm*dp);
            allocator.dealloc(tempA, dn*dm);
        }
    };
}

namespace fmm::detail {
    struct GenStrassenRecursiveLowLevel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 2;
            constexpr int BaseP = 2;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto m1 = allocator.alloc(dn * dp);
            auto m1A = allocator.alloc(dn*dm);
            auto m1B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m1A, dA[0][0], dA[1][1]);
            operationEff<Assign, Add, Add>(dm, dp, dp, effB, m1B, dB[0][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveLowLevel>(m1, m1A, m1B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m1B, dm*dp);
            allocator.dealloc(m1A, dn*dm);

            auto m2 = allocator.alloc(dn * dp);
            auto m2A = allocator.alloc(dn*dm);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m2A, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveLowLevel>(m2, m2A, dB[0][0], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m2A, dn*dm);

            auto m3 = allocator.alloc(dn * dp);
            auto m3B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Sub>(dm, dp, dp, effB, m3B, dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveLowLevel>(m3, dA[0][0], m3B, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m3B, dm*dp);

            auto m4 = allocator.alloc(dn * dp);
            auto m4B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add>(dm, dp, dp, effB, m4B, dB[0][0], dB[1][0]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveLowLevel>(m4, dA[1][1], m4B, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m4B, dm*dp);

            auto m5 = allocator.alloc(dn * dp);
            auto m5A = allocator.alloc(dn*dm);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m5A, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveLowLevel>(m5, m5A, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m5A, dn*dm);

            auto m6 = allocator.alloc(dn * dp);
            auto m6A = allocator.alloc(dn*dm);
            auto m6B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m6A, dA[0][0], dA[1][0]);
            operationEff<Assign, Add, Add>(dm, dp, dp, effB, m6B, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveLowLevel>(m6, m6A, m6B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m6B, dm*dp);
            allocator.dealloc(m6A, dn*dm);

            auto m7 = allocator.alloc(dn * dp);
            auto m7A = allocator.alloc(dn*dm);
            auto m7B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m7A, dA[0][1], dA[1][1]);
            operationEff<Assign, Add, Add>(dm, dp, dp, effB, m7B, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenStrassenRecursiveLowLevel>(m7, m7A, m7B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m7B, dm*dp);
            allocator.dealloc(m7A, dn*dm);

            operationEff<Assign, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[0][0], m1, m4, m5, m7);
            operationEff<Assign, Add, Add>(dn, dp, effC, dp, dC[0][1], m3, m5);
            operationEff<Assign, Add, Add>(dn, dp, effC, dp, dC[1][0], m2, m4);
            operationEff<Assign, Add, Sub, Add, Add>(dn, dp, effC, dp, dC[1][1], m1, m2, m3, m6);

            allocator.dealloc(m7, dn*dp);
            allocator.dealloc(m6, dn*dp);
            allocator.dealloc(m5, dn*dp);
            allocator.dealloc(m4, dn*dp);
            allocator.dealloc(m3, dn*dp);
            allocator.dealloc(m2, dn*dp);
            allocator.dealloc(m1, dn*dp);
        }
    };
}

namespace fmm::detail {
    struct GenStrassenRecursiveLowLevelParallel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 2;
            constexpr int BaseP = 2;

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
            auto m1A = allocator.alloc(dn*dm);
            auto m1B = allocator.alloc(dm*dp);
            auto m2A = allocator.alloc(dn*dm);
            auto m3B = allocator.alloc(dm*dp);
            auto m4B = allocator.alloc(dm*dp);
            auto m5A = allocator.alloc(dn*dm);
            auto m6A = allocator.alloc(dn*dm);
            auto m6B = allocator.alloc(dm*dp);
            auto m7A = allocator.alloc(dn*dm);
            auto m7B = allocator.alloc(dm*dp);

            ThreadPool pool;

            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m1A, dA[0][0], dA[1][1]);
                operationEff<Assign, Add, Add>(dm, dp, dp, effB, m1B, dB[0][0], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, GenStrassenRecursiveMinSpace>(m1, m1A, m1B, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m2A, dA[1][0], dA[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, GenStrassenRecursiveMinSpace>(m2, m2A, dB[0][0], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Sub>(dm, dp, dp, effB, m3B, dB[0][1], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, GenStrassenRecursiveMinSpace>(m3, dA[0][0], m3B, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add>(dm, dp, dp, effB, m4B, dB[0][0], dB[1][0]);
                minSpaceRun<Method, 2, 2, 2, 7, GenStrassenRecursiveMinSpace>(m4, dA[1][1], m4B, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m5A, dA[0][0], dA[0][1]);
                minSpaceRun<Method, 2, 2, 2, 7, GenStrassenRecursiveMinSpace>(m5, m5A, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m6A, dA[0][0], dA[1][0]);
                operationEff<Assign, Add, Add>(dm, dp, dp, effB, m6B, dB[0][0], dB[0][1]);
                minSpaceRun<Method, 2, 2, 2, 7, GenStrassenRecursiveMinSpace>(m6, m6A, m6B, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m7A, dA[0][1], dA[1][1]);
                operationEff<Assign, Add, Add>(dm, dp, dp, effB, m7B, dB[1][0], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, GenStrassenRecursiveMinSpace>(m7, m7A, m7B, dn, dm, dp, dp, dm, dp, steps - 1);
            });

            pool.completeTasksAndStop();

            operationEff<Assign, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[0][0], m1, m4, m5, m7);
            operationEff<Assign, Add, Add>(dn, dp, effC, dp, dC[0][1], m3, m5);
            operationEff<Assign, Add, Add>(dn, dp, effC, dp, dC[1][0], m2, m4);
            operationEff<Assign, Add, Sub, Add, Add>(dn, dp, effC, dp, dC[1][1], m1, m2, m3, m6);

            allocator.dealloc(m7B, dm*dp);
            allocator.dealloc(m7A, dn*dm);
            allocator.dealloc(m6B, dm*dp);
            allocator.dealloc(m6A, dn*dm);
            allocator.dealloc(m5A, dn*dm);
            allocator.dealloc(m4B, dm*dp);
            allocator.dealloc(m3B, dm*dp);
            allocator.dealloc(m2A, dn*dm);
            allocator.dealloc(m1B, dm*dp);
            allocator.dealloc(m1A, dn*dm);
            allocator.dealloc(m7, dn*dp);
            allocator.dealloc(m6, dn*dp);
            allocator.dealloc(m5, dn*dp);
            allocator.dealloc(m4, dn*dp);
            allocator.dealloc(m3, dn*dp);
            allocator.dealloc(m2, dn*dp);
            allocator.dealloc(m1, dn*dp);
        }
    };
}

namespace fmm {
    template<int Method = 0, typename M1, typename M2>
    auto genStrassenMinSpace(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 2, 2, 2, 7, detail::GenStrassenRecursiveMinSpace>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto genStrassenLowLevel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevel>, 2, 2, 2, 7, detail::GenStrassenRecursiveLowLevel>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto genStrassenLowLevelParallel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevelParallel>, 2, 2, 2, 7, detail::GenStrassenRecursiveLowLevelParallel>(a, b, steps, 5, 5);
    }
}

#endif

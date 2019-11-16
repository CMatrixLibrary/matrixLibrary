// this file was generated using fastMatrixMultiplyAlgorithms/generator
#ifndef GEN_FMM_2_3_3_15_H
#define GEN_FMM_2_3_3_15_H
#include "../fmmUtility.h"
#include "../ThreadPool.h"

namespace fmm::detail {
    struct Gen2x3x3RecursiveMinSpace {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 3;
            constexpr int BaseP = 3;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto tempA = allocator.alloc(dn * dm);
            auto tempB = allocator.alloc(dm * dp);
            auto tempC = allocator.alloc(dn * dp);

            operationEff<Assign, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[0][1], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, dA[1][1], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<Assign, Assign>(dn, dp, dp, effC, tempC, std::pair(-1.000000, dC[0][1]), std::pair(1.000000, dC[1][1]));

            operationEff<Assign, Add, Add, Add>(dm, dp, dp, effB, tempB, dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, dA[0][2], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArg<Assign>(dn, dp, dp, effC, tempC, dC[0][2]);

            operationEff<Assign, Add, Add, Add>(dm, dp, dp, effB, tempB, dB[0][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, dA[1][2], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArg<Assign>(dn, dp, dp, effC, tempC, dC[1][2]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[0][1], dA[1][1]);
            operationEff<Assign, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[0][1]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][0]);
            operationEff<Assign, Add, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1], dB[2][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<Assign>(dn, dp, dp, effC, tempC, dC[1][0]);

            operationEff<Assign, Sub, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, dA[0][0], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<Assign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(-1.000000, dC[0][0]), std::pair(1.000000, dC[1][0]));

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, dB[1][0], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][1]);

            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, tempA, dA[0][1], dA[0][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, dB[1][2], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[0][2], dC[1][1], dC[1][2]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, dB[0][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[1][0], dC[1][1]);

            operationEff<Assign, Sub, Add, Sub, Add>(dn, dm, dm, effA, tempA, dA[0][1], dA[0][2], dA[1][1], dA[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, dB[2][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, SubAssign>(dn, dp, dp, effC, tempC, dC[1][1], dC[1][2]);

            operationEff<Assign, Sub, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][1], dA[0][2], dA[1][1]);
            operationEff<Assign, Sub, Add>(dm, dp, dp, effB, tempB, dB[1][2], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[1][1], dC[1][2]);

            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][0], dA[1][2]);
            operationEff<Assign, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][2], dB[2][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][2], dC[1][0]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][1]);
            operationEff<Assign, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][1], dB[1][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<SubAssign, SubAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[1][0]);

            operationEff<Assign, Sub, Add, Sub, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][2], dA[1][0], dA[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, dB[2][0], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, SubAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][2]);

            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, tempA, dA[1][0], dA[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveMinSpace>(tempC, tempA, dB[0][2], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<SubAssign, AddAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][2], dC[1][0], dC[1][2]);

            allocator.dealloc(tempC, dn*dp);
            allocator.dealloc(tempB, dm*dp);
            allocator.dealloc(tempA, dn*dm);
        }
    };
}

namespace fmm::detail {
    struct Gen2x3x3RecursiveLowLevel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 3;
            constexpr int BaseP = 3;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto m1 = allocator.alloc(dn * dp);
            auto m1B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Add>(dm, dp, dp, effB, m1B, dB[0][1], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m1, dA[1][1], m1B, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m1B, dm*dp);

            auto m2 = allocator.alloc(dn * dp);
            auto m2B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add, Add>(dm, dp, dp, effB, m2B, dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m2, dA[0][2], m2B, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m2B, dm*dp);

            auto m3 = allocator.alloc(dn * dp);
            auto m3B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add, Add>(dm, dp, dp, effB, m3B, dB[0][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m3, dA[1][2], m3B, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m3B, dm*dp);

            auto m4 = allocator.alloc(dn * dp);
            auto m4A = allocator.alloc(dn*dm);
            auto m4B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m4A, dA[0][1], dA[1][1]);
            operationEff<Assign, Sub, Add, Add>(dm, dp, dp, effB, m4B, dB[1][0], dB[1][1], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m4, m4A, m4B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m4B, dm*dp);
            allocator.dealloc(m4A, dn*dm);

            auto m5 = allocator.alloc(dn * dp);
            auto m5A = allocator.alloc(dn*dm);
            auto m5B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m5A, dA[0][0], dA[1][0]);
            operationEff<Assign, Add, Sub, Add>(dm, dp, dp, effB, m5B, dB[0][0], dB[0][1], dB[2][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m5, m5A, m5B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m5B, dm*dp);
            allocator.dealloc(m5A, dn*dm);

            auto m6 = allocator.alloc(dn * dp);
            auto m6B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Sub, Add>(dm, dp, dp, effB, m6B, dB[0][0], dB[0][2], dB[1][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m6, dA[0][0], m6B, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m6B, dm*dp);

            auto m7 = allocator.alloc(dn * dp);
            auto m7A = allocator.alloc(dn*dm);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m7A, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m7, m7A, dB[1][0], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m7A, dn*dm);

            auto m8 = allocator.alloc(dn * dp);
            auto m8A = allocator.alloc(dn*dm);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m8A, dA[0][1], dA[0][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m8, m8A, dB[1][2], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m8A, dn*dm);

            auto m9 = allocator.alloc(dn * dp);
            auto m9A = allocator.alloc(dn*dm);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m9A, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m9, m9A, dB[0][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m9A, dn*dm);

            auto m10 = allocator.alloc(dn * dp);
            auto m10A = allocator.alloc(dn*dm);
            operationEff<Assign, Sub, Add, Sub, Add>(dn, dm, dm, effA, m10A, dA[0][1], dA[0][2], dA[1][1], dA[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m10, m10A, dB[2][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m10A, dn*dm);

            auto m11 = allocator.alloc(dn * dp);
            auto m11A = allocator.alloc(dn*dm);
            auto m11B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Sub>(dn, dm, dm, effA, m11A, dA[0][1], dA[0][2], dA[1][1]);
            operationEff<Assign, Sub, Add>(dm, dp, dp, effB, m11B, dB[1][2], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m11, m11A, m11B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m11B, dm*dp);
            allocator.dealloc(m11A, dn*dm);

            auto m12 = allocator.alloc(dn * dp);
            auto m12A = allocator.alloc(dn*dm);
            auto m12B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effA, m12A, dA[0][0], dA[1][0], dA[1][2]);
            operationEff<Assign, Add, Sub>(dm, dp, dp, effB, m12B, dB[0][2], dB[2][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m12, m12A, m12B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m12B, dm*dp);
            allocator.dealloc(m12A, dn*dm);

            auto m13 = allocator.alloc(dn * dp);
            auto m13A = allocator.alloc(dn*dm);
            auto m13B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m13A, dA[0][0], dA[1][1]);
            operationEff<Assign, Sub, Add>(dm, dp, dp, effB, m13B, dB[0][1], dB[1][0]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m13, m13A, m13B, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m13B, dm*dp);
            allocator.dealloc(m13A, dn*dm);

            auto m14 = allocator.alloc(dn * dp);
            auto m14A = allocator.alloc(dn*dm);
            operationEff<Assign, Sub, Add, Sub, Add>(dn, dm, dm, effA, m14A, dA[0][0], dA[0][2], dA[1][0], dA[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m14, m14A, dB[2][0], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m14A, dn*dm);

            auto m15 = allocator.alloc(dn * dp);
            auto m15A = allocator.alloc(dn*dm);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m15A, dA[1][0], dA[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen2x3x3RecursiveLowLevel>(m15, m15A, dB[0][2], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m15A, dn*dm);

            operationEff<Assign, Sub, Add, Sub, Add, Sub>(dn, dp, effC, dp, dC[0][0], m6, m7, m12, m14, m15);
            operationEff<Assign, Sub, Add, Add, Add, Add, Sub>(dn, dp, effC, dp, dC[0][1], m1, m4, m7, m8, m11, m13);
            operationEff<Assign, Add, Sub, Add, Sub, Add>(dn, dp, effC, dp, dC[0][2], m2, m8, m12, m14, m15);
            operationEff<Assign, Add, Add, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[1][0], m5, m6, m9, m12, m13, m15);
            operationEff<Assign, Add, Sub, Add, Add, Sub>(dn, dp, effC, dp, dC[1][1], m1, m8, m9, m10, m11);
            operationEff<Assign, Add, Add, Sub, Add, Sub>(dn, dp, effC, dp, dC[1][2], m3, m8, m10, m11, m15);

            allocator.dealloc(m15, dn*dp);
            allocator.dealloc(m14, dn*dp);
            allocator.dealloc(m13, dn*dp);
            allocator.dealloc(m12, dn*dp);
            allocator.dealloc(m11, dn*dp);
            allocator.dealloc(m10, dn*dp);
            allocator.dealloc(m9, dn*dp);
            allocator.dealloc(m8, dn*dp);
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
    struct Gen2x3x3RecursiveLowLevelParallel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 3;
            constexpr int BaseP = 3;

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
            auto m8 = allocator.alloc(dn * dp);
            auto m9 = allocator.alloc(dn * dp);
            auto m10 = allocator.alloc(dn * dp);
            auto m11 = allocator.alloc(dn * dp);
            auto m12 = allocator.alloc(dn * dp);
            auto m13 = allocator.alloc(dn * dp);
            auto m14 = allocator.alloc(dn * dp);
            auto m15 = allocator.alloc(dn * dp);
            auto m1B = allocator.alloc(dm*dp);
            auto m2B = allocator.alloc(dm*dp);
            auto m3B = allocator.alloc(dm*dp);
            auto m4A = allocator.alloc(dn*dm);
            auto m4B = allocator.alloc(dm*dp);
            auto m5A = allocator.alloc(dn*dm);
            auto m5B = allocator.alloc(dm*dp);
            auto m6B = allocator.alloc(dm*dp);
            auto m7A = allocator.alloc(dn*dm);
            auto m8A = allocator.alloc(dn*dm);
            auto m9A = allocator.alloc(dn*dm);
            auto m10A = allocator.alloc(dn*dm);
            auto m11A = allocator.alloc(dn*dm);
            auto m11B = allocator.alloc(dm*dp);
            auto m12A = allocator.alloc(dn*dm);
            auto m12B = allocator.alloc(dm*dp);
            auto m13A = allocator.alloc(dn*dm);
            auto m13B = allocator.alloc(dm*dp);
            auto m14A = allocator.alloc(dn*dm);
            auto m15A = allocator.alloc(dn*dm);

            ThreadPool pool;

            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Add>(dm, dp, dp, effB, m1B, dB[0][1], dB[1][1], dB[1][2]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m1, dA[1][1], m1B, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add, Add>(dm, dp, dp, effB, m2B, dB[1][2], dB[2][0], dB[2][2]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m2, dA[0][2], m2B, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add, Add>(dm, dp, dp, effB, m3B, dB[0][2], dB[2][1], dB[2][2]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m3, dA[1][2], m3B, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m4A, dA[0][1], dA[1][1]);
                operationEff<Assign, Sub, Add, Add>(dm, dp, dp, effB, m4B, dB[1][0], dB[1][1], dB[2][1]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m4, m4A, m4B, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m5A, dA[0][0], dA[1][0]);
                operationEff<Assign, Add, Sub, Add>(dm, dp, dp, effB, m5B, dB[0][0], dB[0][1], dB[2][0]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m5, m5A, m5B, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Sub, Add>(dm, dp, dp, effB, m6B, dB[0][0], dB[0][2], dB[1][0]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m6, dA[0][0], m6B, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m7A, dA[0][0], dA[0][1]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m7, m7A, dB[1][0], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m8A, dA[0][1], dA[0][2]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m8, m8A, dB[1][2], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m9A, dA[1][0], dA[1][1]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m9, m9A, dB[0][1], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Sub, Add>(dn, dm, dm, effA, m10A, dA[0][1], dA[0][2], dA[1][1], dA[1][2]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m10, m10A, dB[2][1], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Sub>(dn, dm, dm, effA, m11A, dA[0][1], dA[0][2], dA[1][1]);
                operationEff<Assign, Sub, Add>(dm, dp, dp, effB, m11B, dB[1][2], dB[2][1]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m11, m11A, m11B, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effA, m12A, dA[0][0], dA[1][0], dA[1][2]);
                operationEff<Assign, Add, Sub>(dm, dp, dp, effB, m12B, dB[0][2], dB[2][0]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m12, m12A, m12B, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m13A, dA[0][0], dA[1][1]);
                operationEff<Assign, Sub, Add>(dm, dp, dp, effB, m13B, dB[0][1], dB[1][0]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m13, m13A, m13B, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Sub, Add>(dn, dm, dm, effA, m14A, dA[0][0], dA[0][2], dA[1][0], dA[1][2]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m14, m14A, dB[2][0], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m15A, dA[1][0], dA[1][2]);
                minSpaceRun<Method, 2, 3, 3, 15, Gen2x3x3RecursiveMinSpace>(m15, m15A, dB[0][2], dn, dm, dp, dp, dm, effB, steps - 1);
            });

            pool.completeTasksAndStop();

            operationEff<Assign, Sub, Add, Sub, Add, Sub>(dn, dp, effC, dp, dC[0][0], m6, m7, m12, m14, m15);
            operationEff<Assign, Sub, Add, Add, Add, Add, Sub>(dn, dp, effC, dp, dC[0][1], m1, m4, m7, m8, m11, m13);
            operationEff<Assign, Add, Sub, Add, Sub, Add>(dn, dp, effC, dp, dC[0][2], m2, m8, m12, m14, m15);
            operationEff<Assign, Add, Add, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[1][0], m5, m6, m9, m12, m13, m15);
            operationEff<Assign, Add, Sub, Add, Add, Sub>(dn, dp, effC, dp, dC[1][1], m1, m8, m9, m10, m11);
            operationEff<Assign, Add, Add, Sub, Add, Sub>(dn, dp, effC, dp, dC[1][2], m3, m8, m10, m11, m15);

            allocator.dealloc(m15A, dn*dm);
            allocator.dealloc(m14A, dn*dm);
            allocator.dealloc(m13B, dm*dp);
            allocator.dealloc(m13A, dn*dm);
            allocator.dealloc(m12B, dm*dp);
            allocator.dealloc(m12A, dn*dm);
            allocator.dealloc(m11B, dm*dp);
            allocator.dealloc(m11A, dn*dm);
            allocator.dealloc(m10A, dn*dm);
            allocator.dealloc(m9A, dn*dm);
            allocator.dealloc(m8A, dn*dm);
            allocator.dealloc(m7A, dn*dm);
            allocator.dealloc(m6B, dm*dp);
            allocator.dealloc(m5B, dm*dp);
            allocator.dealloc(m5A, dn*dm);
            allocator.dealloc(m4B, dm*dp);
            allocator.dealloc(m4A, dn*dm);
            allocator.dealloc(m3B, dm*dp);
            allocator.dealloc(m2B, dm*dp);
            allocator.dealloc(m1B, dm*dp);
            allocator.dealloc(m15, dn*dp);
            allocator.dealloc(m14, dn*dp);
            allocator.dealloc(m13, dn*dp);
            allocator.dealloc(m12, dn*dp);
            allocator.dealloc(m11, dn*dp);
            allocator.dealloc(m10, dn*dp);
            allocator.dealloc(m9, dn*dp);
            allocator.dealloc(m8, dn*dp);
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
    auto gen2x3x3MinSpace(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 2, 3, 3, 15, detail::Gen2x3x3RecursiveMinSpace>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto gen2x3x3LowLevel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevel>, 2, 3, 3, 15, detail::Gen2x3x3RecursiveLowLevel>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto gen2x3x3LowLevelParallel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevelParallel>, 2, 3, 3, 15, detail::Gen2x3x3RecursiveLowLevelParallel>(a, b, steps, 11, 9);
    }
}

#endif

// this file was generated using fastMatrixMultiplyAlgorithms/generator
#ifndef GEN_FMM_3_3_3_23_H
#define GEN_FMM_3_3_3_23_H
#include "../fmmUtility.h"
#include "../ThreadPool.h"

namespace fmm::detail {
    struct GenFastMul3x3RecursiveMinSpace {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 3;
            constexpr int BaseM = 3;
            constexpr int BaseP = 3;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto tempA = allocator.alloc(dn * dm);
            auto tempB = allocator.alloc(dm * dp);
            auto tempC = allocator.alloc(dn * dp);

            operationEff<Assign, Add, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][1], dA[0][2], dA[1][0], dA[1][1], dA[2][1], dA[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<Assign>(dn, dp, dp, effC, tempC, dC[0][1]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][0]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<Assign, Assign>(dn, dp, dp, effC, tempC, dC[1][0], dC[1][1]);

            operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][1], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[1][1], tempB, dn, dm, dp, dp, effA, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[1][0]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][0], dA[1][1]);
            operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[1][0], dC[1][1]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[1][0], dA[1][1]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[1][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[0][0], dB[0][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<Assign, AddAssign, Assign, AddAssign, AddAssign, Assign, Assign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][1], dC[0][2], dC[1][0], dC[1][1], dC[2][0], dC[2][2]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[2][0], dA[2][1]);
            operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][2], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[2][0], dC[2][2]);

            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[2][0]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, tempB, dB[0][2], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[2][0], dC[2][2]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[2][0], dA[2][1]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[2][2]);

            operationEff<Assign, Add, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][1], dA[0][2], dA[1][1], dA[1][2], dA[2][0], dA[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, dB[1][2], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[0][2]);

            operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][2], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[2][1], tempB, dn, dm, dp, dp, effA, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[2][0]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][2], dA[2][1], dA[2][2]);
            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, Assign>(dn, dp, dp, effC, tempC, dC[0][1], dC[2][0], dC[2][1]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][2], dA[2][2]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][1], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[2][0], dC[2][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[0][2], dB[2][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign, AddAssign, Assign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][1], dC[0][2], dC[1][0], dC[1][2], dC[2][0], dC[2][1]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[2][1], dA[2][2]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[2][1]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][2], dA[1][1], dA[1][2]);
            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[1][0], dC[1][2]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][2], dA[1][2]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][2], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[1][0], dC[1][2]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[1][1], dA[1][2]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[1][2]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[0][1], dB[1][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[0][0]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[1][2], dB[2][1], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[1][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[1][0], dB[0][2], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[1][2]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[2][0], dB[0][1], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[2][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveMinSpace>(tempC, dA[2][2], dB[2][2], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[2][2]);

            allocator.dealloc(tempC, dn*dp);
            allocator.dealloc(tempB, dm*dp);
            allocator.dealloc(tempA, dn*dm);
        }
    };
}

namespace fmm::detail {
    struct GenFastMul3x3RecursiveLowLevel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 3;
            constexpr int BaseM = 3;
            constexpr int BaseP = 3;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto m1 = allocator.alloc(dn * dp);
            auto m1A = allocator.alloc(dn*dm);
            operationEff<Assign, Add, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, m1A, dA[0][0], dA[0][1], dA[0][2], dA[1][0], dA[1][1], dA[2][1], dA[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m1, m1A, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m1A, dn*dm);

            auto m2 = allocator.alloc(dn * dp);
            auto m2A = allocator.alloc(dn*dm);
            auto m2B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m2A, dA[0][0], dA[1][0]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m2B, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m2, m2A, m2B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m2B, dm*dp);
            allocator.dealloc(m2A, dn*dm);

            auto m3 = allocator.alloc(dn * dp);
            auto m3B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, m3B, dB[0][0], dB[0][1], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m3, dA[1][1], m3B, dn, dm, dp, dp, effA, dm, steps - 1, allocator);
            allocator.dealloc(m3B, dm*dp);

            auto m4 = allocator.alloc(dn * dp);
            auto m4A = allocator.alloc(dn*dm);
            auto m4B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m4A, dA[0][0], dA[1][0], dA[1][1]);
            operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, m4B, dB[0][0], dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m4, m4A, m4B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m4B, dm*dp);
            allocator.dealloc(m4A, dn*dm);

            auto m5 = allocator.alloc(dn * dp);
            auto m5A = allocator.alloc(dn*dm);
            auto m5B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m5A, dA[1][0], dA[1][1]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m5B, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m5, m5A, m5B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m5B, dm*dp);
            allocator.dealloc(m5A, dn*dm);

            auto m6 = allocator.alloc(dn * dp);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m6, dA[0][0], dB[0][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);

            auto m7 = allocator.alloc(dn * dp);
            auto m7A = allocator.alloc(dn*dm);
            auto m7B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m7A, dA[0][0], dA[2][0], dA[2][1]);
            operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, m7B, dB[0][0], dB[0][2], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m7, m7A, m7B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m7B, dm*dp);
            allocator.dealloc(m7A, dn*dm);

            auto m8 = allocator.alloc(dn * dp);
            auto m8A = allocator.alloc(dn*dm);
            auto m8B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m8A, dA[0][0], dA[2][0]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, m8B, dB[0][2], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m8, m8A, m8B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m8B, dm*dp);
            allocator.dealloc(m8A, dn*dm);

            auto m9 = allocator.alloc(dn * dp);
            auto m9A = allocator.alloc(dn*dm);
            auto m9B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m9A, dA[2][0], dA[2][1]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m9B, dB[0][0], dB[0][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m9, m9A, m9B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m9B, dm*dp);
            allocator.dealloc(m9A, dn*dm);

            auto m10 = allocator.alloc(dn * dp);
            auto m10A = allocator.alloc(dn*dm);
            operationEff<Assign, Add, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, m10A, dA[0][0], dA[0][1], dA[0][2], dA[1][1], dA[1][2], dA[2][0], dA[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m10, m10A, dB[1][2], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m10A, dn*dm);

            auto m11 = allocator.alloc(dn * dp);
            auto m11B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, m11B, dB[0][0], dB[0][2], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m11, dA[2][1], m11B, dn, dm, dp, dp, effA, dm, steps - 1, allocator);
            allocator.dealloc(m11B, dm*dp);

            auto m12 = allocator.alloc(dn * dp);
            auto m12A = allocator.alloc(dn*dm);
            auto m12B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m12A, dA[0][2], dA[2][1], dA[2][2]);
            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, m12B, dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m12, m12A, m12B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m12B, dm*dp);
            allocator.dealloc(m12A, dn*dm);

            auto m13 = allocator.alloc(dn * dp);
            auto m13A = allocator.alloc(dn*dm);
            auto m13B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m13A, dA[0][2], dA[2][2]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, m13B, dB[1][1], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m13, m13A, m13B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m13B, dm*dp);
            allocator.dealloc(m13A, dn*dm);

            auto m14 = allocator.alloc(dn * dp);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m14, dA[0][2], dB[2][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);

            auto m15 = allocator.alloc(dn * dp);
            auto m15A = allocator.alloc(dn*dm);
            auto m15B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m15A, dA[2][1], dA[2][2]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m15B, dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m15, m15A, m15B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m15B, dm*dp);
            allocator.dealloc(m15A, dn*dm);

            auto m16 = allocator.alloc(dn * dp);
            auto m16A = allocator.alloc(dn*dm);
            auto m16B = allocator.alloc(dm*dp);
            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m16A, dA[0][2], dA[1][1], dA[1][2]);
            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, m16B, dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m16, m16A, m16B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m16B, dm*dp);
            allocator.dealloc(m16A, dn*dm);

            auto m17 = allocator.alloc(dn * dp);
            auto m17A = allocator.alloc(dn*dm);
            auto m17B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m17A, dA[0][2], dA[1][2]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, m17B, dB[1][2], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m17, m17A, m17B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m17B, dm*dp);
            allocator.dealloc(m17A, dn*dm);

            auto m18 = allocator.alloc(dn * dp);
            auto m18A = allocator.alloc(dn*dm);
            auto m18B = allocator.alloc(dm*dp);
            operationEff<Assign, Add, Add>(dn, dm, dm, effA, m18A, dA[1][1], dA[1][2]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m18B, dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m18, m18A, m18B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m18B, dm*dp);
            allocator.dealloc(m18A, dn*dm);

            auto m19 = allocator.alloc(dn * dp);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m19, dA[0][1], dB[1][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);

            auto m20 = allocator.alloc(dn * dp);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m20, dA[1][2], dB[2][1], dn, dm, dp, dp, effA, effB, steps - 1, allocator);

            auto m21 = allocator.alloc(dn * dp);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m21, dA[1][0], dB[0][2], dn, dm, dp, dp, effA, effB, steps - 1, allocator);

            auto m22 = allocator.alloc(dn * dp);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m22, dA[2][0], dB[0][1], dn, dm, dp, dp, effA, effB, steps - 1, allocator);

            auto m23 = allocator.alloc(dn * dp);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3RecursiveLowLevel>(m23, dA[2][2], dB[2][2], dn, dm, dp, dp, effA, effB, steps - 1, allocator);

            operationEff<Assign, Add, Add, Add>(dn, dp, effC, dp, dC[0][0], m6, m14, m19);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[0][1], m1, m4, m5, m6, m12, m14, m15);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[0][2], m6, m7, m9, m10, m14, m16, m18);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[1][0], m2, m3, m4, m6, m14, m16, m17);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[1][1], m2, m4, m5, m6, m20);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[1][2], m14, m16, m17, m18, m21);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][0], m6, m7, m8, m11, m12, m13, m14);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][1], m12, m13, m14, m15, m22);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][2], m6, m7, m8, m9, m23);

            allocator.dealloc(m23, dn*dp);
            allocator.dealloc(m22, dn*dp);
            allocator.dealloc(m21, dn*dp);
            allocator.dealloc(m20, dn*dp);
            allocator.dealloc(m19, dn*dp);
            allocator.dealloc(m18, dn*dp);
            allocator.dealloc(m17, dn*dp);
            allocator.dealloc(m16, dn*dp);
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
    struct GenFastMul3x3RecursiveLowLevelParallel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 3;
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
            auto m16 = allocator.alloc(dn * dp);
            auto m17 = allocator.alloc(dn * dp);
            auto m18 = allocator.alloc(dn * dp);
            auto m19 = allocator.alloc(dn * dp);
            auto m20 = allocator.alloc(dn * dp);
            auto m21 = allocator.alloc(dn * dp);
            auto m22 = allocator.alloc(dn * dp);
            auto m23 = allocator.alloc(dn * dp);
            auto m1A = allocator.alloc(dn*dm);
            auto m2A = allocator.alloc(dn*dm);
            auto m2B = allocator.alloc(dm*dp);
            auto m3B = allocator.alloc(dm*dp);
            auto m4A = allocator.alloc(dn*dm);
            auto m4B = allocator.alloc(dm*dp);
            auto m5A = allocator.alloc(dn*dm);
            auto m5B = allocator.alloc(dm*dp);
            auto m7A = allocator.alloc(dn*dm);
            auto m7B = allocator.alloc(dm*dp);
            auto m8A = allocator.alloc(dn*dm);
            auto m8B = allocator.alloc(dm*dp);
            auto m9A = allocator.alloc(dn*dm);
            auto m9B = allocator.alloc(dm*dp);
            auto m10A = allocator.alloc(dn*dm);
            auto m11B = allocator.alloc(dm*dp);
            auto m12A = allocator.alloc(dn*dm);
            auto m12B = allocator.alloc(dm*dp);
            auto m13A = allocator.alloc(dn*dm);
            auto m13B = allocator.alloc(dm*dp);
            auto m15A = allocator.alloc(dn*dm);
            auto m15B = allocator.alloc(dm*dp);
            auto m16A = allocator.alloc(dn*dm);
            auto m16B = allocator.alloc(dm*dp);
            auto m17A = allocator.alloc(dn*dm);
            auto m17B = allocator.alloc(dm*dp);
            auto m18A = allocator.alloc(dn*dm);
            auto m18B = allocator.alloc(dm*dp);

            ThreadPool pool;

            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, m1A, dA[0][0], dA[0][1], dA[0][2], dA[1][0], dA[1][1], dA[2][1], dA[2][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m1, m1A, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m2A, dA[0][0], dA[1][0]);
                operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m2B, dB[1][0], dB[1][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m2, m2A, m2B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, m3B, dB[0][0], dB[0][1], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m3, dA[1][1], m3B, dn, dm, dp, dp, effA, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m4A, dA[0][0], dA[1][0], dA[1][1]);
                operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, m4B, dB[0][0], dB[0][1], dB[1][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m4, m4A, m4B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m5A, dA[1][0], dA[1][1]);
                operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m5B, dB[0][0], dB[0][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m5, m5A, m5B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m6, dA[0][0], dB[0][0], dn, dm, dp, dp, effA, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m7A, dA[0][0], dA[2][0], dA[2][1]);
                operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, m7B, dB[0][0], dB[0][2], dB[1][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m7, m7A, m7B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add>(dn, dm, dm, effA, m8A, dA[0][0], dA[2][0]);
                operationEff<Assign, Add, Sub>(dn, dm, dm, effB, m8B, dB[0][2], dB[1][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m8, m8A, m8B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m9A, dA[2][0], dA[2][1]);
                operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m9B, dB[0][0], dB[0][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m9, m9A, m9B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, m10A, dA[0][0], dA[0][1], dA[0][2], dA[1][1], dA[1][2], dA[2][0], dA[2][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m10, m10A, dB[1][2], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, m11B, dB[0][0], dB[0][2], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m11, dA[2][1], m11B, dn, dm, dp, dp, effA, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m12A, dA[0][2], dA[2][1], dA[2][2]);
                operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, m12B, dB[1][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m12, m12A, m12B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m13A, dA[0][2], dA[2][2]);
                operationEff<Assign, Add, Sub>(dn, dm, dm, effB, m13B, dB[1][1], dB[2][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m13, m13A, m13B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m14, dA[0][2], dB[2][0], dn, dm, dp, dp, effA, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m15A, dA[2][1], dA[2][2]);
                operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m15B, dB[2][0], dB[2][1]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m15, m15A, m15B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, m16A, dA[0][2], dA[1][1], dA[1][2]);
                operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, m16B, dB[1][2], dB[2][0], dB[2][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m16, m16A, m16B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Sub>(dn, dm, dm, effA, m17A, dA[0][2], dA[1][2]);
                operationEff<Assign, Add, Sub>(dn, dm, dm, effB, m17B, dB[1][2], dB[2][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m17, m17A, m17B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Assign, Add, Add>(dn, dm, dm, effA, m18A, dA[1][1], dA[1][2]);
                operationEff<Assign, Sub, Add>(dn, dm, dm, effB, m18B, dB[2][0], dB[2][2]);
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m18, m18A, m18B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m19, dA[0][1], dB[1][0], dn, dm, dp, dp, effA, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m20, dA[1][2], dB[2][1], dn, dm, dp, dp, effA, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m21, dA[1][0], dB[0][2], dn, dm, dp, dp, effA, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m22, dA[2][0], dB[0][1], dn, dm, dp, dp, effA, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                minSpaceRun<Method, 3, 3, 3, 23, GenFastMul3x3RecursiveMinSpace>(m23, dA[2][2], dB[2][2], dn, dm, dp, dp, effA, effB, steps - 1);
            });

            pool.completeTasksAndStop();

            operationEff<Assign, Add, Add, Add>(dn, dp, effC, dp, dC[0][0], m6, m14, m19);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[0][1], m1, m4, m5, m6, m12, m14, m15);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[0][2], m6, m7, m9, m10, m14, m16, m18);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[1][0], m2, m3, m4, m6, m14, m16, m17);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[1][1], m2, m4, m5, m6, m20);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[1][2], m14, m16, m17, m18, m21);
            operationEff<Assign, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][0], m6, m7, m8, m11, m12, m13, m14);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][1], m12, m13, m14, m15, m22);
            operationEff<Assign, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][2], m6, m7, m8, m9, m23);

            allocator.dealloc(m18B, dm*dp);
            allocator.dealloc(m18A, dn*dm);
            allocator.dealloc(m17B, dm*dp);
            allocator.dealloc(m17A, dn*dm);
            allocator.dealloc(m16B, dm*dp);
            allocator.dealloc(m16A, dn*dm);
            allocator.dealloc(m15B, dm*dp);
            allocator.dealloc(m15A, dn*dm);
            allocator.dealloc(m13B, dm*dp);
            allocator.dealloc(m13A, dn*dm);
            allocator.dealloc(m12B, dm*dp);
            allocator.dealloc(m12A, dn*dm);
            allocator.dealloc(m11B, dm*dp);
            allocator.dealloc(m10A, dn*dm);
            allocator.dealloc(m9B, dm*dp);
            allocator.dealloc(m9A, dn*dm);
            allocator.dealloc(m8B, dm*dp);
            allocator.dealloc(m8A, dn*dm);
            allocator.dealloc(m7B, dm*dp);
            allocator.dealloc(m7A, dn*dm);
            allocator.dealloc(m5B, dm*dp);
            allocator.dealloc(m5A, dn*dm);
            allocator.dealloc(m4B, dm*dp);
            allocator.dealloc(m4A, dn*dm);
            allocator.dealloc(m3B, dm*dp);
            allocator.dealloc(m2B, dm*dp);
            allocator.dealloc(m2A, dn*dm);
            allocator.dealloc(m1A, dn*dm);
            allocator.dealloc(m23, dn*dp);
            allocator.dealloc(m22, dn*dp);
            allocator.dealloc(m21, dn*dp);
            allocator.dealloc(m20, dn*dp);
            allocator.dealloc(m19, dn*dp);
            allocator.dealloc(m18, dn*dp);
            allocator.dealloc(m17, dn*dp);
            allocator.dealloc(m16, dn*dp);
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
    auto genFastMul3x3MinSpace(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 3, 3, 3, 23, detail::GenFastMul3x3RecursiveMinSpace>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto genFastMul3x3LowLevel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevel>, 3, 3, 3, 23, detail::GenFastMul3x3RecursiveLowLevel>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto genFastMul3x3LowLevelParallel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevelParallel>, 3, 3, 3, 23, detail::GenFastMul3x3RecursiveLowLevelParallel>(a, b, steps, 14, 14);
    }
}

#endif

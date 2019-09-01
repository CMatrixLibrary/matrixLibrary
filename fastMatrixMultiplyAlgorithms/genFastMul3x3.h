// this file was generated using fastMatrixMultiplyAlgorithms/generator
#ifndef GEN_FMM_3_3_3_23_H
#define GEN_FMM_3_3_3_23_H
#include "../fmmUtility.h"

namespace fmm::detail {
    struct GenFastMul3x3Recursive {
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
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<Assign>(dn, dp, dp, effC, tempC, dC[0][1]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][0]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<Assign, Assign>(dn, dp, dp, effC, tempC, dC[1][0], dC[1][1]);

            operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][1], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[1][1], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[1][0]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][0], dA[1][1]);
            operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[1][0], dC[1][1]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[1][0], dA[1][1]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[1][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[0][0], dB[0][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<Assign, AddAssign, Assign, AddAssign, AddAssign, Assign, Assign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][1], dC[0][2], dC[1][0], dC[1][1], dC[2][0], dC[2][2]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[2][0], dA[2][1]);
            operationEff<Assign, Add, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][2], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[2][0], dC[2][2]);

            operationEff<Assign, Sub, Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[2][0]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, tempB, dB[0][2], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[2][0], dC[2][2]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[2][0], dA[2][1]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[2][2]);

            operationEff<Assign, Add, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][1], dA[0][2], dA[1][1], dA[1][2], dA[2][0], dA[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, dB[1][2], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[0][2]);

            operationEff<Assign, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effB, tempB, dB[0][0], dB[0][2], dB[1][0], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[2][1], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[2][0]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][2], dA[2][1], dA[2][2]);
            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, Assign>(dn, dp, dp, effC, tempC, dC[0][1], dC[2][0], dC[2][1]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][2], dA[2][2]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][1], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[2][0], dC[2][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[0][2], dB[2][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign, AddAssign, Assign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][1], dC[0][2], dC[1][0], dC[1][2], dC[2][0], dC[2][1]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[2][1], dA[2][2]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][1], dC[2][1]);

            operationEff<Assign, Sub, Add, Add>(dn, dm, dm, effA, tempA, dA[0][2], dA[1][1], dA[1][2]);
            operationEff<Assign, Add, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[1][0], dC[1][2]);

            operationEff<Assign, Add, Sub>(dn, dm, dm, effA, tempA, dA[0][2], dA[1][2]);
            operationEff<Assign, Add, Sub>(dn, dm, dm, effB, tempB, dB[1][2], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[1][0], dC[1][2]);

            operationEff<Assign, Add, Add>(dn, dm, dm, effA, tempA, dA[1][1], dA[1][2]);
            operationEff<Assign, Sub, Add>(dn, dm, dm, effB, tempB, dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][2], dC[1][2]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[0][1], dB[1][0], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[0][0]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[1][2], dB[2][1], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[1][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[1][0], dB[0][2], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[1][2]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[2][0], dB[0][1], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[2][1]);

            nextStep<Method, BaseN, BaseM, BaseP, GenFastMul3x3Recursive>(tempC, dA[2][2], dB[2][2], dn, dm, dp, dp, effA, effB, steps - 1, allocator);
            operationsOnFirstArg<AddAssign>(dn, dp, dp, effC, tempC, dC[2][2]);

            allocator.dealloc(tempC, dn*dp);
            allocator.dealloc(tempB, dm*dp);
            allocator.dealloc(tempA, dn*dm);
        }
    };
}

namespace fmm {
    template<int Method = 0, typename M1, typename M2>
    auto genFastMul3x3MinSpace(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 3, 3, 3, 23, detail::GenFastMul3x3Recursive>(a, b, steps);
    }
}

#endif

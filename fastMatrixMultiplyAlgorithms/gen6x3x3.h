// this file was generated using fastMatrixMultiplyAlgorithms/generator
#ifndef GEN_FMM_6_3_3_40_H
#define GEN_FMM_6_3_3_40_H
#include "../fmmUtility.h"
#include "../ThreadPool.h"

namespace fmm::detail {
    struct Gen6x3x3RecursiveMinSpace {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 6;
            constexpr int BaseM = 3;
            constexpr int BaseP = 3;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto tempA = allocator.alloc(dn * dm);
            auto tempB = allocator.alloc(dm * dp);
            auto tempC = allocator.alloc(dn * dp);

            operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Add, Add, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Add, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<Assign, Assign, Assign, Assign, Assign, Assign, Assign, Assign, Assign>(dn, dp, dp, effC, tempC, std::pair(-1.000000, dC[0][2]), std::pair(-0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(-1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Add, Sub, Add, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Sub, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, Assign, Assign, Assign, Assign, Assign, SubAssign, AddAssign, Assign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(-0.125000, dC[1][2]), std::pair(-0.125000, dC[2][0]), std::pair(-0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(-0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Add, Sub, Add, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Sub, Add, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<Assign, Assign, AddAssign, SubAssign, AddAssign, AddAssign, Assign, SubAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(-1.000000, dC[0][0]), std::pair(-1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Sub, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Add, Sub, Sub, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Add, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Sub, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, SubAssign, AddAssign, SubAssign, AddAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Add, Sub, Sub, Sub, Add, Add, Sub, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Sub, Add, Add, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, SubAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Sub, Add, Sub, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Sub, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Add, Add, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Sub, Add, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Sub, Sub, Add, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Add, Add, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, AddAssign, SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Sub, Sub, Add, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Sub, Add, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, SubAssign, AddAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Add, Sub, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Add, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Add, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Add, Sub, Sub, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, SubAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Sub, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Sub, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, AddAssign, SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Sub, Sub, Add, Sub, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Sub, Add>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, SubAssign, AddAssign, AddAssign, SubAssign, AddAssign, SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Add, Add, Add, Sub, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Add>(dm, dp, dp, effB, tempB, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Sub, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, SubAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Sub, Sub, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Sub, Add>(dm, dp, dp, effB, tempB, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, SubAssign, AddAssign, AddAssign, SubAssign, AddAssign, SubAssign, SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, tempB, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, SubAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Sub, Add, Sub, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, SubAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Sub, Add, Add, Add, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, AddAssign, SubAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Sub, Sub, Sub, Add, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Sub, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Add, Sub, Sub, Sub, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Add, Sub, Add, Add, Add, Add>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, SubAssign, SubAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Sub, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Sub, Add>(dm, dp, dp, effB, tempB, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, SubAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Add, Add, Sub, Add, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Sub, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][2]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][2]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Sub, Sub, Sub, Sub, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Add, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, tempB, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, AddAssign, SubAssign, AddAssign, SubAssign, AddAssign, SubAssign, AddAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, SubAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, AddAssign, SubAssign, SubAssign, SubAssign, AddAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Add, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, AddAssign, SubAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, SubAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            operationEffWithCoeffs<Assign, Add, Sub, Sub, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<AddAssign, SubAssign, AddAssign, AddAssign, AddAssign, AddAssign, SubAssign, AddAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][1]), std::pair(0.125000, dC[1][2]), std::pair(0.125000, dC[2][1]), std::pair(0.125000, dC[2][2]), std::pair(0.125000, dC[3][1]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][1]));

            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, tempA, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveMinSpace>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            operationsOnFirstArgWithCoeffs<SubAssign, AddAssign, AddAssign, AddAssign, SubAssign, SubAssign, SubAssign, SubAssign, AddAssign, AddAssign>(dn, dp, dp, effC, tempC, std::pair(1.000000, dC[0][0]), std::pair(1.000000, dC[0][1]), std::pair(0.125000, dC[1][0]), std::pair(0.125000, dC[2][0]), std::pair(0.125000, dC[3][0]), std::pair(0.125000, dC[3][2]), std::pair(1.000000, dC[4][0]), std::pair(1.000000, dC[4][1]), std::pair(0.125000, dC[5][0]), std::pair(0.125000, dC[5][2]));

            allocator.dealloc(tempC, dn*dp);
            allocator.dealloc(tempB, dm*dp);
            allocator.dealloc(tempA, dn*dm);
        }
    };
}

namespace fmm::detail {
    struct Gen6x3x3RecursiveLowLevel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 6;
            constexpr int BaseM = 3;
            constexpr int BaseP = 3;

            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto m1 = allocator.alloc(dn * dp);
            auto m1A = allocator.alloc(dn*dm);
            auto m1B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Add, Add, Sub, Add>(dn, dm, dm, effA, m1A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Add, Sub, Add>(dm, dp, dp, effB, m1B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m1, m1A, m1B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m1B, dm*dp);
            allocator.dealloc(m1A, dn*dm);

            auto m2 = allocator.alloc(dn * dp);
            auto m2A = allocator.alloc(dn*dm);
            auto m2B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Add, Sub, Add, Sub, Sub>(dn, dm, dm, effA, m2A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Sub, Sub, Add>(dm, dp, dp, effB, m2B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m2, m2A, m2B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m2B, dm*dp);
            allocator.dealloc(m2A, dn*dm);

            auto m3 = allocator.alloc(dn * dp);
            auto m3A = allocator.alloc(dn*dm);
            auto m3B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Sub, Add, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, m3A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Sub, Add, Sub, Add, Add>(dm, dp, dp, effB, m3B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m3, m3A, m3B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m3B, dm*dp);
            allocator.dealloc(m3A, dn*dm);

            auto m4 = allocator.alloc(dn * dp);
            auto m4A = allocator.alloc(dn*dm);
            auto m4B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dm, dm, effA, m4A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Sub, Sub, Add, Add>(dm, dp, dp, effB, m4B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m4, m4A, m4B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m4B, dm*dp);
            allocator.dealloc(m4A, dn*dm);

            auto m5 = allocator.alloc(dn * dp);
            auto m5A = allocator.alloc(dn*dm);
            auto m5B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, m5A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Add, Sub, Sub, Sub>(dm, dp, dp, effB, m5B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m5, m5A, m5B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m5B, dm*dp);
            allocator.dealloc(m5A, dn*dm);

            auto m6 = allocator.alloc(dn * dp);
            auto m6A = allocator.alloc(dn*dm);
            auto m6B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Add, Add, Add>(dn, dm, dm, effA, m6A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Sub, Sub, Add>(dm, dp, dp, effB, m6B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m6, m6A, m6B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m6B, dm*dp);
            allocator.dealloc(m6A, dn*dm);

            auto m7 = allocator.alloc(dn * dp);
            auto m7A = allocator.alloc(dn*dm);
            auto m7B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Sub, Sub, Add, Add, Sub, Add, Sub>(dn, dm, dm, effA, m7A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Sub, Add, Add, Sub, Add>(dm, dp, dp, effB, m7B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m7, m7A, m7B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m7B, dm*dp);
            allocator.dealloc(m7A, dn*dm);

            auto m8 = allocator.alloc(dn * dp);
            auto m8A = allocator.alloc(dn*dm);
            auto m8B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effA, m8A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Sub, Add, Add>(dm, dp, dp, effB, m8B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m8, m8A, m8B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m8B, dm*dp);
            allocator.dealloc(m8A, dn*dm);

            auto m9 = allocator.alloc(dn * dp);
            auto m9A = allocator.alloc(dn*dm);
            auto m9B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Sub, Add, Sub, Add, Add>(dn, dm, dm, effA, m9A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Sub, Add, Sub>(dm, dp, dp, effB, m9B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m9, m9A, m9B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m9B, dm*dp);
            allocator.dealloc(m9A, dn*dm);

            auto m10 = allocator.alloc(dn * dp);
            auto m10A = allocator.alloc(dn*dm);
            auto m10B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dm, dm, effA, m10A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add, Add, Sub>(dm, dp, dp, effB, m10B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m10, m10A, m10B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m10B, dm*dp);
            allocator.dealloc(m10A, dn*dm);

            auto m11 = allocator.alloc(dn * dp);
            auto m11A = allocator.alloc(dn*dm);
            auto m11B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Add, Add, Add, Sub>(dn, dm, dm, effA, m11A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Sub, Add, Add, Add>(dm, dp, dp, effB, m11B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m11, m11A, m11B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m11B, dm*dp);
            allocator.dealloc(m11A, dn*dm);

            auto m12 = allocator.alloc(dn * dp);
            auto m12A = allocator.alloc(dn*dm);
            auto m12B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Sub, Sub, Add, Sub, Add>(dn, dm, dm, effA, m12A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Add, Add, Add, Sub>(dm, dp, dp, effB, m12B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m12, m12A, m12B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m12B, dm*dp);
            allocator.dealloc(m12A, dn*dm);

            auto m13 = allocator.alloc(dn * dp);
            auto m13A = allocator.alloc(dn*dm);
            auto m13B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Sub, Sub, Add, Add, Sub>(dn, dm, dm, effA, m13A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Sub, Add, Sub, Add>(dm, dp, dp, effB, m13B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m13, m13A, m13B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m13B, dm*dp);
            allocator.dealloc(m13A, dn*dm);

            auto m14 = allocator.alloc(dn * dp);
            auto m14A = allocator.alloc(dn*dm);
            auto m14B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, m14A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Add, Sub, Sub>(dm, dp, dp, effB, m14B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m14, m14A, m14B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m14B, dm*dp);
            allocator.dealloc(m14A, dn*dm);

            auto m15 = allocator.alloc(dn * dp);
            auto m15A = allocator.alloc(dn*dm);
            auto m15B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dm, dm, effA, m15A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Add, Add, Add>(dm, dp, dp, effB, m15B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m15, m15A, m15B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m15B, dm*dp);
            allocator.dealloc(m15A, dn*dm);

            auto m16 = allocator.alloc(dn * dp);
            auto m16A = allocator.alloc(dn*dm);
            auto m16B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, m16A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Add, Add, Add>(dm, dp, dp, effB, m16B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m16, m16A, m16B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m16B, dm*dp);
            allocator.dealloc(m16A, dn*dm);

            auto m17 = allocator.alloc(dn * dp);
            auto m17A = allocator.alloc(dn*dm);
            auto m17B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Add, Sub, Sub, Add, Add>(dn, dm, dm, effA, m17A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m17B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m17, m17A, m17B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m17B, dm*dp);
            allocator.dealloc(m17A, dn*dm);

            auto m18 = allocator.alloc(dn * dp);
            auto m18A = allocator.alloc(dn*dm);
            auto m18B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Sub, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, m18A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Sub, Sub>(dm, dp, dp, effB, m18B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m18, m18A, m18B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m18B, dm*dp);
            allocator.dealloc(m18A, dn*dm);

            auto m19 = allocator.alloc(dn * dp);
            auto m19A = allocator.alloc(dn*dm);
            auto m19B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Sub, Sub, Add, Sub, Sub, Sub>(dn, dm, dm, effA, m19A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Sub, Add>(dm, dp, dp, effB, m19B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m19, m19A, m19B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m19B, dm*dp);
            allocator.dealloc(m19A, dn*dm);

            auto m20 = allocator.alloc(dn * dp);
            auto m20A = allocator.alloc(dn*dm);
            auto m20B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Add, Add, Add, Sub, Sub, Sub>(dn, dm, dm, effA, m20A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Add>(dm, dp, dp, effB, m20B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m20, m20A, m20B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m20B, dm*dp);
            allocator.dealloc(m20A, dn*dm);

            auto m21 = allocator.alloc(dn * dp);
            auto m21A = allocator.alloc(dn*dm);
            auto m21B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Sub, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, m21A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m21B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m21, m21A, m21B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m21B, dm*dp);
            allocator.dealloc(m21A, dn*dm);

            auto m22 = allocator.alloc(dn * dp);
            auto m22A = allocator.alloc(dn*dm);
            auto m22B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effA, m22A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m22B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m22, m22A, m22B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m22B, dm*dp);
            allocator.dealloc(m22A, dn*dm);

            auto m23 = allocator.alloc(dn * dp);
            auto m23A = allocator.alloc(dn*dm);
            auto m23B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Sub, Sub, Add, Add>(dn, dm, dm, effA, m23A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Sub, Add>(dm, dp, dp, effB, m23B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m23, m23A, m23B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m23B, dm*dp);
            allocator.dealloc(m23A, dn*dm);

            auto m24 = allocator.alloc(dn * dp);
            auto m24A = allocator.alloc(dn*dm);
            auto m24B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Sub>(dn, dm, dm, effA, m24A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m24B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m24, m24A, m24B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m24B, dm*dp);
            allocator.dealloc(m24A, dn*dm);

            auto m25 = allocator.alloc(dn * dp);
            auto m25A = allocator.alloc(dn*dm);
            auto m25B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Sub, Add, Sub, Add, Add>(dn, dm, dm, effA, m25A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m25B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m25, m25A, m25B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m25B, dm*dp);
            allocator.dealloc(m25A, dn*dm);

            auto m26 = allocator.alloc(dn * dp);
            auto m26A = allocator.alloc(dn*dm);
            auto m26B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Sub, Add, Add, Add, Sub, Sub>(dn, dm, dm, effA, m26A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m26B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m26, m26A, m26B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m26B, dm*dp);
            allocator.dealloc(m26A, dn*dm);

            auto m27 = allocator.alloc(dn * dp);
            auto m27A = allocator.alloc(dn*dm);
            auto m27B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Sub, Sub, Sub, Add, Add, Sub>(dn, dm, dm, effA, m27A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Add, Sub, Sub>(dm, dp, dp, effB, m27B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m27, m27A, m27B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m27B, dm*dp);
            allocator.dealloc(m27A, dn*dm);

            auto m28 = allocator.alloc(dn * dp);
            auto m28A = allocator.alloc(dn*dm);
            auto m28B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Add, Sub, Sub, Sub, Add, Sub>(dn, dm, dm, effA, m28A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m28B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m28, m28A, m28B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m28B, dm*dp);
            allocator.dealloc(m28A, dn*dm);

            auto m29 = allocator.alloc(dn * dp);
            auto m29A = allocator.alloc(dn*dm);
            auto m29B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Add, Sub, Add, Add, Add, Add>(dn, dm, dm, effA, m29A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m29B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m29, m29A, m29B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m29B, dm*dp);
            allocator.dealloc(m29A, dn*dm);

            auto m30 = allocator.alloc(dn * dp);
            auto m30A = allocator.alloc(dn*dm);
            auto m30B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Sub, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, m30A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Sub, Add>(dm, dp, dp, effB, m30B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m30, m30A, m30B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m30B, dm*dp);
            allocator.dealloc(m30A, dn*dm);

            auto m31 = allocator.alloc(dn * dp);
            auto m31A = allocator.alloc(dn*dm);
            auto m31B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m31A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m31B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m31, m31A, m31B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m31B, dm*dp);
            allocator.dealloc(m31A, dn*dm);

            auto m32 = allocator.alloc(dn * dp);
            auto m32A = allocator.alloc(dn*dm);
            auto m32B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Add, Add, Sub, Add, Sub, Sub>(dn, dm, dm, effA, m32A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Sub, Add>(dm, dp, dp, effB, m32B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m32, m32A, m32B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m32B, dm*dp);
            allocator.dealloc(m32A, dn*dm);

            auto m33 = allocator.alloc(dn * dp);
            auto m33A = allocator.alloc(dn*dm);
            auto m33B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Sub, Sub, Sub, Sub, Add, Sub>(dn, dm, dm, effA, m33A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m33B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m33, m33A, m33B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m33B, dm*dp);
            allocator.dealloc(m33A, dn*dm);

            auto m34 = allocator.alloc(dn * dp);
            auto m34A = allocator.alloc(dn*dm);
            auto m34B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dm, dm, effA, m34A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m34B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m34, m34A, m34B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m34B, dm*dp);
            allocator.dealloc(m34A, dn*dm);

            auto m35 = allocator.alloc(dn * dp);
            auto m35A = allocator.alloc(dn*dm);
            auto m35B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Add, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, m35A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m35B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m35, m35A, m35B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m35B, dm*dp);
            allocator.dealloc(m35A, dn*dm);

            auto m36 = allocator.alloc(dn * dp);
            auto m36A = allocator.alloc(dn*dm);
            auto m36B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, m36A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Add, Add, Add, Add>(dm, dp, dp, effB, m36B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m36, m36A, m36B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m36B, dm*dp);
            allocator.dealloc(m36A, dn*dm);

            auto m37 = allocator.alloc(dn * dp);
            auto m37A = allocator.alloc(dn*dm);
            auto m37B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m37A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m37B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m37, m37A, m37B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m37B, dm*dp);
            allocator.dealloc(m37A, dn*dm);

            auto m38 = allocator.alloc(dn * dp);
            auto m38A = allocator.alloc(dn*dm);
            auto m38B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Add, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, m38A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m38B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m38, m38A, m38B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m38B, dm*dp);
            allocator.dealloc(m38A, dn*dm);

            auto m39 = allocator.alloc(dn * dp);
            auto m39A = allocator.alloc(dn*dm);
            auto m39B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Sub, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m39A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m39B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m39, m39A, m39B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m39B, dm*dp);
            allocator.dealloc(m39A, dn*dm);

            auto m40 = allocator.alloc(dn * dp);
            auto m40A = allocator.alloc(dn*dm);
            auto m40B = allocator.alloc(dm*dp);
            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m40A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
            operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m40B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
            nextStep<Method, BaseN, BaseM, BaseP, Gen6x3x3RecursiveLowLevel>(m40, m40A, m40B, dn, dm, dp, dp, dm, dm, steps - 1, allocator);
            allocator.dealloc(m40B, dm*dp);
            allocator.dealloc(m40A, dn*dm);

            operationEff<Assign, Sub, Sub, Sub, Sub, Sub, Sub, Add, Add, Sub, Sub, Sub, Add, Add, Add, Add, Sub, Sub, Add, Sub, Sub, Sub, Add, Add, Sub>(dn, dp, effC, dp, dC[0][0], m3, m6, m7, m9, m10, m11, m14, m15, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Sub, Add, Add, Add, Add, Sub, Add, Add, Sub, Sub, Add, Add, Sub, Add, Sub, Add, Add, Add, Sub, Sub, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[0][1], m3, m6, m7, m9, m10, m11, m14, m15, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Sub, Sub, Sub, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Add, Sub, Sub, Add, Add>(dn, dp, effC, dp, dC[0][2], m1, m2, m4, m5, m8, m12, m13, m16, m23, m24, m27, m28, m29, m30, m31, m32);
            operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Add, Sub, Add, Sub, Sub, Add, Sub, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[1][0], std::pair(0.125000, m1), std::pair(0.125000, m4), std::pair(0.125000, m6), std::pair(0.125000, m8), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m13), std::pair(0.125000, m15), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Add, Add, Sub, Sub, Add, Add, Add, Add, Add, Add, Sub, Sub, Sub, Add, Sub, Sub, Sub, Sub, Sub, Sub, Sub, Add, Sub, Add>(dn, dp, effC, dp, dC[1][1], std::pair(0.125000, m2), std::pair(0.125000, m3), std::pair(0.125000, m5), std::pair(0.125000, m7), std::pair(0.125000, m9), std::pair(0.125000, m12), std::pair(0.125000, m14), std::pair(0.125000, m16), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Sub, Add, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub, Sub, Add, Sub, Add, Sub, Sub, Sub, Add>(dn, dp, effC, dp, dC[1][2], std::pair(0.125000, m2), std::pair(0.125000, m3), std::pair(0.125000, m5), std::pair(0.125000, m7), std::pair(0.125000, m9), std::pair(0.125000, m12), std::pair(0.125000, m14), std::pair(0.125000, m16), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Sub, Add, Sub, Add, Sub, Sub, Add, Add, Add, Add, Add, Sub, Add, Sub, Add, Add>(dn, dp, effC, dp, dC[2][0], std::pair(0.125000, m2), std::pair(0.125000, m3), std::pair(0.125000, m5), std::pair(0.125000, m7), std::pair(0.125000, m9), std::pair(0.125000, m12), std::pair(0.125000, m14), std::pair(0.125000, m16), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Add, Add, Add, Sub, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][1], std::pair(0.125000, m1), std::pair(0.125000, m4), std::pair(0.125000, m6), std::pair(0.125000, m8), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m13), std::pair(0.125000, m15), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Add, Sub, Sub, Sub, Sub, Sub, Sub, Add, Sub, Add, Sub, Add, Add, Sub, Add, Sub, Sub, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dp, effC, dp, dC[2][2], std::pair(0.125000, m1), std::pair(0.125000, m4), std::pair(0.125000, m6), std::pair(0.125000, m8), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m13), std::pair(0.125000, m15), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Sub, Add, Sub, Sub, Add, Sub, Sub, Add, Add, Add, Add, Add, Add, Add, Add, Add, Sub, Sub, Sub>(dn, dp, effC, dp, dC[3][0], std::pair(0.125000, m2), std::pair(0.125000, m5), std::pair(0.125000, m6), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m12), std::pair(0.125000, m15), std::pair(0.125000, m16), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Add, Add, Add, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dp, effC, dp, dC[3][1], std::pair(0.125000, m1), std::pair(0.125000, m3), std::pair(0.125000, m4), std::pair(0.125000, m7), std::pair(0.125000, m8), std::pair(0.125000, m9), std::pair(0.125000, m13), std::pair(0.125000, m14), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Add, Add, Sub, Sub, Add, Sub, Sub, Add, Add, Sub, Add, Sub, Sub, Sub, Add, Add, Sub, Add, Add, Sub, Sub, Sub, Sub, Sub>(dn, dp, effC, dp, dC[3][2], std::pair(0.125000, m2), std::pair(0.125000, m5), std::pair(0.125000, m6), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m12), std::pair(0.125000, m15), std::pair(0.125000, m16), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEff<Assign, Add, Sub, Sub, Add, Sub, Add, Sub, Add, Add, Sub, Add, Add, Add, Add, Add, Sub, Sub, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dp, effC, dp, dC[4][0], m1, m2, m4, m5, m8, m12, m13, m16, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Sub, Add, Sub, Add, Sub, Sub, Add, Add, Add, Sub, Sub, Sub, Add, Sub, Sub, Sub, Add, Sub, Add, Sub, Sub, Sub, Add, Sub>(dn, dp, effC, dp, dC[4][1], m1, m2, m4, m5, m8, m12, m13, m16, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Add, Sub, Add, Sub, Add, Add, Add, Add, Sub, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dp, effC, dp, dC[4][2], m3, m6, m7, m9, m10, m11, m14, m15, m23, m24, m27, m28, m29, m30, m31, m32);
            operationEffWithCoeffs<Assign, Add, Sub, Sub, Add, Add, Sub, Add, Sub, Sub, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dp, effC, dp, dC[5][0], std::pair(0.125000, m1), std::pair(0.125000, m3), std::pair(0.125000, m4), std::pair(0.125000, m7), std::pair(0.125000, m8), std::pair(0.125000, m9), std::pair(0.125000, m13), std::pair(0.125000, m14), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Sub, Sub, Sub, Add, Add, Sub, Add, Add, Add, Sub, Sub, Add>(dn, dp, effC, dp, dC[5][1], std::pair(0.125000, m2), std::pair(0.125000, m5), std::pair(0.125000, m6), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m12), std::pair(0.125000, m15), std::pair(0.125000, m16), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Add, Add, Sub, Add, Sub, Add, Sub, Add, Add, Add, Add, Add, Sub, Sub, Sub, Sub, Add, Add, Add, Add>(dn, dp, effC, dp, dC[5][2], std::pair(0.125000, m1), std::pair(0.125000, m3), std::pair(0.125000, m4), std::pair(0.125000, m7), std::pair(0.125000, m8), std::pair(0.125000, m9), std::pair(0.125000, m13), std::pair(0.125000, m14), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));

            allocator.dealloc(m40, dn*dp);
            allocator.dealloc(m39, dn*dp);
            allocator.dealloc(m38, dn*dp);
            allocator.dealloc(m37, dn*dp);
            allocator.dealloc(m36, dn*dp);
            allocator.dealloc(m35, dn*dp);
            allocator.dealloc(m34, dn*dp);
            allocator.dealloc(m33, dn*dp);
            allocator.dealloc(m32, dn*dp);
            allocator.dealloc(m31, dn*dp);
            allocator.dealloc(m30, dn*dp);
            allocator.dealloc(m29, dn*dp);
            allocator.dealloc(m28, dn*dp);
            allocator.dealloc(m27, dn*dp);
            allocator.dealloc(m26, dn*dp);
            allocator.dealloc(m25, dn*dp);
            allocator.dealloc(m24, dn*dp);
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
    struct Gen6x3x3RecursiveLowLevelParallel {
        template<int Method, typename T>
        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 6;
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
            auto m24 = allocator.alloc(dn * dp);
            auto m25 = allocator.alloc(dn * dp);
            auto m26 = allocator.alloc(dn * dp);
            auto m27 = allocator.alloc(dn * dp);
            auto m28 = allocator.alloc(dn * dp);
            auto m29 = allocator.alloc(dn * dp);
            auto m30 = allocator.alloc(dn * dp);
            auto m31 = allocator.alloc(dn * dp);
            auto m32 = allocator.alloc(dn * dp);
            auto m33 = allocator.alloc(dn * dp);
            auto m34 = allocator.alloc(dn * dp);
            auto m35 = allocator.alloc(dn * dp);
            auto m36 = allocator.alloc(dn * dp);
            auto m37 = allocator.alloc(dn * dp);
            auto m38 = allocator.alloc(dn * dp);
            auto m39 = allocator.alloc(dn * dp);
            auto m40 = allocator.alloc(dn * dp);
            auto m1A = allocator.alloc(dn*dm);
            auto m1B = allocator.alloc(dm*dp);
            auto m2A = allocator.alloc(dn*dm);
            auto m2B = allocator.alloc(dm*dp);
            auto m3A = allocator.alloc(dn*dm);
            auto m3B = allocator.alloc(dm*dp);
            auto m4A = allocator.alloc(dn*dm);
            auto m4B = allocator.alloc(dm*dp);
            auto m5A = allocator.alloc(dn*dm);
            auto m5B = allocator.alloc(dm*dp);
            auto m6A = allocator.alloc(dn*dm);
            auto m6B = allocator.alloc(dm*dp);
            auto m7A = allocator.alloc(dn*dm);
            auto m7B = allocator.alloc(dm*dp);
            auto m8A = allocator.alloc(dn*dm);
            auto m8B = allocator.alloc(dm*dp);
            auto m9A = allocator.alloc(dn*dm);
            auto m9B = allocator.alloc(dm*dp);
            auto m10A = allocator.alloc(dn*dm);
            auto m10B = allocator.alloc(dm*dp);
            auto m11A = allocator.alloc(dn*dm);
            auto m11B = allocator.alloc(dm*dp);
            auto m12A = allocator.alloc(dn*dm);
            auto m12B = allocator.alloc(dm*dp);
            auto m13A = allocator.alloc(dn*dm);
            auto m13B = allocator.alloc(dm*dp);
            auto m14A = allocator.alloc(dn*dm);
            auto m14B = allocator.alloc(dm*dp);
            auto m15A = allocator.alloc(dn*dm);
            auto m15B = allocator.alloc(dm*dp);
            auto m16A = allocator.alloc(dn*dm);
            auto m16B = allocator.alloc(dm*dp);
            auto m17A = allocator.alloc(dn*dm);
            auto m17B = allocator.alloc(dm*dp);
            auto m18A = allocator.alloc(dn*dm);
            auto m18B = allocator.alloc(dm*dp);
            auto m19A = allocator.alloc(dn*dm);
            auto m19B = allocator.alloc(dm*dp);
            auto m20A = allocator.alloc(dn*dm);
            auto m20B = allocator.alloc(dm*dp);
            auto m21A = allocator.alloc(dn*dm);
            auto m21B = allocator.alloc(dm*dp);
            auto m22A = allocator.alloc(dn*dm);
            auto m22B = allocator.alloc(dm*dp);
            auto m23A = allocator.alloc(dn*dm);
            auto m23B = allocator.alloc(dm*dp);
            auto m24A = allocator.alloc(dn*dm);
            auto m24B = allocator.alloc(dm*dp);
            auto m25A = allocator.alloc(dn*dm);
            auto m25B = allocator.alloc(dm*dp);
            auto m26A = allocator.alloc(dn*dm);
            auto m26B = allocator.alloc(dm*dp);
            auto m27A = allocator.alloc(dn*dm);
            auto m27B = allocator.alloc(dm*dp);
            auto m28A = allocator.alloc(dn*dm);
            auto m28B = allocator.alloc(dm*dp);
            auto m29A = allocator.alloc(dn*dm);
            auto m29B = allocator.alloc(dm*dp);
            auto m30A = allocator.alloc(dn*dm);
            auto m30B = allocator.alloc(dm*dp);
            auto m31A = allocator.alloc(dn*dm);
            auto m31B = allocator.alloc(dm*dp);
            auto m32A = allocator.alloc(dn*dm);
            auto m32B = allocator.alloc(dm*dp);
            auto m33A = allocator.alloc(dn*dm);
            auto m33B = allocator.alloc(dm*dp);
            auto m34A = allocator.alloc(dn*dm);
            auto m34B = allocator.alloc(dm*dp);
            auto m35A = allocator.alloc(dn*dm);
            auto m35B = allocator.alloc(dm*dp);
            auto m36A = allocator.alloc(dn*dm);
            auto m36B = allocator.alloc(dm*dp);
            auto m37A = allocator.alloc(dn*dm);
            auto m37B = allocator.alloc(dm*dp);
            auto m38A = allocator.alloc(dn*dm);
            auto m38B = allocator.alloc(dm*dp);
            auto m39A = allocator.alloc(dn*dm);
            auto m39B = allocator.alloc(dm*dp);
            auto m40A = allocator.alloc(dn*dm);
            auto m40B = allocator.alloc(dm*dp);

            ThreadPool pool;

            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Add, Add, Sub, Add>(dn, dm, dm, effA, m1A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
                operationEff<Assign, Add, Add, Add, Add, Sub, Add>(dm, dp, dp, effB, m1B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m1, m1A, m1B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Add, Sub, Add, Sub, Sub>(dn, dm, dm, effA, m2A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Sub, Add, Sub, Sub, Add>(dm, dp, dp, effB, m2B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m2, m2A, m2B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Sub, Add, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, m3A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Sub, Add, Sub, Add, Add>(dm, dp, dp, effB, m3B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m3, m3A, m3B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dm, dm, effA, m4A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Add, Sub, Sub, Add, Add>(dm, dp, dp, effB, m4B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m4, m4A, m4B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, m5A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Add, Add, Sub, Sub, Sub>(dm, dp, dp, effB, m5B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m5, m5A, m5B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Add, Add, Add>(dn, dm, dm, effA, m6A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Add, Add, Sub, Sub, Add>(dm, dp, dp, effB, m6B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m6, m6A, m6B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Sub, Sub, Add, Add, Sub, Add, Sub>(dn, dm, dm, effA, m7A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Sub, Add, Add, Sub, Add>(dm, dp, dp, effB, m7B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m7, m7A, m7B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effA, m8A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
                operationEff<Assign, Add, Add, Add, Sub, Add, Add>(dm, dp, dp, effB, m8B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m8, m8A, m8B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Sub, Add, Sub, Add, Add>(dn, dm, dm, effA, m9A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Add, Add, Sub, Add, Sub>(dm, dp, dp, effB, m9B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m9, m9A, m9B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dm, dm, effA, m10A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Add, Add, Add, Sub>(dm, dp, dp, effB, m10B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m10, m10A, m10B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Add, Add, Add, Sub>(dn, dm, dm, effA, m11A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Add, Sub, Add, Add, Add>(dm, dp, dp, effB, m11B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m11, m11A, m11B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Sub, Sub, Add, Sub, Add>(dn, dm, dm, effA, m12A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Add, Add, Add, Add, Sub>(dm, dp, dp, effB, m12B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m12, m12A, m12B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Sub, Sub, Add, Add, Sub>(dn, dm, dm, effA, m13A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Add, Sub, Add, Sub, Add>(dm, dp, dp, effB, m13B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m13, m13A, m13B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, m14A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Add, Add, Add, Sub, Sub>(dm, dp, dp, effB, m14B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m14, m14A, m14B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dm, dm, effA, m15A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][1]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Add, Add, Add, Add, Add>(dm, dp, dp, effB, m15B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m15, m15A, m15B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, m16A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Sub, Add, Add, Add, Add>(dm, dp, dp, effB, m16B, dB[0][0], dB[0][2], dB[1][1], dB[1][2], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m16, m16A, m16B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Add, Sub, Sub, Add, Add>(dn, dm, dm, effA, m17A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m17B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m17, m17A, m17B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Sub, Add, Add, Sub, Sub, Add>(dn, dm, dm, effA, m18A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Add, Add, Sub, Sub>(dm, dp, dp, effB, m18B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m18, m18A, m18B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Sub, Sub, Add, Sub, Sub, Sub>(dn, dm, dm, effA, m19A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Add, Sub, Add>(dm, dp, dp, effB, m19B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m19, m19A, m19B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Add, Sub, Add, Add, Add, Sub, Sub, Sub>(dn, dm, dm, effA, m20A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Add, Add, Add, Add>(dm, dp, dp, effB, m20B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m20, m20A, m20B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Sub, Add, Sub, Add, Add, Add>(dn, dm, dm, effA, m21A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m21B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m21, m21A, m21B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dm, dm, effA, m22A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m22B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m22, m22A, m22B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Sub, Sub, Add, Add>(dn, dm, dm, effA, m23A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Add, Sub, Add>(dm, dp, dp, effB, m23B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m23, m23A, m23B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Sub>(dn, dm, dm, effA, m24A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m24B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m24, m24A, m24B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Sub, Add, Sub, Add, Add>(dn, dm, dm, effA, m25A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m25B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m25, m25A, m25B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Sub, Add, Add, Add, Sub, Sub>(dn, dm, dm, effA, m26A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m26B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m26, m26A, m26B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Sub, Sub, Sub, Add, Add, Sub>(dn, dm, dm, effA, m27A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Add, Sub, Sub>(dm, dp, dp, effB, m27B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m27, m27A, m27B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Add, Sub, Sub, Sub, Add, Sub>(dn, dm, dm, effA, m28A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Sub, Sub>(dm, dp, dp, effB, m28B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m28, m28A, m28B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Add, Add, Add, Sub, Add, Add, Add, Add>(dn, dm, dm, effA, m29A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m29B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m29, m29A, m29B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Add, Add, Sub, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, m30A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Sub, Sub, Add>(dm, dp, dp, effB, m30B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m30, m30A, m30B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m31A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m31B, dB[1][1], dB[1][2], dB[2][1], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m31, m31A, m31B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Add, Add, Sub, Add, Sub, Sub>(dn, dm, dm, effA, m32A, std::pair(0.125, dA[0][2]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][2]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Sub, Sub, Add>(dm, dp, dp, effB, m32B, dB[0][0], dB[0][2], dB[2][0], dB[2][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m32, m32A, m32B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Sub, Sub, Sub, Sub, Add, Sub>(dn, dm, dm, effA, m33A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m33B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m33, m33A, m33B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Add, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dm, dm, effA, m34A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m34B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m34, m34A, m34B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Sub, Add, Add, Add, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, m35A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Add, Sub, Add, Sub>(dm, dp, dp, effB, m35B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m35, m35A, m35B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Sub, Sub>(dn, dm, dm, effA, m36A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Add, Add, Add, Add>(dm, dp, dp, effB, m36B, dB[0][0], dB[0][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m36, m36A, m36B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m37A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Add, Add, Sub>(dm, dp, dp, effB, m37B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m37, m37A, m37B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Add, Sub, Add, Sub, Add, Add, Add, Sub>(dn, dm, dm, effA, m38A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m38B, dB[1][0], dB[1][1], dB[2][0], dB[2][1]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m38, m38A, m38B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Sub, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m39A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][0]), std::pair(1, dA[1][2]), std::pair(1, dA[2][0]), std::pair(1, dA[2][2]), std::pair(1, dA[3][0]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][0]));
                operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m39B, dB[0][1], dB[0][2], dB[1][1], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m39, m39A, m39B, dn, dm, dp, dp, dm, dm, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dm, dm, effA, m40A, std::pair(0.125, dA[0][0]), std::pair(0.125, dA[0][1]), std::pair(1, dA[1][1]), std::pair(1, dA[2][1]), std::pair(1, dA[3][1]), std::pair(1, dA[3][2]), std::pair(0.125, dA[4][0]), std::pair(0.125, dA[4][1]), std::pair(1, dA[5][1]), std::pair(1, dA[5][2]));
                operationEff<Assign, Sub, Sub, Add, Add>(dm, dp, dp, effB, m40B, dB[0][0], dB[0][2], dB[1][0], dB[1][2]);
                minSpaceRun<Method, 6, 3, 3, 40, Gen6x3x3RecursiveMinSpace>(m40, m40A, m40B, dn, dm, dp, dp, dm, dm, steps - 1);
            });

            pool.completeTasksAndStop();

            operationEff<Assign, Sub, Sub, Sub, Sub, Sub, Sub, Add, Add, Sub, Sub, Sub, Add, Add, Add, Add, Sub, Sub, Add, Sub, Sub, Sub, Add, Add, Sub>(dn, dp, effC, dp, dC[0][0], m3, m6, m7, m9, m10, m11, m14, m15, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Sub, Add, Add, Add, Add, Sub, Add, Add, Sub, Sub, Add, Add, Sub, Add, Sub, Add, Add, Add, Sub, Sub, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[0][1], m3, m6, m7, m9, m10, m11, m14, m15, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Sub, Sub, Sub, Sub, Add, Sub, Sub, Add, Add, Sub, Sub, Add, Sub, Sub, Add, Add>(dn, dp, effC, dp, dC[0][2], m1, m2, m4, m5, m8, m12, m13, m16, m23, m24, m27, m28, m29, m30, m31, m32);
            operationEffWithCoeffs<Assign, Sub, Sub, Sub, Add, Add, Add, Sub, Add, Sub, Sub, Add, Sub, Add, Add, Sub, Add>(dn, dp, effC, dp, dC[1][0], std::pair(0.125000, m1), std::pair(0.125000, m4), std::pair(0.125000, m6), std::pair(0.125000, m8), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m13), std::pair(0.125000, m15), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Add, Add, Sub, Sub, Add, Add, Add, Add, Add, Add, Sub, Sub, Sub, Add, Sub, Sub, Sub, Sub, Sub, Sub, Sub, Add, Sub, Add>(dn, dp, effC, dp, dC[1][1], std::pair(0.125000, m2), std::pair(0.125000, m3), std::pair(0.125000, m5), std::pair(0.125000, m7), std::pair(0.125000, m9), std::pair(0.125000, m12), std::pair(0.125000, m14), std::pair(0.125000, m16), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Sub, Add, Add, Add, Sub, Add, Add, Sub, Add, Sub, Add, Sub, Sub, Add, Sub, Add, Sub, Sub, Sub, Add>(dn, dp, effC, dp, dC[1][2], std::pair(0.125000, m2), std::pair(0.125000, m3), std::pair(0.125000, m5), std::pair(0.125000, m7), std::pair(0.125000, m9), std::pair(0.125000, m12), std::pair(0.125000, m14), std::pair(0.125000, m16), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Sub, Add, Sub, Add, Sub, Sub, Add, Add, Add, Add, Add, Sub, Add, Sub, Add, Add>(dn, dp, effC, dp, dC[2][0], std::pair(0.125000, m2), std::pair(0.125000, m3), std::pair(0.125000, m5), std::pair(0.125000, m7), std::pair(0.125000, m9), std::pair(0.125000, m12), std::pair(0.125000, m14), std::pair(0.125000, m16), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Sub, Add, Add, Add, Sub, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add>(dn, dp, effC, dp, dC[2][1], std::pair(0.125000, m1), std::pair(0.125000, m4), std::pair(0.125000, m6), std::pair(0.125000, m8), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m13), std::pair(0.125000, m15), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Add, Sub, Sub, Sub, Sub, Sub, Sub, Add, Sub, Add, Sub, Add, Add, Sub, Add, Sub, Sub, Sub, Add, Add, Sub, Sub, Sub, Add>(dn, dp, effC, dp, dC[2][2], std::pair(0.125000, m1), std::pair(0.125000, m4), std::pair(0.125000, m6), std::pair(0.125000, m8), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m13), std::pair(0.125000, m15), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Sub, Add, Add, Sub, Add, Sub, Add, Sub, Sub, Add, Sub, Sub, Add, Add, Add, Add, Add, Add, Add, Add, Add, Sub, Sub, Sub>(dn, dp, effC, dp, dC[3][0], std::pair(0.125000, m2), std::pair(0.125000, m5), std::pair(0.125000, m6), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m12), std::pair(0.125000, m15), std::pair(0.125000, m16), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Add, Add, Add, Add, Sub, Sub, Add, Add, Add, Add, Add, Sub, Add, Add, Sub, Sub>(dn, dp, effC, dp, dC[3][1], std::pair(0.125000, m1), std::pair(0.125000, m3), std::pair(0.125000, m4), std::pair(0.125000, m7), std::pair(0.125000, m8), std::pair(0.125000, m9), std::pair(0.125000, m13), std::pair(0.125000, m14), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Add, Add, Sub, Sub, Add, Sub, Sub, Add, Add, Sub, Add, Sub, Sub, Sub, Add, Add, Sub, Add, Add, Sub, Sub, Sub, Sub, Sub>(dn, dp, effC, dp, dC[3][2], std::pair(0.125000, m2), std::pair(0.125000, m5), std::pair(0.125000, m6), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m12), std::pair(0.125000, m15), std::pair(0.125000, m16), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEff<Assign, Add, Sub, Sub, Add, Sub, Add, Sub, Add, Add, Sub, Add, Add, Add, Add, Add, Sub, Sub, Add, Sub, Sub, Add, Sub, Add, Sub>(dn, dp, effC, dp, dC[4][0], m1, m2, m4, m5, m8, m12, m13, m16, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Sub, Add, Sub, Add, Sub, Sub, Add, Add, Add, Sub, Sub, Sub, Add, Sub, Sub, Sub, Add, Sub, Add, Sub, Sub, Sub, Add, Sub>(dn, dp, effC, dp, dC[4][1], m1, m2, m4, m5, m8, m12, m13, m16, m17, m18, m19, m20, m21, m22, m25, m26, m33, m34, m35, m36, m37, m38, m39, m40);
            operationEff<Assign, Add, Sub, Add, Sub, Add, Add, Add, Add, Sub, Add, Sub, Add, Sub, Add, Sub, Add>(dn, dp, effC, dp, dC[4][2], m3, m6, m7, m9, m10, m11, m14, m15, m23, m24, m27, m28, m29, m30, m31, m32);
            operationEffWithCoeffs<Assign, Add, Sub, Sub, Add, Add, Sub, Add, Sub, Sub, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Add, Sub, Sub, Sub, Add>(dn, dp, effC, dp, dC[5][0], std::pair(0.125000, m1), std::pair(0.125000, m3), std::pair(0.125000, m4), std::pair(0.125000, m7), std::pair(0.125000, m8), std::pair(0.125000, m9), std::pair(0.125000, m13), std::pair(0.125000, m14), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));
            operationEffWithCoeffs<Assign, Sub, Sub, Add, Sub, Sub, Sub, Sub, Add, Add, Sub, Add, Add, Add, Sub, Sub, Add>(dn, dp, effC, dp, dC[5][1], std::pair(0.125000, m2), std::pair(0.125000, m5), std::pair(0.125000, m6), std::pair(0.125000, m10), std::pair(0.125000, m11), std::pair(0.125000, m12), std::pair(0.125000, m15), std::pair(0.125000, m16), std::pair(0.125000, m18), std::pair(0.125000, m20), std::pair(0.125000, m25), std::pair(0.125000, m26), std::pair(0.125000, m33), std::pair(0.125000, m35), std::pair(0.125000, m36), std::pair(0.125000, m39));
            operationEffWithCoeffs<Assign, Add, Sub, Add, Add, Add, Add, Sub, Add, Sub, Add, Sub, Add, Add, Add, Add, Add, Sub, Sub, Sub, Sub, Add, Add, Add, Add>(dn, dp, effC, dp, dC[5][2], std::pair(0.125000, m1), std::pair(0.125000, m3), std::pair(0.125000, m4), std::pair(0.125000, m7), std::pair(0.125000, m8), std::pair(0.125000, m9), std::pair(0.125000, m13), std::pair(0.125000, m14), std::pair(0.125000, m17), std::pair(0.125000, m19), std::pair(0.125000, m21), std::pair(0.125000, m22), std::pair(0.125000, m23), std::pair(0.125000, m24), std::pair(0.125000, m27), std::pair(0.125000, m28), std::pair(0.125000, m29), std::pair(0.125000, m30), std::pair(0.125000, m31), std::pair(0.125000, m32), std::pair(0.125000, m34), std::pair(0.125000, m37), std::pair(0.125000, m38), std::pair(0.125000, m40));

            allocator.dealloc(m40B, dm*dp);
            allocator.dealloc(m40A, dn*dm);
            allocator.dealloc(m39B, dm*dp);
            allocator.dealloc(m39A, dn*dm);
            allocator.dealloc(m38B, dm*dp);
            allocator.dealloc(m38A, dn*dm);
            allocator.dealloc(m37B, dm*dp);
            allocator.dealloc(m37A, dn*dm);
            allocator.dealloc(m36B, dm*dp);
            allocator.dealloc(m36A, dn*dm);
            allocator.dealloc(m35B, dm*dp);
            allocator.dealloc(m35A, dn*dm);
            allocator.dealloc(m34B, dm*dp);
            allocator.dealloc(m34A, dn*dm);
            allocator.dealloc(m33B, dm*dp);
            allocator.dealloc(m33A, dn*dm);
            allocator.dealloc(m32B, dm*dp);
            allocator.dealloc(m32A, dn*dm);
            allocator.dealloc(m31B, dm*dp);
            allocator.dealloc(m31A, dn*dm);
            allocator.dealloc(m30B, dm*dp);
            allocator.dealloc(m30A, dn*dm);
            allocator.dealloc(m29B, dm*dp);
            allocator.dealloc(m29A, dn*dm);
            allocator.dealloc(m28B, dm*dp);
            allocator.dealloc(m28A, dn*dm);
            allocator.dealloc(m27B, dm*dp);
            allocator.dealloc(m27A, dn*dm);
            allocator.dealloc(m26B, dm*dp);
            allocator.dealloc(m26A, dn*dm);
            allocator.dealloc(m25B, dm*dp);
            allocator.dealloc(m25A, dn*dm);
            allocator.dealloc(m24B, dm*dp);
            allocator.dealloc(m24A, dn*dm);
            allocator.dealloc(m23B, dm*dp);
            allocator.dealloc(m23A, dn*dm);
            allocator.dealloc(m22B, dm*dp);
            allocator.dealloc(m22A, dn*dm);
            allocator.dealloc(m21B, dm*dp);
            allocator.dealloc(m21A, dn*dm);
            allocator.dealloc(m20B, dm*dp);
            allocator.dealloc(m20A, dn*dm);
            allocator.dealloc(m19B, dm*dp);
            allocator.dealloc(m19A, dn*dm);
            allocator.dealloc(m18B, dm*dp);
            allocator.dealloc(m18A, dn*dm);
            allocator.dealloc(m17B, dm*dp);
            allocator.dealloc(m17A, dn*dm);
            allocator.dealloc(m16B, dm*dp);
            allocator.dealloc(m16A, dn*dm);
            allocator.dealloc(m15B, dm*dp);
            allocator.dealloc(m15A, dn*dm);
            allocator.dealloc(m14B, dm*dp);
            allocator.dealloc(m14A, dn*dm);
            allocator.dealloc(m13B, dm*dp);
            allocator.dealloc(m13A, dn*dm);
            allocator.dealloc(m12B, dm*dp);
            allocator.dealloc(m12A, dn*dm);
            allocator.dealloc(m11B, dm*dp);
            allocator.dealloc(m11A, dn*dm);
            allocator.dealloc(m10B, dm*dp);
            allocator.dealloc(m10A, dn*dm);
            allocator.dealloc(m9B, dm*dp);
            allocator.dealloc(m9A, dn*dm);
            allocator.dealloc(m8B, dm*dp);
            allocator.dealloc(m8A, dn*dm);
            allocator.dealloc(m7B, dm*dp);
            allocator.dealloc(m7A, dn*dm);
            allocator.dealloc(m6B, dm*dp);
            allocator.dealloc(m6A, dn*dm);
            allocator.dealloc(m5B, dm*dp);
            allocator.dealloc(m5A, dn*dm);
            allocator.dealloc(m4B, dm*dp);
            allocator.dealloc(m4A, dn*dm);
            allocator.dealloc(m3B, dm*dp);
            allocator.dealloc(m3A, dn*dm);
            allocator.dealloc(m2B, dm*dp);
            allocator.dealloc(m2A, dn*dm);
            allocator.dealloc(m1B, dm*dp);
            allocator.dealloc(m1A, dn*dm);
            allocator.dealloc(m40, dn*dp);
            allocator.dealloc(m39, dn*dp);
            allocator.dealloc(m38, dn*dp);
            allocator.dealloc(m37, dn*dp);
            allocator.dealloc(m36, dn*dp);
            allocator.dealloc(m35, dn*dp);
            allocator.dealloc(m34, dn*dp);
            allocator.dealloc(m33, dn*dp);
            allocator.dealloc(m32, dn*dp);
            allocator.dealloc(m31, dn*dp);
            allocator.dealloc(m30, dn*dp);
            allocator.dealloc(m29, dn*dp);
            allocator.dealloc(m28, dn*dp);
            allocator.dealloc(m27, dn*dp);
            allocator.dealloc(m26, dn*dp);
            allocator.dealloc(m25, dn*dp);
            allocator.dealloc(m24, dn*dp);
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
    auto gen6x3x3MinSpace(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 6, 3, 3, 40, detail::Gen6x3x3RecursiveMinSpace>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto gen6x3x3LowLevel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevel>, 6, 3, 3, 40, detail::Gen6x3x3RecursiveLowLevel>(a, b, steps, 0, 0);
    }
    template<int Method = 0, typename M1, typename M2>
    auto gen6x3x3LowLevelParallel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevelParallel>, 6, 3, 3, 40, detail::Gen6x3x3RecursiveLowLevelParallel>(a, b, steps, 40, 40);
    }
}

#endif

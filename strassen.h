#ifndef STRASSEN_H
#define STRASSEN_H
#include "fmmUtility.h"

// Min-Space
namespace fmm::detail {
    struct StrassenMinSpaceRecursive {
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

            operationEff<Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[1][1]);
            operationEff<Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenMinSpaceRecursive>(dC[0][0], tempA, tempB, dn, dm, dp, effC, dm, dp, steps - 1, allocator);
            operationEff<Assign>(dn, dp, effC, effC, dC[1][1], dC[0][0]);

            operationEff<Add>(dn, dm, dm, effA, tempA, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenMinSpaceRecursive>(dC[1][0], tempA, dB[0][0], dn, dm, dp, effC, dm, effB, steps - 1, allocator);
            operationEff<SubAssign>(dn, dp, effC, effC, dC[1][1], dC[1][0]);

            operationEff<Sub>(dm, dp, dp, effB, tempB, dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenMinSpaceRecursive>(dC[0][1], dA[0][0], tempB, dn, dm, dp, effC, effA, dp, steps - 1, allocator);
            operationEff<AddAssign>(dn, dp, effC, effC, dC[1][1], dC[0][1]);

            operationEff<Sub>(dm, dp, dp, effB, tempB, dB[1][0], dB[0][0]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenMinSpaceRecursive>(tempC, dA[1][1], tempB, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            operationsOnFirstArg<AddAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[1][0]);

            operationEff<Add>(dn, dm, dm, effA, tempA, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenMinSpaceRecursive>(tempC, tempA, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            operationsOnFirstArg<SubAssign, AddAssign>(dn, dp, dp, effC, tempC, dC[0][0], dC[0][1]);

            operationEff<Sub>(dn, dm, dm, effA, tempA, dA[1][0], dA[0][0]);
            operationEff<Add>(dm, dp, dp, effB, tempB, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenMinSpaceRecursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationEff<AddAssign>(dn, dp, effC, dp, dC[1][1], tempC);

            operationEff<Sub>(dn, dm, dm, effA, tempA, dA[0][1], dA[1][1]);
            operationEff<Add>(dm, dp, dp, effB, tempB, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenMinSpaceRecursive>(tempC, tempA, tempB, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            operationEff<AddAssign>(dn, dp, effC, dp, dC[0][0], tempC);

            allocator.dealloc(tempC, dn*dp);
            allocator.dealloc(tempB, dm*dp);
            allocator.dealloc(tempA, dn*dm);
        }
    };
}

// Low-Level
namespace fmm::detail {
    struct StrassenLowLevelRecursive {
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
            auto m1_a = operationEffWithAlloc<Add>(allocator, dn, dm, effA, dA[0][0], dA[1][1]);
            auto m1_b = operationEffWithAlloc<Add>(allocator, dm, dp, effB, dB[0][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenLowLevelRecursive>(m1, m1_a, m1_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m1_b, dm*dp);
            allocator.dealloc(m1_a, dn*dm);

            auto m2 = allocator.alloc(dn * dp);
            auto m2_a = operationEffWithAlloc<Add>(allocator, dn, dm, effA, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenLowLevelRecursive>(m2, m2_a, dB[0][0], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m2_a, dn*dm);

            auto m3 = allocator.alloc(dn * dp);
            auto m3_b = operationEffWithAlloc<Sub>(allocator, dm, dp, effB, dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenLowLevelRecursive>(m3, dA[0][0], m3_b, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m3_b, dm*dp);

            auto m4 = allocator.alloc(dn * dp);
            auto m4_b = operationEffWithAlloc<Sub>(allocator, dm, dp, effB, dB[1][0], dB[0][0]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenLowLevelRecursive>(m4, dA[1][1], m4_b, dn, dm, dp, dp, effA, dp, steps - 1, allocator);
            allocator.dealloc(m4_b, dm*dp);

            auto m5 = allocator.alloc(dn * dp);
            auto m5_a = operationEffWithAlloc<Add>(allocator, dn, dm, effA, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenLowLevelRecursive>(m5, m5_a, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1, allocator);
            allocator.dealloc(m5_a, dn*dm);

            auto m6 = allocator.alloc(dn * dp);
            auto m6_a = operationEffWithAlloc<Sub>(allocator, dn, dm, effA, dA[1][0], dA[0][0]);
            auto m6_b = operationEffWithAlloc<Add>(allocator, dm, dp, effB, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenLowLevelRecursive>(m6, m6_a, m6_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
            allocator.dealloc(m6_b, dm*dp);
            allocator.dealloc(m6_a, dn*dm);

            auto m7 = allocator.alloc(dn * dp);
            auto m7_a = operationEffWithAlloc<Sub>(allocator, dn, dm, effA, dA[0][1], dA[1][1]);
            auto m7_b = operationEffWithAlloc<Add>(allocator, dm, dp, effB, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, StrassenLowLevelRecursive>(m7, m7_a, m7_b, dn, dm, dp, dp, dm, dp, steps - 1, allocator);
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
        }
    };
}


// Parallel Low-Level
namespace fmm::detail {
    struct StrassenLowLevelParallelRecursive {
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

            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Add>(dn, dm, dm, effA, m1_a, dA[0][0], dA[1][1]);
                operationEff<Add>(dm, dp, dp, effB, m1_b, dB[0][0], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, StrassenMinSpaceRecursive>(m1, m1_a, m1_b, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Add>(dn, dm, dm, effA, m2_a, dA[1][0], dA[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, StrassenMinSpaceRecursive>(m2, m2_a, dB[0][0], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Sub>(dm, dp, dp, effB, m3_b, dB[0][1], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, StrassenMinSpaceRecursive>(m3, dA[0][0], m3_b, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Sub>(dm, dp, dp, effB, m4_b, dB[1][0], dB[0][0]);
                minSpaceRun<Method, 2, 2, 2, 7, StrassenMinSpaceRecursive>(m4, dA[1][1], m4_b, dn, dm, dp, dp, effA, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Add>(dn, dm, dm, effA, m5_a, dA[0][0], dA[0][1]);
                minSpaceRun<Method, 2, 2, 2, 7, StrassenMinSpaceRecursive>(m5, m5_a, dB[1][1], dn, dm, dp, dp, dm, effB, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Sub>(dn, dm, dm, effA, m6_a, dA[1][0], dA[0][0]);
                operationEff<Add>(dm, dp, dp, effB, m6_b, dB[0][0], dB[0][1]);
                minSpaceRun<Method, 2, 2, 2, 7, StrassenMinSpaceRecursive>(m6, m6_a, m6_b, dn, dm, dp, dp, dm, dp, steps - 1);
            });
            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {
                operationEff<Sub>(dn, dm, dm, effA, m7_a, dA[0][1], dA[1][1]);
                operationEff<Add>(dm, dp, dp, effB, m7_b, dB[1][0], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, StrassenMinSpaceRecursive>(m7, m7_a, m7_b, dn, dm, dp, dp, dm, dp, steps - 1);
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
    };
}



// Static Min-Space
namespace fmm::detail {
     struct StrassenMinSpaceRecursiveStatic {
         template<int Method, int n, int m, int p, int effC, int effA, int effB, int steps, typename T>
         static void Run(T* c, T* a, T* b, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 2;
            constexpr int BaseP = 2;

            constexpr int dn = n / BaseN;
            constexpr int dm = m / BaseM;
            constexpr int dp = p / BaseP;

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto tempA = allocator.alloc(dn * dm);
            auto tempB = allocator.alloc(dm * dp);
            auto tempC = allocator.alloc(dn * dp);

            operationEff<dn, dm, dm, effA, Add>(tempA, dA[0][0], dA[1][1]);
            operationEff<dm, dp, dp, effB, Add>(tempB, dB[0][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, effC, dm, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(dC[0][0], tempA, tempB, allocator);
            operationEff<dn, dp, effC, effC, Assign>(dC[1][1], dC[0][0]);

            operationEff<dn, dm, dm, effA, Add>(tempA, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, effC, dm, effB, steps - 1, StrassenMinSpaceRecursiveStatic>(dC[1][0], tempA, dB[0][0], allocator);
            operationEff<dn, dp, effC, effC, SubAssign>(dC[1][1], dC[1][0]);

            operationEff<dm, dp, dp, effB, Sub>(tempB, dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, effC, effA, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(dC[0][1], dA[0][0], tempB, allocator);
            operationEff<dn, dp, effC, effC, AddAssign>(dC[1][1], dC[0][1]);

            operationEff<dm, dp, dp, effB, Sub>(tempB, dB[1][0], dB[0][0]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, effA, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(tempC, dA[1][1], tempB, allocator);
            operationEff<dn, dp, effC, dp, AddAssign>(dC[0][0], tempC);
            operationEff<dn, dp, effC, dp, AddAssign>(dC[1][0], tempC);

            operationEff<dn, dm, dm, effA, Add>(tempA, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, effB, steps - 1, StrassenMinSpaceRecursiveStatic>(tempC, tempA, dB[1][1], allocator);
            operationEff<dn, dp, effC, dp, SubAssign>(dC[0][0], tempC);
            operationEff<dn, dp, effC, dp, AddAssign>(dC[0][1], tempC);

            operationEff<dn, dm, dm, effA, Sub>(tempA, dA[1][0], dA[0][0]);
            operationEff<dm, dp, dp, effB, Add>(tempB, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(tempC, tempA, tempB, allocator);
            operationEff<dn, dp, effC, dp, AddAssign>(dC[1][1], tempC);

            operationEff<dn, dm, dm, effA, Sub>(tempA, dA[0][1], dA[1][1]);
            operationEff<dm, dp, dp, effB, Add>(tempB, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(tempC, tempA, tempB, allocator);
            operationEff<dn, dp, effC, dp, AddAssign>(dC[0][0], tempC);

            allocator.dealloc(tempC, dn*dp);
            allocator.dealloc(tempB, dm*dp);
            allocator.dealloc(tempA, dn*dm);
        }
    };
}


// Static Low-Level
namespace fmm::detail {
     struct StrassenLowLevelRecursiveStatic {
         template<int Method, int n, int m, int p, int effC, int effA, int effB, int steps, typename T>
         static void Run(T* c, T* a, T* b, StackAllocator<T>& allocator) {
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 2;
            constexpr int BaseP = 2;

            constexpr int dn = n / BaseN;
            constexpr int dm = m / BaseM;
            constexpr int dp = p / BaseP;

            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);
            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);
            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);

            auto m1 = allocator.alloc(dn * dp);
            auto m1_a = operationEffWithAlloc<dn, dm, effA, Add>(allocator, dA[0][0], dA[1][1]);
            auto m1_b = operationEffWithAlloc<dm, dp, effB, Add>(allocator, dB[0][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, dp, steps - 1, StrassenLowLevelRecursiveStatic>(m1, m1_a, m1_b, allocator);
            allocator.dealloc(m1_b, dm*dp);
            allocator.dealloc(m1_a, dn*dm);

            auto m2 = allocator.alloc(dn * dp);
            auto m2_a = operationEffWithAlloc<dn, dm, effA, Add>(allocator, dA[1][0], dA[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, effB, steps - 1, StrassenLowLevelRecursiveStatic>(m2, m2_a, dB[0][0], allocator);
            allocator.dealloc(m2_a, dn*dm);

            auto m3 = allocator.alloc(dn * dp);
            auto m3_b = operationEffWithAlloc<dm, dp, effB, Sub>(allocator, dB[0][1], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, effA, dp, steps - 1, StrassenLowLevelRecursiveStatic>(m3, dA[0][0], m3_b, allocator);
            allocator.dealloc(m3_b, dm*dp);

            auto m4 = allocator.alloc(dn * dp);
            auto m4_b = operationEffWithAlloc<dm, dp, effB, Sub>(allocator, dB[1][0], dB[0][0]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, effA, dp, steps - 1, StrassenLowLevelRecursiveStatic>(m4, dA[1][1], m4_b, allocator);
            allocator.dealloc(m4_b, dm*dp);

            auto m5 = allocator.alloc(dn * dp);
            auto m5_a = operationEffWithAlloc<dn, dm, effA, Add>(allocator, dA[0][0], dA[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, effB, steps - 1, StrassenLowLevelRecursiveStatic>(m5, m5_a, dB[1][1], allocator);
            allocator.dealloc(m5_a, dn*dm);

            auto m6 = allocator.alloc(dn * dp);
            auto m6_a = operationEffWithAlloc<dn, dm, effA, Sub>(allocator, dA[1][0], dA[0][0]);
            auto m6_b = operationEffWithAlloc<dm, dp, effB, Add>(allocator, dB[0][0], dB[0][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, dp, steps - 1, StrassenLowLevelRecursiveStatic>(m6, m6_a, m6_b, allocator);
            allocator.dealloc(m6_b, dm*dp);
            allocator.dealloc(m6_a, dn*dm);

            auto m7 = allocator.alloc(dn * dp);
            auto m7_a = operationEffWithAlloc<dn, dm, effA, Sub>(allocator, dA[0][1], dA[1][1]);
            auto m7_b = operationEffWithAlloc<dm, dp, effB, Add>(allocator, dB[1][0], dB[1][1]);
            nextStep<Method, BaseN, BaseM, BaseP, dn, dm, dp, dp, dm, dp, steps - 1, StrassenLowLevelRecursiveStatic>(m7, m7_a, m7_b, allocator);
            allocator.dealloc(m7_b, dm*dp);
            allocator.dealloc(m7_a, dn*dm);

            operationEff<dn, dp, effC, dp, Assign, Add, Sub, Add>(dC[0][0], m1, m4, m5, m7);
            operationEff<dn, dp, effC, dp, Assign, Add>(dC[0][1], m3, m5);
            operationEff<dn, dp, effC, dp, Assign, Add>(dC[1][0], m2, m4);
            operationEff<dn, dp, effC, dp, Assign, Add, Sub, Add>(dC[1][1], m1, m3, m2, m6);

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


// Static Parallel Low-Level
namespace fmm::detail {
    struct StrassenLowLevelParallelRecursiveStatic {
        template<int Method, int n, int m, int p, int efffC, int efffA, int efffB, int steps, typename T>
        static void Run(T* c, T* a, T* b, StackAllocator<T>& allocator) {
            static constexpr int effC = efffC;
            static constexpr int effA = efffA;
            static constexpr int effB = efffB;
            
            using namespace ArithmeticOperation;

            constexpr int BaseN = 2;
            constexpr int BaseM = 2;
            constexpr int BaseP = 2;

            static constexpr int dn = n / BaseN;
            static constexpr int dm = m / BaseM;
            static constexpr int dp = p / BaseP;

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
                operationEff<dn, dm, dm, effA, Add>(m1_a, dA[0][0], dA[1][1]);
                operationEff<dm, dp, dp, effB, Add>(m1_b, dB[0][0], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, dn, dm, dp, dp, dm, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(m1, m1_a, m1_b);
            });
            pool.addTask([=]() {
                operationEff<dn, dm, dm, effA, Add>(m2_a, dA[1][0], dA[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, dn, dm, dp, dp, dm, effB, steps - 1, StrassenMinSpaceRecursiveStatic>(m2, m2_a, dB[0][0]);
            });
            pool.addTask([=]() {
                operationEff<dm, dp, dp, effB, Sub>(m3_b, dB[0][1], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, dn, dm, dp, dp, effA, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(m3, dA[0][0], m3_b);
            });
            pool.addTask([=]() {
                operationEff<dm, dp, dp, effB, Sub>(m4_b, dB[1][0], dB[0][0]);
                minSpaceRun<Method, 2, 2, 2, 7, dn, dm, dp, dp, effA, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(m4, dA[1][1], m4_b);
            });
            pool.addTask([=]() {
                operationEff<dn, dm, dm, effA, Add>(m5_a, dA[0][0], dA[0][1]);
                minSpaceRun<Method, 2, 2, 2, 7, dn, dm, dp, dp, dm, effB, steps - 1, StrassenMinSpaceRecursiveStatic>(m5, m5_a, dB[1][1]);
            });
            pool.addTask([=]() {
                operationEff<dn, dm, dm, effA, Sub>( m6_a, dA[1][0], dA[0][0]);
                operationEff<dm, dp, dp, effB, Add>(m6_b, dB[0][0], dB[0][1]);
                minSpaceRun<Method, 2, 2, 2, 7, dn, dm, dp, dp, dm, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(m6, m6_a, m6_b);
            });
            pool.addTask([=]() {
                operationEff<dn, dm, dm, effA, Sub>(m7_a, dA[0][1], dA[1][1]);
                operationEff<dm, dp, dp, effB, Add>(m7_b, dB[1][0], dB[1][1]);
                minSpaceRun<Method, 2, 2, 2, 7, dn, dm, dp, dp, dm, dp, steps - 1, StrassenMinSpaceRecursiveStatic>(m7, m7_a, m7_b);
            });

            pool.completeTasksAndStop();

            operationEff<dn, dp, effC, dp, Assign, Add, Sub, Add>(dC[0][0], m1, m4, m5, m7);
            operationEff<dn, dp, effC, dp, Assign, Add>(dC[0][1], m3, m5);
            operationEff<dn, dp, effC, dp, Assign, Add>(dC[1][0], m2, m4);
            operationEff<dn, dp, effC, dp, Assign, Add, Sub, Add>(dC[1][1], m1, m3, m2, m6);

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
    };
}



namespace fmm {
    template<int Method=0, typename M1, typename M2>
    auto strassenMinSpace(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 2, 2, 2, 7, detail::StrassenMinSpaceRecursive>(a, b, steps);
    }

    template<int Method=0, typename M1, typename M2>
    auto strassenLowLevel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevel>, 2, 2, 2, 7, detail::StrassenLowLevelRecursive>(a, b, steps);
    }

    template<int Method=0, typename M1, typename M2>
    auto strassenParallelLowLevel(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {
        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::LowLevelParallel>, 2, 2, 2, 7, detail::StrassenLowLevelParallelRecursive>(a, b, steps);
    }

    template<int steps, int Method=0, typename M1, typename M2>
    auto strassenMinSpaceStatic(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        return detail::runAlgorithmStatic<getNewWithAlgorithm<Method, Algorithm::MinSpace>, 2, 2, 2, 7, steps, detail::StrassenMinSpaceRecursiveStatic>(a, b);
    }

    template<int steps, int Method=0, typename M1, typename M2>
    auto strassenLowLevelStatic(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        return detail::runAlgorithmStatic<getNewWithAlgorithm<Method, Algorithm::LowLevel>, 2, 2, 2, 7, steps, detail::StrassenLowLevelRecursiveStatic>(a, b);
    }

    template<int steps, int Method=0, typename M1, typename M2>
    auto strassenParallelLowLevelStatic(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        return detail::runAlgorithmStatic<getNewWithAlgorithm<Method, Algorithm::LowLevelParallel>, 2, 2, 2, 7, steps, detail::StrassenLowLevelParallelRecursiveStatic>(a, b);
    }
}



#endif

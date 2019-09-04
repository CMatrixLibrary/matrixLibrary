#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#if __has_include("mkl.h")
    #include "mkl.h"
    #define USE_BLAS
#endif
#include "blasMul.h"
#include "Matrix.h"
#include "MatrixView.h"
#include "naiveBasicOperations.h"
#include "matrixOperators.h"
#include "Range.h"
#include "RangeZip.h"
#include "MatrixExtendedFunctions.h"
#include "benchmark.h"
#include "strassen.h"
#include "blockMul.h"
#include "parallelMul.h"
#include "parallelBlockMul.h"
#include "fmmUtility.h"

template<typename T, int Reps, int N, int NEnd, int Inc, bool UseBlas = true, bool UseAvx = true> void testBaseMul(std::ostream& out);

template<typename T, int Reps, int N, int NEnd, int Inc, bool UseBlas, bool UseAvx>
struct TestBaseMulImpl {
    static void TestBaseMul(std::ostream& out) {
        std::cout << "\rN = " << N << "/" << NEnd-1;
        HeapMatrix<T> dynamicA(N, N);
        HeapMatrix<T> dynamicB(N, N);
        StaticHeapMatrix<T, N, N> staticA;
        StaticHeapMatrix<T, N, N> staticB;
        out << N;
        out << ',' << benchmark(Reps, [&]() { naiveMul(dynamicA, dynamicB); });
        out << ',' << benchmark(Reps, [&]() { blockMul(dynamicA, dynamicB); });
        out << ',' << benchmark(Reps, [&]() { parallelMul(dynamicA, dynamicB); });
        out << ',' << benchmark(Reps, [&]() { parallelBlockMul(dynamicA, dynamicB); });
        if constexpr (UseAvx) out << ',' << benchmark(Reps, [&]() { avx::mul(dynamicA, dynamicB); });
        if constexpr (UseAvx) out << ',' << benchmark(Reps, [&]() { avx::parallelMul(dynamicA, dynamicB); });
        if constexpr (UseBlas) out << ',' << benchmark(Reps, [&]() { blas::mul(dynamicA, dynamicB); });
        out << ',' << benchmark(Reps, [&]() { naiveMul(staticA, staticB); });
        out << ',' << benchmark(Reps, [&]() { blockMul(staticA, staticB); });
        out << ',' << benchmark(Reps, [&]() { parallelMul(staticA, staticB); });
        out << ',' << benchmark(Reps, [&]() { parallelBlockMul(staticA, staticB); });
        if constexpr (UseAvx) out << ',' << benchmark(Reps, [&]() { avx::mul(staticA, staticB); });
        if constexpr (UseAvx) out << ',' << benchmark(Reps, [&]() { avx::parallelMul(staticA, staticB); });
        if constexpr (UseBlas) out << ',' << benchmark(Reps, [&]() { blas::mul(staticA, staticB); });
        out << std::endl;
        ::testBaseMul<T, Reps, N+Inc, NEnd, Inc, UseBlas, UseAvx>(out);
    }
};

template<typename T, int Reps, int N, int Inc, bool UseBlas, bool UseAvx>
struct TestBaseMulImpl<T, Reps, N, N, Inc, UseBlas, UseAvx> {
    static void TestBaseMul(std::ostream& out) {}
};

template<typename T, int Reps, int N, int NEnd, int Inc, bool UseBlas, bool UseAvx>
void testBaseMul(std::ostream& out) {
    TestBaseMulImpl<T, Reps, N, NEnd, Inc, UseBlas, UseAvx>::TestBaseMul(out);
}
template<typename T, int Reps, int N, int NEnd, int Inc, bool UseBlas=true, bool UseAvx=true>
void testBaseMul(bool append = false, const std::string& id="") {
    std::ofstream outFile("BaseMulTest" + id + ".csv", append ? std::ios_base::app : std::ios_base::out);
    if (!append) {
        outFile << "N";
        outFile << ",D_Naive,D_Block,D_Parallel,D_ParallelBlock";
        if constexpr (UseAvx) outFile << ",D_Avx,D_ParallelAvx";
        if constexpr (UseBlas) outFile << ",D_Blas";
        outFile << ",S_Naive,S_Block,S_Parallel,S_ParallelBlock";
        if constexpr (UseAvx) outFile << ",S_Avx,S_ParallelAvx";
        if constexpr (UseBlas) outFile << ",S_Blas";
        outFile << std::endl;
    }
    std::cout << '\n';
    testBaseMul<T, Reps, N, NEnd, Inc, UseBlas, UseAvx>(outFile);
    std::cout << "\nDone\n";
}



void strassenBenchmarkPrintProgress(int n, int steps, int current, int max) {
    std::cout << "\rProcessing [" << n << ' ' << steps << "] " << current << "/" << max;
}

template<int N, int Steps, int Method, typename T> void strassenBenchmarkMethod(int id, int reps, std::ostream& out, HeapMatrix<T>& a, HeapMatrix<T>& b) {
    strassenBenchmarkPrintProgress(N, Steps, id, 42);
    out << ',';
    if constexpr (fmm::contains<Method, fmm::LowLevel>) 
        out << benchmark(reps, [&]() {fmm::strassenLowLevel<Method>(a, b, Steps); });
    if constexpr (fmm::contains<Method, fmm::MinSpace>)
        out << benchmark(reps, [&]() {fmm::strassenMinSpace<Method>(a, b, Steps); });
    if constexpr (fmm::contains<Method, fmm::LowLevelParallel>)
        out << benchmark(reps, [&]() {fmm::strassenParallelLowLevel<Method>(a, b, Steps); });
}

int strassenBenchmarkId = 0;
template<typename T, int N, int Steps, bool UseNaive=true> void strassenBenchmark(std::ostream& out, int repetitions) {
    HeapMatrix<T> a(N, N);
    HeapMatrix<T> b(N, N);
    for (int i = 0; i < a.size(); ++i) {
        a.data()[i] = i;
        b.data()[i] = i * 2;
    }
    std::cout << '\n';

    constexpr int BaseMethod = fmm::DynamicPeeling;

    constexpr int E_Method = BaseMethod | fmm::Effective;
    constexpr int N_Method = BaseMethod | fmm::Normal;

    constexpr int L_E_Method = E_Method | fmm::LowLevel;
    constexpr int M_E_Method = E_Method | fmm::MinSpace;
    constexpr int P_E_Method = E_Method | fmm::LowLevelParallel;
    constexpr int L_N_Method = N_Method | fmm::LowLevel;
    constexpr int M_N_Method = N_Method | fmm::MinSpace;
    constexpr int P_N_Method = N_Method | fmm::LowLevelParallel;

    int id = 0;
    strassenBenchmarkPrintProgress(N, Steps, id, 42);
    
    out << strassenBenchmarkId++ << ',' << N << ',' << Steps;

    if constexpr (Steps == 0) {
        out << ',' << benchmark(repetitions, [&]() {avx::mul(a, b); });
        out << ',' << benchmark(repetitions, [&]() {avx::parallelMul(a, b); });
        //out << ',' << benchmark(repetitions, [&]() {blas::mul(a, b); });
        out << ',';
        return;
    }

    if constexpr (UseNaive) {
        strassenBenchmarkMethod<N, Steps, fmm::Naive | L_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Naive | L_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Naive | M_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Naive | M_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Naive | P_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Naive | P_N_Method>(++id, repetitions, out, a, b);

        strassenBenchmarkMethod<N, Steps, fmm::Block | L_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Block | L_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Block | M_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Block | M_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Block | P_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Block | P_N_Method>(++id, repetitions, out, a, b);

        strassenBenchmarkMethod<N, Steps, fmm::Parallel | L_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Parallel | L_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Parallel | M_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Parallel | M_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Parallel | P_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::Parallel | P_N_Method>(++id, repetitions, out, a, b);

        strassenBenchmarkMethod<N, Steps, fmm::ParallelBlock | L_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::ParallelBlock | L_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::ParallelBlock | M_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::ParallelBlock | M_N_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::ParallelBlock | P_E_Method>(++id, repetitions, out, a, b);
        strassenBenchmarkMethod<N, Steps, fmm::ParallelBlock | P_N_Method>(++id, repetitions, out, a, b);
    } else {
        out << ",-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-";
    }
    strassenBenchmarkMethod<N, Steps, fmm::Avx | L_E_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::Avx | L_N_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::Avx | M_E_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::Avx | M_N_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::Avx | P_E_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::Avx | P_N_Method>(++id, repetitions, out, a, b);

    strassenBenchmarkMethod<N, Steps, fmm::ParallelAvx | L_E_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::ParallelAvx | L_N_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::ParallelAvx | M_E_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::ParallelAvx | M_N_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::ParallelAvx | P_E_Method>(++id, repetitions, out, a, b);
    strassenBenchmarkMethod<N, Steps, fmm::ParallelAvx | P_N_Method>(++id, repetitions, out, a, b);

    //strassenBenchmarkMethod<N, Steps, fmm::Blas | L_E_Method>(++id, repetitions, out, a, b);
    //strassenBenchmarkMethod<N, Steps, fmm::Blas | L_N_Method>(++id, repetitions, out, a, b);
    //strassenBenchmarkMethod<N, Steps, fmm::Blas | M_E_Method>(++id, repetitions, out, a, b);
    //strassenBenchmarkMethod<N, Steps, fmm::Blas | M_N_Method>(++id, repetitions, out, a, b);
    //strassenBenchmarkMethod<N, Steps, fmm::Blas | P_E_Method>(++id, repetitions, out, a, b);
    //strassenBenchmarkMethod<N, Steps, fmm::Blas | P_N_Method>(++id, repetitions, out, a, b);
    
    //out << std::endl;
    out << ',';
}

template<typename T> void strassenBenchmark() {
    std::ofstream outFile("strassenBenchmark.csv");
    outFile << "ID,N,Steps"; 
    for (auto baseMul : std::array<std::string, 7>{"Naive", "Block", "Parallel", "ParallelBlock", "Avx", "ParallelAvx", "Blas"}) {
        for (auto algorithm : std::array<std::string, 3>{"Low-Level", "Min-Space", "Parallel-Low-Level"}) {
            for (auto baseSize : std::array<std::string, 2>{"Effective", "Normals"}) {
                outFile << "," << algorithm << "_" << baseMul << "_" << baseSize;
            }
        }
    }
    outFile << '\n';
    //strassenBenchmark<T, 2, 0>(outFile, 100);
    //strassenBenchmark<T, 2, 1>(outFile, 100);
    //strassenBenchmark<T, 4, 0>(outFile, 100);
    //strassenBenchmark<T, 4, 1>(outFile, 100);
    //strassenBenchmark<T, 4, 2>(outFile, 100);
    //strassenBenchmark<T, 8, 0>(outFile, 10);
    //strassenBenchmark<T, 8, 1>(outFile, 10);
    //strassenBenchmark<T, 8, 2>(outFile, 10);
    //strassenBenchmark<T, 16, 0>(outFile, 10);
    //strassenBenchmark<T, 16, 1>(outFile, 10);
    //strassenBenchmark<T, 16, 2>(outFile, 10);
    //strassenBenchmark<T, 32, 0>(outFile, 10);
    //strassenBenchmark<T, 32, 1>(outFile, 10);
    //strassenBenchmark<T, 32, 2>(outFile, 10);
    //strassenBenchmark<T, 64, 0>(outFile, 10);
    //strassenBenchmark<T, 64, 1>(outFile, 10);
    //strassenBenchmark<T, 64, 2>(outFile, 10);
    //strassenBenchmark<T, 128, 0>(outFile, 10);
    //strassenBenchmark<T, 128, 1>(outFile, 10);
    //strassenBenchmark<T, 128, 2>(outFile, 10);
    //strassenBenchmark<T, 128, 3>(outFile, 10);
    //strassenBenchmark<T, 256, 0>(outFile, 10);
    //strassenBenchmark<T, 256, 1>(outFile, 10);
    //strassenBenchmark<T, 256, 2>(outFile, 10);
    //strassenBenchmark<T, 256, 3>(outFile, 10);
    //strassenBenchmark<T, 256, 4>(outFile, 10);
    //strassenBenchmark<T, 512, 0>(outFile, 5);
    //strassenBenchmark<T, 512, 1>(outFile, 5);
    //strassenBenchmark<T, 512, 2>(outFile, 5);
    //strassenBenchmark<T, 512, 3>(outFile, 5);
    //strassenBenchmark<T, 512, 4>(outFile, 5);
    //strassenBenchmark<T, 512, 5>(outFile, 5);
    //strassenBenchmark<T, 1024, 0>(outFile, 1);
    //strassenBenchmark<T, 1024, 1>(outFile, 2);
    //strassenBenchmark<T, 1024, 2>(outFile, 2);
    //strassenBenchmark<T, 1024, 3>(outFile, 2);
    //strassenBenchmark<T, 1024, 4>(outFile, 2);
    //strassenBenchmark<T, 1024, 5>(outFile, 2);
    //strassenBenchmark<T, 1024, 6>(outFile, 2);
    strassenBenchmark<T, 4096, 0, false>(outFile, 2);
    strassenBenchmark<T, 4096, 1, false>(outFile, 2);
    strassenBenchmark<T, 4096, 2, false>(outFile, 2);
    strassenBenchmark<T, 4096, 3, false>(outFile, 2);
    strassenBenchmark<T, 4096, 3, false>(outFile, 2);
    strassenBenchmark<T, 4096, 4, false>(outFile, 2);
    strassenBenchmark<T, 4096, 5, false>(outFile, 2);

    std::cout << "\nDONE\n";
}

int main() {
    using T = float;
    constexpr int N = 128;
    constexpr int Steps = 2;

    StaticHeapMatrix<T, N, N> a;
    StaticHeapMatrix<T, N, N> b;
    fmm::strassenLowLevelStatic<Steps>(a, b);
    fmm::strassenMinSpaceStatic<Steps>(a, b);
    fmm::strassenParallelLowLevelStatic<Steps>(a, b);

    strassenBenchmark<float>();
    std::cin.get();
    return 0;
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
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

template<typename T, typename Function> void strassenDynamicBenchmark(std::string label, int n, int steps, Function fun) {
    HeapMatrix<T> a(n, n);
    HeapMatrix<T> b(n, n);
    for (int i = 0; i < a.size(); ++i) {
        a.data()[i] = i;
        b.data()[i] = i * 2;
    }
    auto start = MyTime();
    auto c = fun(a, b, steps);
    auto end = MyTime();
    std::cout << "Dynamic " << label << ' ' << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
}
template<typename T, int N, typename Function> void strassenStaticBenchmark(std::string label, Function fun) {
    StaticHeapMatrix<T, N, N> a;
    StaticHeapMatrix<T, N, N> b;
    for (int i = 0; i < a.size(); ++i) {
        a.data()[i] = i;
        b.data()[i] = i * 2;
    }
    auto start = MyTime();
    auto c = fun(a, b);
    auto end = MyTime();
    std::cout << "Static  " << label << ' ' << std::setw(10) << end - start << "; data[14] = " << c.data()[14] << '\n';
}
template<typename T, int N, int Steps, bool OnlyAvx=false, bool OnlyDynamic=false> void strassenTest(int n) {
    using StaticMatrix = StaticHeapMatrix<T, N, N>;
    using DynamicMatrix = HeapMatrix<T>;
    
    if (n != N) {
        std::cout << "incorrect stdin\n";
        return;
    }
    std::cout << std::left;
    std::cout << "=================================================================\n";
    std::cout << "N               : " << N << '\n';
    std::cout << "Recursive Steps : " << Steps << '\n';
    std::cout << "-------------------------------\n";
    if constexpr (!OnlyAvx) {
        strassenDynamicBenchmark<T>("High-level noAVX :", n, Steps, [](const auto& a, const auto& b, int steps) {
            return strassen(a, b, steps);
        });
        if constexpr (!OnlyDynamic) strassenStaticBenchmark<T, N>("High-level noAVX :", [](const auto& a, const auto& b) {
            return strassen<Steps>(a, b);
        });
        strassenDynamicBenchmark<T>("Low-level  noAVX :", n, Steps, [](const auto& a, const auto& b, int steps) {
            return lowLevelStrassen(a, b, steps);
        });
        if constexpr (!OnlyDynamic)  strassenStaticBenchmark<T, N>("Low-level  noAVX :", [](const auto& a, const auto& b) {
            return lowLevelStrassen<Steps>(a, b);
        });
        strassenDynamicBenchmark<T>("Min-Space  noAVX :", n, Steps, [](const auto& a, const auto& b, int steps) {
            return minSpaceStrassen(a, b, steps);
        });
        if constexpr (!OnlyDynamic)  strassenStaticBenchmark<T, N>("Min-Space  noAVX :", [](const auto& a, const auto& b) {
            return minSpaceStrassen<Steps>(a, b);
        });
    }
    strassenDynamicBenchmark<T>("High-level AVX   :", n, Steps, [](const auto& a, const auto& b, int steps) {
        return strassenAvx(a, b, steps);
    });
    if constexpr (!OnlyDynamic) strassenStaticBenchmark<T, N>("High-level AVX   :", [](const auto& a, const auto& b) {
        return strassenAvx<Steps>(a, b);
    });
    strassenDynamicBenchmark<T>("Low-level  AVX   :", n, Steps, [](const auto& a, const auto& b, int steps) {
        return lowLevelAvxStrassen(a, b, steps);
    });
    if constexpr (!OnlyDynamic) strassenStaticBenchmark<T, N>("Low-level  AVX   :", [](const auto& a, const auto& b) {
        return lowLevelAvxStrassen<Steps>(a, b);
    });
    strassenDynamicBenchmark<T>("Min-Space  AVX   :", n, Steps, [](const auto& a, const auto& b, int steps) {
        return minSpaceAvxStrassen(a, b, steps);
    });
    if constexpr (!OnlyDynamic) strassenStaticBenchmark<T, N>("Min-Space  AVX   :", [](const auto& a, const auto& b) {
        return minSpaceAvxStrassen<Steps>(a, b);
    });
    std::cout << "=================================================================\n\n";
}

template<typename T> void strassenTestForType() {
    std::array<int, 6> vec;
    std::cout << "copy line below and click enter: \n";
    std::cout << "480 992 1024 1920 2048 3840\n";
    std::cin >> vec[0] >> vec[1] >> vec[2] >> vec[3] >> vec[4] >> vec[5];
    strassenTest<T, 480, 0>(vec[0]);
    strassenTest<T, 480, 1>(vec[0]);
    strassenTest<T, 480, 2>(vec[0]);
    strassenTest<T, 992, 0>(vec[1]);
    strassenTest<T, 992, 1>(vec[1]);
    strassenTest<T, 992, 2>(vec[1]);
    strassenTest<T, 992, 3>(vec[1]);
    strassenTest<T, 992, 4>(vec[1]);
    strassenTest<T, 992, 5>(vec[1]);
    strassenTest<T, 1024, 0>(vec[2]);
    strassenTest<T, 1024, 1>(vec[2]);
    strassenTest<T, 1024, 2>(vec[2]);
    strassenTest<T, 1920, 0>(vec[3]);
    strassenTest<T, 1920, 1>(vec[3]);
    strassenTest<T, 1920, 2>(vec[3]);
    strassenTest<T, 1920, 3>(vec[3]);
    strassenTest<T, 1920, 4>(vec[3]);
    strassenTest<T, 1920, 5>(vec[3]);
    strassenTest<T, 2048, 2>(vec[4]);
    strassenTest<T, 2048, 3>(vec[4]);
    strassenTest<T, 3840, 0, true>(vec[5]);
    strassenTest<T, 3840, 1, true>(vec[5]);
    strassenTest<T, 3840, 2, true>(vec[5]);
}
template<typename T> void strassenTestWithStaticPadding() {
    std::array<int, 8> vec;
    std::cout << "copy line below and click enter: \n";
    std::cout << "479 480 991 992 1913 1920 3833 3840\n";
    std::cin >> vec[0] >> vec[1] >> vec[2] >> vec[3] >> vec[4] >> vec[5] >> vec[6] >> vec[7];
    strassenTest<T, 479, 1, false, true>(vec[0]);
    strassenTest<T, 480, 1, false, true>(vec[1]);
    strassenTest<T, 479, 2, false, true>(vec[0]);
    strassenTest<T, 480, 2, false, true>(vec[1]);

    strassenTest<T, 991, 1, false, true>(vec[2]);
    strassenTest<T, 992, 1, false, true>(vec[3]);
    strassenTest<T, 991, 2, false, true>(vec[2]);
    strassenTest<T, 992, 2, false, true>(vec[3]);
    strassenTest<T, 991, 3, false, true>(vec[2]);
    strassenTest<T, 992, 3, false, true>(vec[3]);
    strassenTest<T, 991, 4, false, true>(vec[2]);
    strassenTest<T, 992, 4, false, true>(vec[3]);
    strassenTest<T, 991, 5, false, true>(vec[2]);
    strassenTest<T, 992, 5, false, true>(vec[3]);

    strassenTest<T, 1913, 1, false, true>(vec[4]);
    strassenTest<T, 1920, 1, false, true>(vec[5]);
    strassenTest<T, 1913, 2, false, true>(vec[4]);
    strassenTest<T, 1920, 2, false, true>(vec[5]);
    strassenTest<T, 1913, 3, false, true>(vec[4]);
    strassenTest<T, 1920, 3, false, true>(vec[5]);
    strassenTest<T, 1913, 4, false, true>(vec[4]);
    strassenTest<T, 1920, 4, false, true>(vec[5]);
    strassenTest<T, 1913, 5, false, true>(vec[4]);
    strassenTest<T, 1920, 5, false, true>(vec[5]);

    strassenTest<T, 3833, 1, true, true>(vec[6]);
    strassenTest<T, 3840, 1, true, true>(vec[7]);
    strassenTest<T, 3833, 2, true, true>(vec[6]);
    strassenTest<T, 3840, 2, true, true>(vec[7]);
}


void strassenSquereBenchmarkPrintProgress(int n, int steps, int current, int max) {
    std::cout << "\rProcessing [" << n << ' ' << steps << "] " << current << "/" << max;
}
int strassenSquereBenchmarkId = 0;
template<typename T, int N, int Steps, bool UseNaive=true> void strassenSquereBenchmark(std::ostream& out, int repetitions) {
    HeapMatrix<T> dynamicA(N, N);
    HeapMatrix<T> dynamicB(N, N);
    for (int i = 0; i < dynamicA.size(); ++i) {
        dynamicA.data()[i] = i;
        dynamicB.data()[i] = i * 2;
    }
    StaticHeapMatrix<T, N, N> staticA;
    StaticHeapMatrix<T, N, N> staticB;
    for (int i = 0; i < staticA.size(); ++i) {
        staticA.data()[i] = i;
        staticB.data()[i] = i * 2;
    }
    std::cout << '\n';
    out << strassenSquereBenchmarkId++ << ',' << N << ',' << Steps << ',';
    if constexpr (UseNaive) {
        strassenSquereBenchmarkPrintProgress(N, Steps, 0, 18);
        out << benchmark(repetitions, [&]() { strassen(dynamicA, dynamicB, Steps); }) << ',';
        strassenSquereBenchmarkPrintProgress(N, Steps, 1, 18);
        out << benchmark(repetitions, [&]() { lowLevelStrassen(dynamicA, dynamicB, Steps); }) << ',';
        strassenSquereBenchmarkPrintProgress(N, Steps, 2, 18);
        out << benchmark(repetitions, [&]() { minSpaceStrassen(dynamicA, dynamicB, Steps); }) << ',';
        strassenSquereBenchmarkPrintProgress(N, Steps, 3, 18);
        out << benchmark(repetitions, [&]() { strassen<Steps>(staticA, staticB); }) << ',';
        strassenSquereBenchmarkPrintProgress(N, Steps, 4, 18);
        out << benchmark(repetitions, [&]() { lowLevelStrassen<Steps>(staticA, staticB); }) << ',';
        strassenSquereBenchmarkPrintProgress(N, Steps, 5, 18);
        out << benchmark(repetitions, [&]() { minSpaceStrassen<Steps>(staticA, staticB); }) << ',';
    } else {
        out << "-,-,-,-,-,-";
    }
    strassenSquereBenchmarkPrintProgress(N, Steps, 6, 18);
    out << benchmark(repetitions, [&]() { strassenAvx(dynamicA, dynamicB, Steps); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 7, 18);
    out << benchmark(repetitions, [&]() { lowLevelAvxStrassen(dynamicA, dynamicB, Steps); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 8, 18);
    out << benchmark(repetitions, [&]() { minSpaceAvxStrassen(dynamicA, dynamicB, Steps); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 9, 18);
    out << benchmark(repetitions, [&]() { strassenAvx<Steps>(staticA, staticB); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 10, 18);
    out << benchmark(repetitions, [&]() { lowLevelAvxStrassen<Steps>(staticA, staticB); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 11, 18);
    out << benchmark(repetitions, [&]() { minSpaceAvxStrassen<Steps>(staticA, staticB); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 12, 18);
    out << benchmark(repetitions, [&]() { strassenAuto(dynamicA, dynamicB, Steps); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 13, 18);
    out << benchmark(repetitions, [&]() { lowLevelStrassen<BaseMulType::Blas>(dynamicA, dynamicB, Steps); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 14, 18);
    out << benchmark(repetitions, [&]() { minSpaceStrassen<BaseMulType::Blas>(dynamicA, dynamicB, Steps); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 15, 18);
    out << benchmark(repetitions, [&]() { strassenAuto<Steps>(staticA, staticB); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 16, 18);
    out << benchmark(repetitions, [&]() { lowLevelStrassen<BaseMulType::Blas, Steps>(staticA, staticB); }) << ',';
    strassenSquereBenchmarkPrintProgress(N, Steps, 17, 18);
    out << benchmark(repetitions, [&]() { minSpaceStrassen<Steps, BaseMulType::Blas>(staticA, staticB); });
    strassenSquereBenchmarkPrintProgress(N, Steps, 18, 18);
    out << std::endl;
}
template<typename T> void strassenSquereBenchmark() {
    std::ofstream outFile("StrassenSquereBenchmark.csv");
    outFile << "ID,N,Steps,HDN,LDN,MDN,HSN,LSN,MSN,HDA,LDA,MDA,HSA,LSA,MSA,HDB,LDB,MDB,HSB,LSB,MSB\n";
    strassenSquereBenchmark<T, 2, 0>(outFile, 10000);
    strassenSquereBenchmark<T, 2, 1>(outFile, 10000);
    strassenSquereBenchmark<T, 4, 0>(outFile, 10000);
    strassenSquereBenchmark<T, 4, 1>(outFile, 10000);
    strassenSquereBenchmark<T, 4, 2>(outFile, 10000);
    strassenSquereBenchmark<T, 8, 0>(outFile, 1000);
    strassenSquereBenchmark<T, 8, 1>(outFile, 1000);
    strassenSquereBenchmark<T, 8, 2>(outFile, 1000);
    strassenSquereBenchmark<T, 8, 3>(outFile, 1000);
    strassenSquereBenchmark<T, 16, 0>(outFile, 100);
    strassenSquereBenchmark<T, 16, 1>(outFile, 100);
    strassenSquereBenchmark<T, 16, 2>(outFile, 100);
    strassenSquereBenchmark<T, 16, 3>(outFile, 100);
    strassenSquereBenchmark<T, 16, 4>(outFile, 100);
    strassenSquereBenchmark<T, 32, 0>(outFile, 100);
    strassenSquereBenchmark<T, 32, 1>(outFile, 100);
    strassenSquereBenchmark<T, 32, 2>(outFile, 100);
    strassenSquereBenchmark<T, 32, 3>(outFile, 100);
    strassenSquereBenchmark<T, 32, 4>(outFile, 100);
    strassenSquereBenchmark<T, 32, 5>(outFile, 100);
    strassenSquereBenchmark<T, 64, 0>(outFile, 10);
    strassenSquereBenchmark<T, 64, 1>(outFile, 10);
    strassenSquereBenchmark<T, 64, 2>(outFile, 10);
    strassenSquereBenchmark<T, 64, 3>(outFile, 10);
    strassenSquereBenchmark<T, 64, 4>(outFile, 10);
    strassenSquereBenchmark<T, 64, 5>(outFile, 10);
    strassenSquereBenchmark<T, 64, 6>(outFile, 10);
    strassenSquereBenchmark<T, 128, 0>(outFile, 10);
    strassenSquereBenchmark<T, 128, 1>(outFile, 10);
    strassenSquereBenchmark<T, 128, 2>(outFile, 10);
    strassenSquereBenchmark<T, 128, 3>(outFile, 10);
    strassenSquereBenchmark<T, 128, 4>(outFile, 10);
    strassenSquereBenchmark<T, 128, 5>(outFile, 10);
    strassenSquereBenchmark<T, 128, 6>(outFile, 10);
    strassenSquereBenchmark<T, 256, 0>(outFile, 10);
    strassenSquereBenchmark<T, 256, 1>(outFile, 10);
    strassenSquereBenchmark<T, 256, 2>(outFile, 10);
    strassenSquereBenchmark<T, 256, 3>(outFile, 10);
    strassenSquereBenchmark<T, 256, 4>(outFile, 10);
    strassenSquereBenchmark<T, 256, 5>(outFile, 10);
    strassenSquereBenchmark<T, 512, 0>(outFile, 5);
    strassenSquereBenchmark<T, 512, 1>(outFile, 5);
    strassenSquereBenchmark<T, 512, 2>(outFile, 5);
    strassenSquereBenchmark<T, 512, 3>(outFile, 5);
    strassenSquereBenchmark<T, 512, 4>(outFile, 5);
    strassenSquereBenchmark<T, 512, 5>(outFile, 5);
    strassenSquereBenchmark<T, 1024, 0>(outFile, 1);
    strassenSquereBenchmark<T, 1024, 1>(outFile, 2);
    strassenSquereBenchmark<T, 1024, 2>(outFile, 2);
    strassenSquereBenchmark<T, 1024, 3>(outFile, 2);
    strassenSquereBenchmark<T, 1024, 4>(outFile, 2);
    strassenSquereBenchmark<T, 1024, 5>(outFile, 2);
    strassenSquereBenchmark<T, 1024, 6>(outFile, 2);
    strassenSquereBenchmark<T, 2048, 0, false>(outFile, 1);
    strassenSquereBenchmark<T, 2048, 1>(outFile, 1);
    strassenSquereBenchmark<T, 2048, 2>(outFile, 1);
    strassenSquereBenchmark<T, 2048, 3>(outFile, 1);
    strassenSquereBenchmark<T, 2048, 4>(outFile, 1);
    strassenSquereBenchmark<T, 2048, 5>(outFile, 1);
    strassenSquereBenchmark<T, 2048, 6>(outFile, 1);
    strassenSquereBenchmark<T, 2048, 7>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 0, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 1, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 2, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 3, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 4, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 5, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 6, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 7, false>(outFile, 1);
    strassenSquereBenchmark<T, 4096, 8, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 0, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 1, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 2, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 3, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 4, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 5, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 6, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 7, false>(outFile, 1);
    strassenSquereBenchmark<T, 8192, 8, false>(outFile, 1);
    std::cout << "\nDONE\n";
}

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

int main() {
    strassenSquereBenchmark<float>();
    //strassenTestWithStaticPadding<float>();
    std::cin.get();
    return 0;
}

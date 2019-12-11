#ifndef MATRIX_CSV_BENCHMARK_H
#define MATRIX_CSV_BENCHMARK_H

#include "Matrix.h"
#include "benchmark.h"
#include <iostream>
#include <fstream>
#include <iomanip>

enum class IncrementType {
    Add,
    Mul
};

namespace detail {
    template<typename T, int Reps, int N, int NInc, int M, int MInc, int P, int PInc, int step, int endStep, IncrementType incType, int InputMatrixCount, bool onlyDynamic, typename... BenchmarkFunctions>
    void matrixBenchmark(std::ostream& out);

    template<typename T, int Reps, int N, int NInc, int M, int MInc, int P, int PInc, int step, int endStep, IncrementType incType, int InputMatrixCount, bool onlyDynamic, typename... BenchmarkFunctions>
    struct MatrixBenchmarkImpl {
        HeapMatrix<T> D1;
        HeapMatrix<T> D2;
        StaticHeapMatrix<T, onlyDynamic ? 1 : N, onlyDynamic ? 1 : M> S1;
        StaticHeapMatrix<T, onlyDynamic ? 1 : M, onlyDynamic ? 1 : P> S2;
        StackMatrix<T, onlyDynamic ? 1 : N, onlyDynamic ? 1 : M> M1;
        StackMatrix<T, onlyDynamic ? 1 : M, onlyDynamic ? 1 : P> M2;
        CachePaddedHeapMatrix<T> DC1;
        CachePaddedHeapMatrix<T> DC2;
        CachePaddedStaticHeapMatrix<T, onlyDynamic ? 1 : N, onlyDynamic ? 1 : M> SC1;
        CachePaddedStaticHeapMatrix<T, onlyDynamic ? 1 : M, onlyDynamic ? 1 : P> SC2;
        int functionCount;

        MatrixBenchmarkImpl() : functionCount(CalculateFunctionCount<BenchmarkFunctions...>()) {
            if constexpr(!onlyDynamic) D1 = HeapMatrix<T>(N, M);
            if constexpr(!onlyDynamic) D2 = HeapMatrix<T>(M, P);
            DC1 = CachePaddedHeapMatrix<T>(N, M);
            DC2 = CachePaddedHeapMatrix<T>(M, P);
        }

        template<typename Fun> static int CalculateFunctionCount() {
            return Fun::RunDynamic + Fun::RunStatic + Fun::RunStack;
        }
        template<typename Fun, typename Fun2, typename... Rest> static int CalculateFunctionCount() {
            return CalculateFunctionCount<Fun>() + CalculateFunctionCount<Fun2, Rest...>();
        }

        template<typename Fun> int runBenchmarks(std::ostream& out, int index = 0) {
            for (int i = 0; i < Fun::SkipCount; ++i) {
                out << ',';
            }
            int benchmarksRunCount = 0;
            std::cout << "\rStep = (" << std::setw(2) << step << "/" << std::setw(2) << endStep << ") [" << std::setw(2) << index << "/" << functionCount << "]";
            if constexpr (Fun::RunDynamic) {
                benchmarksRunCount += 1;
                if constexpr (Fun::UseCache) {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(DC1); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(DC1, DC2); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(DC1, DC2, step); });
                } else {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(D1); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(D1, D2); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(D1, D2, step); });
                }
                std::cout << "\rStep = (" << std::setw(2) << step << "/" << std::setw(2) << endStep << ") [" << std::setw(2) << index + benchmarksRunCount << "/" << functionCount << "]";
            }
            if constexpr (Fun::RunStatic) {
                benchmarksRunCount += 1;
                if constexpr (Fun::UseCache) {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(SC1); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(SC1, SC2); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(SC1, SC2, step); });
                } else {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(S1); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(S1, S2); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(S1, S2, step); });
                }
                std::cout << "\rStep = (" << std::setw(2) << step << "/" << std::setw(2) << endStep << ") [" << std::setw(2) << index + benchmarksRunCount << "/" << functionCount << "]";
            }
            if constexpr (Fun::RunStack) {
                benchmarksRunCount += 1;
                if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(M1); });
                if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(M1, M2); });
                std::cout << "\rStep = (" << std::setw(2) << step << "/" << std::setw(2) << endStep << ") [" << std::setw(2) << index + benchmarksRunCount << "/" << functionCount << "]";
            }
            return benchmarksRunCount;
        }
        template<typename Fun, typename Fun2, typename... Rest> void runBenchmarks(std::ostream& out, int index = 0) {
            int benchmarksRunCount = runBenchmarks<Fun>(out, index);
            runBenchmarks<Fun2, Rest...>(out, index + benchmarksRunCount);
        }

        void run(std::ostream& out) {
            InputMatrixCount == 3 ? (out << step) : (NInc > 0 ? (out << N) : (out << M));
            runBenchmarks<BenchmarkFunctions...>(out);
            out << std::endl;
            matrixBenchmark<T, Reps, 
                (incType == IncrementType::Add ? (N + NInc) : (N * NInc)), NInc, 
                (incType == IncrementType::Add ? (M + MInc) : (M * MInc)), MInc,
                (incType == IncrementType::Add ? (P + PInc) : (P * PInc)), PInc,
                step + 1, endStep, 
                incType, 
                InputMatrixCount, onlyDynamic, BenchmarkFunctions...
            >(out);
        }
    };

    template<typename T, int Reps, int N, int NInc, int M, int MInc, int P, int PInc, int endStep, IncrementType incType, int InputMatrixCount, bool onlyDynamic, typename... BenchmarkFunctions>
    struct MatrixBenchmarkImpl<T, Reps, N, NInc, M, MInc, P, PInc, endStep, endStep, incType, InputMatrixCount, onlyDynamic, BenchmarkFunctions...> {
        void run(std::ostream& out) {}
    };

    template<typename T, int Reps, int N, int NInc, int M, int MInc, int P, int PInc, int step, int endStep, IncrementType incType, int InputMatrixCount, bool onlyDynamic, typename... BenchmarkFunctions>
    void matrixBenchmark(std::ostream& out) {
        MatrixBenchmarkImpl<T, Reps, N, NInc, M, MInc, P, PInc, step, endStep, incType, InputMatrixCount, onlyDynamic, BenchmarkFunctions...>{}.run(out);
    }

    template<typename Fun> void printBenchmarkFunctionNames(std::ostream& outFile) {
        if (Fun::RunDynamic) outFile << "," << Fun::Name;
        if (Fun::RunStatic) outFile << "," << Fun::StaticName;
        if (Fun::RunStack) outFile << "," << Fun::StackName;
    }
    template<typename Fun, typename Fun2, typename... Rest> void printBenchmarkFunctionNames(std::ostream& outFile) {
        printBenchmarkFunctionNames<Fun>(outFile);
        printBenchmarkFunctionNames<Fun2, Rest...>(outFile);
    }
}

template<typename T, int Reps, int N, int NInc, int M, int MInc, int P, int PInc, int step, int endStep, IncrementType incType, int InputMatrixCount, bool onlyDynamic, typename... BenchmarkFunctions>
void matrixCsvBenchmark(std::string fileName, bool append=false) {
    std::ofstream outFile(fileName + ".csv", append ? std::ios_base::app : std::ios_base::out);
    if (!append) {
        InputMatrixCount == 3 ? (outFile <<  "Steps") : (outFile << "N");
        detail::printBenchmarkFunctionNames<BenchmarkFunctions...>(outFile);
        outFile << std::endl;
    }
    std::cout << '\n';
    detail::matrixBenchmark<T, Reps, N, NInc, M, MInc, P, PInc, step, endStep, incType, InputMatrixCount, onlyDynamic, BenchmarkFunctions...>(outFile);
}

#define CsvBenchmarkName(name) CSVB_ ## name

#define CreateCsvBenchmark(name, function, useCache, runDynamic, runStatic, runStack) \
struct CsvBenchmarkName(name) {\
    static constexpr bool RunDynamic = runDynamic;\
    static constexpr bool RunStatic = runStatic;\
    static constexpr bool RunStack = runStack;\
    static constexpr bool UseCache = useCache;\
    static constexpr auto Name = "D_" #name;\
    static constexpr auto StaticName = "S_" #name;\
    static constexpr auto StackName = "M_" #name;\
    static constexpr auto SkipCount = 0;\
    template<typename... Args> static void Run(const Args&... args) {\
        function(args...);\
    }\
};
#define CreateCsvBenchmarkFullFunction(name, fullFunction, useCache, runDynamic, runStatic, runStack) \
struct CsvBenchmarkName(name) {\
    static constexpr bool RunDynamic = runDynamic;\
    static constexpr bool RunStatic = runStatic;\
    static constexpr bool RunStack = runStack;\
    static constexpr bool UseCache = useCache;\
    static constexpr auto Name = "D_" #name;\
    static constexpr auto StaticName = "S_" #name;\
    static constexpr auto StackName = "M_" #name;\
    static constexpr auto SkipCount = 0;\
    template<typename... Args> static void Run(const Args&... args) {\
        fullFunction;\
    }\
};

template<int Count> struct CsvBenchmarkName(skip) {
    static constexpr bool RunDynamic = false;
    static constexpr bool RunStatic = false;
    static constexpr bool RunStack = false;
    static constexpr bool UseCache = false;
    static constexpr auto Name = "Skip";
    static constexpr auto StaticName = "Skip";
    static constexpr auto StackName = "Skip";
    static constexpr auto SkipCount = Count;
    template<typename... Args> static void Run(const Args&... args) {}
};

#define CsvBenchmarkNames_1(n1) CsvBenchmarkName(n1)
#define CsvBenchmarkNames_2(n1, n2) CsvBenchmarkNames_1(n1), CsvBenchmarkName(n2)
#define CsvBenchmarkNames_3(n1, n2, n3) CsvBenchmarkNames_2(n1, n2), CsvBenchmarkName(n3)
#define CsvBenchmarkNames_4(n1, n2, n3, n4) CsvBenchmarkNames_3(n1, n2, n3), CsvBenchmarkName(n4)
#define CsvBenchmarkNames_5(n1, n2, n3, n4, n5) CsvBenchmarkNames_4(n1, n2, n3, n4), CsvBenchmarkName(n5)
#define CsvBenchmarkNames_6(n1, n2, n3, n4, n5, n6) CsvBenchmarkNames_5(n1, n2, n3, n4, n5), CsvBenchmarkName(n6)
#define CsvBenchmarkNames_7(n1, n2, n3, n4, n5, n6, n7) CsvBenchmarkNames_6(n1, n2, n3, n4, n5, n6), CsvBenchmarkName(n7)
#define CsvBenchmarkNames_8(n1, n2, n3, n4, n5, n6, n7, n8) CsvBenchmarkNames_7(n1, n2, n3, n4, n5, n6, n7), CsvBenchmarkName(n8)
#define CsvBenchmarkNames_9(n1, n2, n3, n4, n5, n6, n7, n8, n9) CsvBenchmarkNames_8(n1, n2, n3, n4, n5, n6, n7, n8), CsvBenchmarkName(n9)
#define CsvBenchmarkNames_10(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) CsvBenchmarkNames_9(n1, n2, n3, n4, n5, n6, n7, n8, n9), CsvBenchmarkName(n10)

#endif

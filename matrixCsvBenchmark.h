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
    template<typename T, int Reps, int N, int NEnd, int Inc, IncrementType incType, int InputMatrixCount, typename... BenchmarkFunctions>
    void matrixBenchmark(std::ostream& out);

    template<typename T, int Reps, int N, int NEnd, int Inc, IncrementType incType, int InputMatrixCount, typename... BenchmarkFunctions>
    struct MatrixBenchmarkImpl {
        HeapMatrix<T> D[InputMatrixCount];
        StaticHeapMatrix<T, N, N> S[InputMatrixCount];
        CachePaddedHeapMatrix<T> DC[InputMatrixCount];
        CachePaddedStaticHeapMatrix<T, N, N> SC[InputMatrixCount];
        int functionCount;

        MatrixBenchmarkImpl() : functionCount(CalculateFunctionCount<BenchmarkFunctions...>()) {
            for (int i = 0; i < InputMatrixCount; ++i) {
                D[i] = HeapMatrix<T>(N, N);
                DC[i] = CachePaddedHeapMatrix<T>(N, N);
            }
        }

        template<typename Fun> static int CalculateFunctionCount() {
            if constexpr (Fun::RunStatic && Fun::RunDynamic) return 2;
            if constexpr (Fun::RunStatic || Fun::RunDynamic) return 1;
            return 0;
        }
        template<typename Fun, typename Fun2, typename... Rest> static int CalculateFunctionCount() {
            return CalculateFunctionCount<Fun>() + CalculateFunctionCount<Fun2, Rest...>();
        }

        template<typename Fun> int runBenchmarks(std::ostream& out, int index = 0) {
            for (int i = 0; i < Fun::SkipCount; ++i) {
                out << ',';
            }
            int benchmarksRunCount = 0;
            std::cout << "\rN = (" << std::setw(5) << N << "/" << std::setw(5) << (incType == IncrementType::Add ? (NEnd - Inc) : (NEnd / Inc)) << ") [" << std::setw(2) << index << "/" << functionCount << "]";
            if constexpr (Fun::RunDynamic) {
                benchmarksRunCount += 1;
                if constexpr (Fun::UseCache) {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(DC[0]); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(DC[0], DC[1]); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(DC[0], DC[1], DC[2]); });
                } else {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(D[0]); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(D[0], D[1]); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(D[0], D[1], D[2]); });
                }
                std::cout << "\rN = (" << std::setw(5) << N << "/" << std::setw(5) << (incType == IncrementType::Add ? (NEnd - Inc) : (NEnd / Inc)) << ") [" << std::setw(2) << index + benchmarksRunCount << "/" << functionCount << "]";
            }
            if constexpr (Fun::RunStatic) {
                benchmarksRunCount += 1;
                if constexpr (Fun::UseCache) {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(SC[0]); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(SC[0], SC[1]); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(SC[0], SC[1], SC[2]); });
                } else {
                    if constexpr (InputMatrixCount == 1) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(S[0]); });
                    if constexpr (InputMatrixCount == 2) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(S[0], S[1]); });
                    if constexpr (InputMatrixCount == 3) out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(S[0], S[1], S[2]); });
                }
                std::cout << "\rN = (" << std::setw(5) << N << "/" << std::setw(5) << (incType == IncrementType::Add ? (NEnd - Inc) : (NEnd / Inc)) << ") [" << std::setw(2) << index + benchmarksRunCount << "/" << functionCount << "]";
            }
            return benchmarksRunCount;
        }
        template<typename Fun, typename Fun2, typename... Rest> void runBenchmarks(std::ostream& out, int index = 0) {
            int benchmarksRunCount = runBenchmarks<Fun>(out, index);
            runBenchmarks<Fun2, Rest...>(out, index + benchmarksRunCount);
        }

        void run(std::ostream& out) {
            out << N;
            runBenchmarks<BenchmarkFunctions...>(out);
            out << std::endl;
            matrixBenchmark<T, Reps, (incType == IncrementType::Add ? (N + Inc) : (N * Inc)), NEnd, Inc, incType, InputMatrixCount, BenchmarkFunctions...>(out);
        }
    };

    template<typename T, int Reps, int N, int Inc, IncrementType incType, int InputMatrixCount, typename... BenchmarkFunctions>
    struct MatrixBenchmarkImpl<T, Reps, N, N, Inc, incType, InputMatrixCount, BenchmarkFunctions...> {
        void run(std::ostream& out) {}
    };

    template<typename T, int Reps, int N, int NEnd, int Inc, IncrementType incType, int InputMatrixCount, typename... BenchmarkFunctions>
    void matrixBenchmark(std::ostream& out) {
        MatrixBenchmarkImpl<T, Reps, N, NEnd, Inc, incType, InputMatrixCount, BenchmarkFunctions...>{}.run(out);
    }

    template<typename Fun> void printBenchmarkFunctionNames(std::ostream& outFile) {
        if (Fun::RunDynamic) outFile << "," << Fun::Name;
        if (Fun::RunStatic) outFile << "," << Fun::StaticName;
    }
    template<typename Fun, typename Fun2, typename... Rest> void printBenchmarkFunctionNames(std::ostream& outFile) {
        printBenchmarkFunctionNames<Fun>(outFile);
        printBenchmarkFunctionNames<Fun2, Rest...>(outFile);
    }
}

template<typename T, int Reps, int N, int NEnd, int Inc, IncrementType incType, int InputMatrixCount, typename... BenchmarkFunctions>
void matrixCsvBenchmark(std::string fileName, bool append=false) {
    std::ofstream outFile(fileName + ".csv", append ? std::ios_base::app : std::ios_base::out);
    if (!append) {
        outFile << "N";
        detail::printBenchmarkFunctionNames<BenchmarkFunctions...>(outFile);
        outFile << std::endl;
    }
    std::cout << '\n';
    detail::matrixBenchmark<T, Reps, N, NEnd, Inc, incType, InputMatrixCount, BenchmarkFunctions...>(outFile);
}

#define CsvBenchmarkName(name) CSVB_ ## name

#define CreateCsvBenchmark(name, function, useCache, runDynamic, runStatic) \
struct CsvBenchmarkName(name) {\
    static constexpr bool RunDynamic = runDynamic;\
    static constexpr bool RunStatic = runStatic;\
    static constexpr bool UseCache = useCache;\
    static constexpr auto Name = "D_" #name;\
    static constexpr auto StaticName = "S_" #name;\
    static constexpr auto SkipCount = 0;\
    template<typename... Args> static void Run(const Args&... args) {\
        function(args...);\
    }\
};
#define CreateCsvBenchmarkFullFunction(name, fullFunction, useCache, runDynamic, runStatic) \
struct CsvBenchmarkName(name) {\
    static constexpr bool RunDynamic = runDynamic;\
    static constexpr bool RunStatic = runStatic;\
    static constexpr bool UseCache = useCache;\
    static constexpr auto Name = "D_" #name;\
    static constexpr auto StaticName = "S_" #name;\
    static constexpr auto SkipCount = 0;\
    template<typename... Args> static void Run(const Args&... args) {\
        fullFunction;\
    }\
};

template<int Count> struct CsvBenchmarkName(skip) {
    static constexpr bool RunDynamic = false;
    static constexpr bool RunStatic = false;
    static constexpr bool UseCache = false;
    static constexpr auto Name = "Skip";
    static constexpr auto StaticName = "Skip";
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

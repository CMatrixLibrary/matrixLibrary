#ifndef MATRIX_CSV_BENCHMARK_H
#define MATRIX_CSV_BENCHMARK_H

#include "Matrix.h"
#include "benchmark.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace detail {
    template<typename T, int Reps, int N, int NEnd, int Inc, typename... BenchmarkFunctions>
    void matrixBenchmark(std::ostream& out);

    template<typename T, int Reps, int N, int NEnd, int Inc, typename... BenchmarkFunctions>
    struct MatrixBenchmarkImpl {
        CachePaddedHeapMatrix<T> da = CachePaddedHeapMatrix<T>(N,N);
        CachePaddedHeapMatrix<T> db = CachePaddedHeapMatrix<T>(N,N);
        CachePaddedStaticHeapMatrix<T, N, N> sa;
        CachePaddedStaticHeapMatrix<T, N, N> sb;
        int functionCount;

        MatrixBenchmarkImpl() : functionCount(CalculateFunctionCount<BenchmarkFunctions...>()) {}

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
            std::cout << "\rN = (" << std::setw(5) << N << "/" << std::setw(5) << NEnd - Inc << ") [" << std::setw(2) << index << "/" << functionCount << "]";
            if constexpr (Fun::RunDynamic) {
                benchmarksRunCount += 1;
                out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(da, db); }).getNano() / ((long long)N*(long long)N*(long long)N);
                std::cout << "\rN = (" << std::setw(5) << N << "/" << std::setw(5) << NEnd - Inc << ") [" << std::setw(2) << index + benchmarksRunCount << "/" << functionCount << "]";
            }
            if constexpr (Fun::RunStatic) {
                benchmarksRunCount += 1;
                out << ',' << benchmarkMedian(Reps, [&]() { Fun::Run(sa, sb); }).getNano() / ((long long)N*(long long)N*(long long)N);
                std::cout << "\rN = (" << std::setw(5) << N << "/" << std::setw(5) << NEnd - Inc << ") [" << std::setw(2) << index + benchmarksRunCount << "/" << functionCount << "]";
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
            matrixBenchmark<T, Reps, N + Inc, NEnd, Inc, BenchmarkFunctions...>(out);
        }
    };

    template<typename T, int Reps, int N, int Inc, typename... BenchmarkFunctions>
    struct MatrixBenchmarkImpl<T, Reps, N, N, Inc, BenchmarkFunctions...> {
        void run(std::ostream& out) {}
    };

    template<typename T, int Reps, int N, int NEnd, int Inc, typename... BenchmarkFunctions>
    void matrixBenchmark(std::ostream& out) {
        MatrixBenchmarkImpl<T, Reps, N, NEnd, Inc, BenchmarkFunctions...>{}.run(out);
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

template<typename T, int Reps, int N, int NEnd, int Inc, typename... BenchmarkFunctions>
void matrixCsvBenchmark(std::string fileName, bool append=false) {
    std::ofstream outFile(fileName + ".csv", append ? std::ios_base::app : std::ios_base::out);
    if (!append) {
        outFile << "N";
        detail::printBenchmarkFunctionNames<BenchmarkFunctions...>(outFile);
        outFile << std::endl;
    }
    std::cout << '\n';
    detail::matrixBenchmark<T, Reps, N, NEnd, Inc, BenchmarkFunctions...>(outFile);
}

#define CsvBenchmarkName(name) CSVB_ ## name

#define CreateCsvBenchmark(name, function, runDynamic, runStatic) \
struct CsvBenchmarkName(name) {\
    static constexpr bool RunDynamic = runDynamic;\
    static constexpr bool RunStatic = runStatic;\
    static constexpr auto Name = "D_" #name;\
    static constexpr auto StaticName = "S_" #name;\
    static constexpr auto SkipCount = 0;\
    template<typename... Args> static void Run(const Args&... args) {\
        function(args...);\
    }\
};
#define CreateCsvBenchmarkFullFunction(name, fullFunction, runDynamic, runStatic) \
struct CsvBenchmarkName(name) {\
    static constexpr bool RunDynamic = runDynamic;\
    static constexpr bool RunStatic = runStatic;\
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

#endif

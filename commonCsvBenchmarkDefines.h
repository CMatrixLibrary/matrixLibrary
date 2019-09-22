#ifndef COMMON_CSV_BENCHMARK_DEFINES_H
#define COMMON_CSV_BENCHMARK_DEFINES_H

#include "matrixCsvBenchmark.h"
#include "naiveMul.h"
#include "parallelMul.h"
#include "blockMul.h"
#include "parallelBlockMul.h"
#include "avxMul.h"
#include "blasMul.h"
#include "fmmUtility.h"
#include "fastMatrixMultiplyAlgorithms/genStrassen.h"

#define CreateCsvBenchmarkFmm(name, function, steps, runDynamic, runStatic) \
    CreateCsvBenchmarkFullFunction(name, [&]() {function(args..., steps); }(), runDynamic, runStatic);

#define CreateCsvBenchmarkFmm5Steps(name, function, runDynamic, runStatic) \
    CreateCsvBenchmarkFmm(name##_1, function, 1, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_2, function, 2, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_3, function, 3, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_4, function, 4, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_5, function, 5, runDynamic, runStatic);

#define CreateCsvBenchmarkFmmAllAlgorithms5Steps(name, function, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Naive, function<fmm::Naive>, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Parallel, function<fmm::Parallel>, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Block, function<fmm::Block>, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_ParallelBlock, function<fmm::ParallelBlock>, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Avx, function<fmm::Avx>, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_ParallelAvx, function<fmm::ParallelAvx>, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Blas, function<fmm::Blas>, runDynamic, runStatic)

CreateCsvBenchmark(naiveMul, naiveMul, true, true);
CreateCsvBenchmark(parallelMul, parallelMul, true, true);
CreateCsvBenchmark(blockMul, blockMul, true, true);
CreateCsvBenchmark(parallelBlockMul, parallelBlockMul, true, true);
CreateCsvBenchmark(avxMul, avx::mul, true, true);
CreateCsvBenchmark(avxParallelMul, avx::parallelMul, true, true);
CreateCsvBenchmark(blasMul, blas::mul, true, true);

CreateCsvBenchmarkFmmAllAlgorithms5Steps(strassenMinSpace, fmm::genStrassenMinSpace, true, true)
CreateCsvBenchmarkFmmAllAlgorithms5Steps(strassenLowLevel, fmm::genStrassenLowLevel, true, true)
CreateCsvBenchmarkFmmAllAlgorithms5Steps(strassenParallel, fmm::genStrassenLowLevelParallel, true, true)

#define CsvBenchmarkBaseMulAll CsvBenchmarkNames_7(naiveMul, parallelMul, blockMul, parallelBlockMul, avxMul, avxParallelMul, blasMul)
#define CsvBenchmarkBaseMulNoBlas CsvBenchmarkNames_6(naiveMul, parallelMul, blockMul, parallelBlockMul, avxMul, avxParallelMul)
#define CsvBenchmarkBaseMulNoBlasAvx CsvBenchmarkNames_4(naiveMul, parallelMul, blockMul, parallelBlockMul)

#endif

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
#include "strassen.h"

#define CreateCsvBenchmarkFmm(name, function, steps, useCache, runDynamic, runStatic) \
    CreateCsvBenchmarkFullFunction(name, [&]() {function(args..., steps); }(), useCache, runDynamic, runStatic, false);

#define CreateCsvBenchmarkFmm5Steps(name, function, useCache, runDynamic, runStatic) \
    CreateCsvBenchmark(name, function, useCache, runDynamic, runStatic, false);\
    CreateCsvBenchmarkFmm(name##_1, function, 1, useCache, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_2, function, 2, useCache, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_3, function, 3, useCache, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_4, function, 4, useCache, runDynamic, runStatic);\
    CreateCsvBenchmarkFmm(name##_5, function, 5, useCache, runDynamic, runStatic);

#define CreateCsvBenchmarkFmmAllAlgorithms5Steps(name, function, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Empty, function<fmm::Empty>, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Naive, function<fmm::Naive>, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Parallel, function<fmm::Parallel>, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Block, function<fmm::Block>, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_ParallelBlock, function<fmm::ParallelBlock>, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Avx, function<fmm::Avx>, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_ParallelAvx, function<fmm::ParallelAvx>, useCache, runDynamic, runStatic)\
    CreateCsvBenchmarkFmm5Steps(name##_Blas, function<fmm::Blas>, useCache, runDynamic, runStatic)

CreateCsvBenchmark(naiveMul, naiveMul, true, true, false, false);
CreateCsvBenchmark(parallelMul, parallelMul, true, true, false, false);
CreateCsvBenchmark(blockMul, blockMul, true, true, false, false);
CreateCsvBenchmark(parallelBlockMul, parallelBlockMul, true, true, false, false);
CreateCsvBenchmark(avxMul, avx::mul, true, true, false, false);
CreateCsvBenchmark(avxParallelMul, avx::parallelMul, true, true, false, false);
CreateCsvBenchmark(blasMul, blas::mul, true, true, false, false);

CreateCsvBenchmarkFmmAllAlgorithms5Steps(strassenHighLevel, fmm::strassenHighLevel, true, true, false)
CreateCsvBenchmarkFmmAllAlgorithms5Steps(strassenMinSpace, fmm::strassenMinSpace, true, true, false)
CreateCsvBenchmarkFmmAllAlgorithms5Steps(strassenLowLevel, fmm::strassenLowLevel, true, true, false)
CreateCsvBenchmarkFmmAllAlgorithms5Steps(strassenParallel, fmm::genStrassenLowLevelParallel, true, true, false)

#define CsvBenchmarkBaseMulAll CsvBenchmarkNames_7(naiveMul, parallelMul, blockMul, parallelBlockMul, avxMul, avxParallelMul, blasMul)
#define CsvBenchmarkBaseMulNoBlas CsvBenchmarkNames_6(naiveMul, parallelMul, blockMul, parallelBlockMul, avxMul, avxParallelMul)
#define CsvBenchmarkBaseMulNoBlasAvx CsvBenchmarkNames_4(naiveMul, parallelMul, blockMul, parallelBlockMul)

#endif

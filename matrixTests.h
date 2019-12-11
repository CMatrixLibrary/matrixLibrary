#ifndef MATRIX_TESTS_H
#define MATRIX_TESTS_H

#include "matrixCsvBenchmark.h"
#include "commonCsvBenchmarkDefines.h"
#include "fastMatrixMultiplyAlgorithms/genStrassen.h"
#include "fastMatrixMultiplyAlgorithms/genFastMul3x3.h"
#include "fastMatrixMultiplyAlgorithms/gen2x3x3.h"
#include "fastMatrixMultiplyAlgorithms/gen6x3x3.h"

CreateCsvBenchmarkFullFunction(naiveMulCache, naiveMul(args...), true, true, false, false);
CreateCsvBenchmarkFullFunction(naiveMulNoCache, naiveMul(args...), false, true, false, false);

void heapMatrixVsCachePaddedMatrix() {
    matrixCsvBenchmark<
        float,  // matrix element type
        10,     // number of executed runs to take median from (repetitions)
        64,     // starting N (matrix of size NxN)
        64,     // increment N value
        64,64,  // starting M and increment M
        64,64,  // starting P and increment P
        1,33,   // do 32 times for increasing N,M,P (so up to and including 2048) 
        IncrementType::Add, // increment type: Add means increase N by adding it with increment value
        2,      // number of matrix arguments. multiplication function requires 2 matricies
        false,
        CsvBenchmarkNames_2(naiveMulCache, naiveMulNoCache) // tested functions
    >("HeapMatrixVsCachePaddedHeapMatrix");
}

template<typename T> void baseMulTestAllDynamic() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<T, 10, 256, 256, 256, 256, 256, 256, 1, 9, IncrementType::Add, 2, false,
            CsvBenchmarkBaseMulNoBlas
        >("baseMulTestAllDynamic" + std::string(typeid(T).name()));
    }
}

template<typename T> void baseMulTestAvxDynamic() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<T, 10, 1024, 1024, 1024, 1024, 1024, 1024, 1, 9, IncrementType::Add, 2, false,
            CsvBenchmarkNames_2(avxMul, avxParallelMul)
        >("baseMulTestAvxBlasDynamic" + std::string(typeid(T).name()));
    }
}


CreateCsvBenchmarkFullFunction(avxMulCache, avx::parallelMul(args...), true, true, false, false);
CreateCsvBenchmarkFullFunction(avxMulNoCache, avx::parallelMul(args...), false, true, false, false);

void heapMatrixVsCachePaddedMatrixAvx() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<float, 10, 2048, 2048, 2048, 2048, 2048, 2048, 1, 5, IncrementType::Add, 2, false,
            CsvBenchmarkNames_2(avxMulCache, avxMulNoCache)
        >("HeapMatrixVsCachePaddedHeapMatrixAvx");
    }
}


CreateCsvBenchmark(strassenHighLevel_StaticPadding_Empty, fmm::strassenHighLevel<fmm::Empty | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenLowLevel_StaticPadding_Empty, fmm::strassenLowLevel<fmm::Empty | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenMinSpace_StaticPadding_Empty, fmm::strassenMinSpace<fmm::Empty | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenParallel_StaticPadding_Empty, fmm::strassenParallelLowLevel<fmm::Empty | fmm::StaticPadding>, true, true, false, false);

template<typename T, int N, int EndStep> void emptyStrassenByStep() {
    matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 0, EndStep, IncrementType::Add, 3, false,
        CsvBenchmarkNames_8(
            strassenHighLevel_Empty, 
            strassenHighLevel_StaticPadding_Empty,
            strassenLowLevel_Empty, 
            strassenLowLevel_StaticPadding_Empty,
            strassenMinSpace_Empty, 
            strassenMinSpace_StaticPadding_Empty,
            strassenParallel_Empty,
            strassenParallel_StaticPadding_Empty
        )
    >("emptyStrassenByStep" + std::string(typeid(T).name()) + std::to_string(N));
}
template<typename T, int N, int EndStep> void emptyStrassenByStepOnlyPeeling() {
    matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 0, EndStep, IncrementType::Add, 3, false,
        CsvBenchmarkNames_4(
            strassenHighLevel_Empty, 
            strassenLowLevel_Empty, 
            strassenMinSpace_Empty, 
            strassenParallel_Empty
        )
    >("emptyStrassenByStepOnlyPeeling" + std::string(typeid(T).name()) + std::to_string(N));
}
template<typename T, int N, int EndStep> void emptyStrassenByStepNoHigh() {
    matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 0, EndStep, IncrementType::Add, 3, false,
        CsvBenchmarkNames_6(
            strassenLowLevel_Empty, 
            strassenLowLevel_StaticPadding_Empty,
            strassenMinSpace_Empty, 
            strassenMinSpace_StaticPadding_Empty,
            strassenParallel_Empty,
            strassenParallel_StaticPadding_Empty
        )
    >("emptyStrassenByStep" + std::string(typeid(T).name()) + std::to_string(N));
}
template<typename T, int N, int EndStep> void emptyStrassenByStepNoHighOnlyPeeling() {
    matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 0, EndStep, IncrementType::Add, 3, false,
        CsvBenchmarkNames_3(
            strassenLowLevel_Empty, 
            strassenMinSpace_Empty, 
            strassenParallel_Empty
        )
    >("emptyStrassenByStepOnlyPeeling" + std::string(typeid(T).name()) + std::to_string(N));
}


CreateCsvBenchmark(strassenHighLevel_StaticPadding_Block, fmm::strassenHighLevel<fmm::Block | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenLowLevel_StaticPadding_Block, fmm::strassenLowLevel<fmm::Block | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenMinSpace_StaticPadding_Block, fmm::strassenMinSpace<fmm::Block | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenParallel_StaticPadding_Block, fmm::strassenParallelLowLevel<fmm::Block | fmm::StaticPadding>, true, true, false, false);

template<typename T, int N, int EndStep> void blockStrassenByStep() {
    matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 1, EndStep, IncrementType::Add, 3, false,
        CsvBenchmarkNames_8(
            strassenHighLevel_Block,
            strassenHighLevel_StaticPadding_Block,
            strassenLowLevel_Block,
            strassenLowLevel_StaticPadding_Block,
            strassenMinSpace_Block,
            strassenMinSpace_StaticPadding_Block,
            strassenParallel_Block,
            strassenParallel_StaticPadding_Block
        )
    >("blockStrassenByStep" + std::string(typeid(T).name()) + std::to_string(N));
}
template<typename T, int N, int EndStep> void blockStrassenByStepOnlyPeeling() {
    matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 1, EndStep, IncrementType::Add, 3, false,
        CsvBenchmarkNames_4(
            strassenHighLevel_Block,
            strassenLowLevel_Block,
            strassenMinSpace_Block,
            strassenParallel_Block
        )
    >("blockStrassenByStep" + std::string(typeid(T).name()) + std::to_string(N));
}


CreateCsvBenchmark(strassenHighLevel_StaticPadding_ParallelAvx, fmm::strassenHighLevel<fmm::ParallelAvx | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenLowLevel_StaticPadding_ParallelAvx, fmm::strassenLowLevel<fmm::ParallelAvx | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenMinSpace_StaticPadding_ParallelAvx, fmm::strassenMinSpace<fmm::ParallelAvx | fmm::StaticPadding>, true, true, false, false);
CreateCsvBenchmark(strassenParallel_StaticPadding_ParallelAvx, fmm::strassenParallelLowLevel<fmm::ParallelAvx | fmm::StaticPadding>, true, true, false, false);

template<typename T, int N, int EndStep> void parallelAvxStrassenByStep() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 1, EndStep, IncrementType::Add, 3, false,
            CsvBenchmarkNames_6(
                strassenLowLevel_ParallelAvx,
                strassenLowLevel_StaticPadding_ParallelAvx,
                strassenMinSpace_ParallelAvx,
                strassenMinSpace_StaticPadding_ParallelAvx,
                strassenParallel_ParallelAvx,
                strassenParallel_StaticPadding_ParallelAvx
            )
        >("ParallelAvxStrassenByStep" + std::string(typeid(T).name()) + std::to_string(N));
    }
}
template<typename T, int N, int EndStep> void parallelAvxStrassenByStepOnlyPeeling() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<T, 10, N, 0, N, 0, N, 0, 1, EndStep, IncrementType::Add, 3, false,
            CsvBenchmarkNames_3(
                strassenLowLevel_ParallelAvx,
                strassenMinSpace_ParallelAvx,
                strassenParallel_ParallelAvx
            )
        >("ParallelAvxStrassenByStepOnlyPeeling" + std::string(typeid(T).name()) + std::to_string(N));
    }
}

CreateCsvBenchmarkFmm5Steps(strassen_ParallelBlock, fmm::strassenMinSpace<fmm::ParallelBlock>, true, true, false);
CreateCsvBenchmarkFmm5Steps(gen2x3x3_ParallelBlock, fmm::gen2x3x3MinSpace<fmm::ParallelBlock>, true, true, false);
CreateCsvBenchmarkFmm5Steps(gen6x3x3_ParallelBlock, fmm::gen6x3x3MinSpace<fmm::ParallelBlock>, true, true, false);
CreateCsvBenchmarkFmm5Steps(gen3x3x3_ParallelBlock, fmm::genFastMul3x3MinSpace<fmm::ParallelBlock>, true, true, false);

template<typename T> void differentFmm1() {
    matrixCsvBenchmark<T, 10, 108, 108, 108, 108, 108, 108, 0, 18, IncrementType::Add, 2, false,
        CSVB_parallelBlockMul,
        CSVB_strassen_ParallelBlock_1,
        CSVB_strassen_ParallelBlock_2,
        CSVB_strassen_ParallelBlock_3,
        CSVB_gen2x3x3_ParallelBlock_1,
        CSVB_gen2x3x3_ParallelBlock_2,
        CSVB_gen2x3x3_ParallelBlock_3,
        CSVB_gen6x3x3_ParallelBlock_1,
        CSVB_gen6x3x3_ParallelBlock_2,
        CSVB_gen3x3x3_ParallelBlock_1,
        CSVB_gen3x3x3_ParallelBlock_2
    >("differentFmm1" + std::string(typeid(T).name()));
}

template<typename T> void differentFmm2() {
    matrixCsvBenchmark<T, 10, 216, 216, 1296, 0, 1296, 0, 0, 18, IncrementType::Add, 2, false,
        CSVB_parallelBlockMul,
        CSVB_strassen_ParallelBlock_1,
        CSVB_strassen_ParallelBlock_2,
        CSVB_strassen_ParallelBlock_3,
        CSVB_gen2x3x3_ParallelBlock_1,
        CSVB_gen2x3x3_ParallelBlock_2,
        CSVB_gen2x3x3_ParallelBlock_3,
        CSVB_gen6x3x3_ParallelBlock_1,
        CSVB_gen6x3x3_ParallelBlock_2,
        CSVB_gen6x3x3_ParallelBlock_3,
        CSVB_gen3x3x3_ParallelBlock_1,
        CSVB_gen3x3x3_ParallelBlock_2,
        CSVB_gen3x3x3_ParallelBlock_3
    >("differentFmm2" + std::string(typeid(T).name()));
}

template<typename T, int endStep> void strassenPivotPoints() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<T, 10, 512, 512, 512, 512, 512, 512, 0, endStep, IncrementType::Add, 2, true,
            CSVB_avxParallelMul,
            CSVB_strassenMinSpace_ParallelAvx_1,
            CSVB_strassenMinSpace_ParallelAvx_2,
            CSVB_strassenMinSpace_ParallelAvx_3,
            CSVB_strassenMinSpace_ParallelAvx_4
        >("strassenPivotPoints" + std::string(typeid(T).name()));
    }
}


template<typename M> auto add100Times(const MatrixInterface<M>& m) {
    return m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m + m+m+m+m+m+m+m+m+m+m;
}
CreateCsvBenchmark(add100, add100Times, true, true, true, true);
template<typename T> void smallMatrix100Adds() {
    matrixCsvBenchmark<T, 10000, 1, 1, 1, 1, 1, 1, 0, 25, IncrementType::Add, 1, false,
        CsvBenchmarkName(add100)
    >("smallMatrix100Adds" + std::string(typeid(T).name()));
}

template<typename T> void bigMatrix100Adds() {
    matrixCsvBenchmark<T, 10, 64, 64, 64, 64, 64, 64, 0, 20, IncrementType::Add, 1, false,
        CsvBenchmarkName(add100)
    >("bigMatrix100Adds" + std::string(typeid(T).name()));
}


CreateCsvBenchmark(inverse, inverse, true, true, false, false);
CreateCsvBenchmark(blockInverse, blockInverse, true, true, false, false);
CreateCsvBenchmark(fastUnstableBlockInverse, fastUnstableBlockInverse, true, true, false, false);
template<typename T> void smallInverseTest() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<T, 100, 16, 16, 16, 16, 16, 16, 0, 20, IncrementType::Add, 1, false,
            CsvBenchmarkNames_3(
                inverse,
                blockInverse,
                fastUnstableBlockInverse
            )
        >("smallInverseTest" + std::string(typeid(T).name()));
    }
}


template<typename T> void bigInverseTest() {
    if constexpr (avx::IsAvailable) {
        matrixCsvBenchmark<T, 10, 256, 256, 256, 256, 256, 256, 0, 8, IncrementType::Add, 1, false,
            CsvBenchmarkNames_3(
                inverse,
                blockInverse,
                fastUnstableBlockInverse
            )
        >("bigInverseTest" + std::string(typeid(T).name()));
    }
}

#endif

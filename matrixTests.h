#ifndef MATRIX_TESTS_H
#define MATRIX_TESTS_H

#include "matrixCsvBenchmark.h"
#include "commonCsvBenchmarkDefines.h"

CreateCsvBenchmarkFullFunction(naiveMulCache, naiveMul(args...), true, true, false);
CreateCsvBenchmarkFullFunction(naiveMulNoCache, naiveMul(args...), false, true, false);

void heapMatrixVsCachePaddedMatrix() {
    matrixCsvBenchmark<
        double, // matrix element type
        3,      // number of executed runs to take median from (repetitions)
        8,      // starting N (matrix of size NxN)
        4096,   // ending N (not included, so up to 2048)
        2,      // increment value
        IncrementType::Mul, // increment type: Mul means increase N by multiplying it with increment value
        TimeFormat::Absolute, // Absolute -> show time in seconds
        2,      // number of matrix arguments. multiplication function requires 2 matricies
        CsvBenchmarkNames_2(naiveMulCache, naiveMulNoCache) // tested functions
    >("HeapMatrixVsCachePaddedHeapMatrix");
}

template<typename T> void baseMulTestAllDynamic() {
    matrixCsvBenchmark<T, 3, 64, 4096, 2, IncrementType::Mul, TimeFormat::DivByNCubed, 2,
        CsvBenchmarkBaseMulAll
    >("baseMulTestAllDynamic");
}

template<typename T> void baseMulTestAvxBlasDynamic() {
    matrixCsvBenchmark<T, 3, 256, 32768, 2, IncrementType::Mul, TimeFormat::DivByNCubed, 2,
        CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)
    >("baseMulTestAvxBlasDynamic");
}


/*
    old ones (not correct)
*/
template<typename T> void baseMulTest_all(std::string fileName="baseMulTest") {
    matrixCsvBenchmark<T, 20, 65, 129, 16, CsvBenchmarkBaseMulAll>(fileName);
    matrixCsvBenchmark<T, 10, 129, 257, 32, CsvBenchmarkBaseMulAll>(fileName, true);
    matrixCsvBenchmark<T, 4, 257, 513, 64, CsvBenchmarkBaseMulAll>(fileName, true);
    matrixCsvBenchmark<T, 4, 513, 1025, 128, CsvBenchmarkBaseMulAll>(fileName, true);
    matrixCsvBenchmark<T, 2, 1025, 2049, 256, CsvBenchmarkBaseMulAll>(fileName, true);
    matrixCsvBenchmark<T, 2, 2049, 3073, 512, CsvBenchmarkBaseMulAll>(fileName, true);
    matrixCsvBenchmark<T, 2, 3073, 4097, 512, CsvBenchmarkNames_4(skip<8>, avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 2, 4097, 8193, 1024, CsvBenchmarkNames_4(skip<8>, avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 1, 8193, 13313, 1024, CsvBenchmarkNames_3(skip<10>, avxParallelMul, blasMul)>(fileName, true);
}
template<typename T> void baseMulTest_noBlas(std::string fileName="baseMulTest") {
    matrixCsvBenchmark<T, 20, 65, 129, 16, CsvBenchmarkBaseMulNoBlas>(fileName);
    matrixCsvBenchmark<T, 10, 129, 257, 32, CsvBenchmarkBaseMulNoBlas>(fileName, true);
    matrixCsvBenchmark<T, 4, 257, 513, 64, CsvBenchmarkBaseMulNoBlas>(fileName, true);
    matrixCsvBenchmark<T, 4, 513, 1025, 128, CsvBenchmarkBaseMulNoBlas>(fileName, true);
    matrixCsvBenchmark<T, 2, 1025, 2049, 256, CsvBenchmarkBaseMulNoBlas>(fileName, true);
    matrixCsvBenchmark<T, 2, 2049, 3073, 512, CsvBenchmarkBaseMulNoBlas>(fileName, true);
    matrixCsvBenchmark<T, 2, 3073, 4097, 512, CsvBenchmarkNames_3(skip<8>, avxMul, avxParallelMul)>(fileName, true);
    matrixCsvBenchmark<T, 2, 4097, 9217, 1024, CsvBenchmarkNames_3(skip<8>, avxMul, avxParallelMul)>(fileName, true);
}
template<typename T> void baseMulTest_noBlasAvx(std::string fileName="baseMulTest") {
    matrixCsvBenchmark<T, 10, 65, 129, 16, CsvBenchmarkBaseMulNoBlasAvx>(fileName);
    matrixCsvBenchmark<T, 6, 129, 257, 32, CsvBenchmarkBaseMulNoBlasAvx>(fileName, true);
    matrixCsvBenchmark<T, 2, 257, 513, 64, CsvBenchmarkBaseMulNoBlasAvx>(fileName, true);
    matrixCsvBenchmark<T, 2, 513, 1153, 128, CsvBenchmarkBaseMulNoBlasAvx>(fileName, true);
}
template<typename T> void baseMulTest_avxAndBlasOnly(std::string fileName="baseMulTest") {
    matrixCsvBenchmark<T, 20, 65, 129, 16, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName);
    matrixCsvBenchmark<T, 10, 129, 257, 32, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 4, 257, 513, 64, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 4, 513, 1025, 128, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 2, 1025, 2049, 256, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 2, 2049, 3073, 512, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 2, 3073, 4097, 512, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 2, 4097, 8193, 1024, CsvBenchmarkNames_3(avxMul, avxParallelMul, blasMul)>(fileName, true);
    matrixCsvBenchmark<T, 1, 8193, 13313, 1024, CsvBenchmarkNames_3(skip<2>, avxParallelMul, blasMul)>(fileName, true);
}

#endif

#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H

#include "MatrixInterface.h"
#include "Matrix.h"
#include "MatrixView.h"

#include "genericArithmeticOperations.h"
#include "MatrixExtendedFunctions.h"
#include "matrixOperators.h"
#include "naiveBasicOperations.h"

#include "fmmUtility.h"
#include "strassen.h"
#include "fastMatrixMultiplyAlgorithms/gen2x3x3.h"
#include "fastMatrixMultiplyAlgorithms/gen6x3x3.h"
#include "fastMatrixMultiplyAlgorithms/genFastMul3x3.h"
#include "fastMatrixMultiplyAlgorithms/genStrassen.h"

#include "avxSimd.h"
#include "avxMul.h"
#include "blasMul.h"
#include "naiveMul.h"
#include "parallelMul.h"
#include "blockMul.h"
#include "parallelBlockMul.h"

#include "inverse.h"

#include "ArrayConstStaticView.h"
#include "ArrayConstView.h"
#include "ArrayStaticView.h"
#include "ArrayView.h"

#include "alignedAllocation.h"
#include "StackAllocator.h"

#include "compilerMacros.h"

#include "Range.h"
#include "RangeZip.h"

#include "benchmark.h"
#include "matrixCsvBenchmark.h"
#include "matrixTests.h"
#include "ThreadPool.h"

#endif

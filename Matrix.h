#ifndef MATRIX_H
#define MATRIX_H
#include "StackMatrix.h"
#include "CachePaddedHeapMatrix.h"
#include "HeapMatrix.h"
#include "StaticHeapMatrix.h"
#include "CachePaddedStaticHeapMatrix.h"
#include <type_traits>

template<typename T, mtl::size_t Rows = 0, mtl::size_t Columns = 0>
using Matrix = 
    typename std::conditional <(Rows != 0 && Columns != 0 && Rows * Columns <= 1024),
        StackMatrix<T, Rows, Columns>,
    typename std::conditional<(Rows == 0 && Columns == 0),
        HeapMatrix<T>,
    // else
        StaticHeapMatrix<T, Rows, Columns>
    >::type
>::type;

#endif

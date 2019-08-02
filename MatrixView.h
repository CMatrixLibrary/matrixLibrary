#pragma once
#include "MatrixDynamicView.h"
#include "MatrixStaticSizeView.h"
#include "MatrixStaticEffectiveColumnView.h"
#include "MatrixStaticView.h"

#include "MatrixConstDynamicView.h"
#include "MatrixConstStaticSizeView.h"
#include "MatrixConstStaticEffectiveColumnView.h"
#include "MatrixConstStaticView.h"

#include <type_traits>

template<typename T, mtl::size_t Rows = 0, mtl::size_t Columns = 0, mtl::size_t EffColumns = 0>
using MatrixView =
    typename std::conditional <(Rows != 0 && Columns != 0 && EffColumns != 0),
        MatrixStaticView<T, Rows, Columns, EffColumns>,
    typename std::conditional<(Rows != 0 && Columns != 0),
        MatrixStaticSizeView<T, Rows, Columns>,
    typename std::conditional<(Rows != 0),
        MatrixStaticEffectiveColumnView<T, Rows>,
    // else
        MatrixDynamicView<T>
    >::type>::type
>::type;

template<typename T, mtl::size_t Rows = 0, mtl::size_t Columns = 0, mtl::size_t EffColumns = 0>
using MatrixConstView =
    typename std::conditional <(Rows != 0 && Columns != 0 && EffColumns != 0),
        MatrixConstStaticView<T, Rows, Columns, EffColumns>,
    typename std::conditional<(Rows != 0 && Columns != 0),
        MatrixConstStaticSizeView<T, Rows, Columns>,
    typename std::conditional<(Rows != 0),
        MatrixConstStaticEffectiveColumnView<T, Rows>,
    // else
        MatrixConstDynamicView<T>
    >::type>::type
>::type;
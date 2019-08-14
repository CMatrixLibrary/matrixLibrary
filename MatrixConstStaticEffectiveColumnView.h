#ifndef MATRIX_CONST_STATIC_EFFECTIVE_COLUMN_VIEW_H
#define MATRIX_CONST_STATIC_EFFECTIVE_COLUMN_VIEW_H
#include "MatrixInterface.h"
#include "HeapMatrix.h"
#include "StackMatrix.h"
#include "StaticHeapMatrix.h"

template<typename T, mtl::size_t EffCol> class MatrixConstStaticEffectiveColumnView : public MatrixInterface<MatrixConstStaticEffectiveColumnView<T, EffCol>> {
    using Interface = MatrixInterface<MatrixConstStaticEffectiveColumnView<T, EffCol>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = Interface::DynamicValue;
    static const mtl::size_t Cols = Interface::DynamicValue;
    static const mtl::size_t EffCols = EffCol;

private:
    const T* data_;
    mtl::size_t rowCount_;
    mtl::size_t columnCount_;

public:
    MatrixConstStaticEffectiveColumnView() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0)
    {}
    template<typename MT> MatrixConstStaticEffectiveColumnView(const MatrixInterface<MT>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount())
    {}
    template<typename MT> MatrixConstStaticEffectiveColumnView(const MatrixInterface<MT>& matrix, mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount) :
        data_(matrix.data() + startColumn + startRow * EffCol),
        rowCount_(rowCount),
        columnCount_(columnCount)
    {}
    MatrixConstStaticEffectiveColumnView(HeapMatrix<T>&&) = delete;
    template<int RCount, int CCount> MatrixConstStaticEffectiveColumnView(StaticHeapMatrix<T, RCount, CCount>&&) = delete;
    template<int RCount, int CCount> MatrixConstStaticEffectiveColumnView(StackMatrix<T, RCount, CCount>&&) = delete;
};
#endif

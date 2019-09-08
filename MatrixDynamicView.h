#ifndef MATRIX_DYNAMIC_VIEW_H
#define MATRIX_DYNAMIC_VIEW_H
#include "MatrixInterface.h"

template<typename T> class HeapMatrix;

template<typename T> class MatrixDynamicView : public MatrixInterface<MatrixDynamicView<T>> {
    using Interface = MatrixInterface<MatrixDynamicView<T>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = Interface::DynamicValue;
    static const mtl::size_t Cols = Interface::DynamicValue;
    static const mtl::size_t EffCols = Interface::DynamicValue;

private:
    T* data_;
    mtl::size_t rowCount_;
    mtl::size_t columnCount_;
    mtl::size_t effectiveColumnCount_;

public:
    MatrixDynamicView() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0),
        effectiveColumnCount_(0)
    {}
    template<typename MT> MatrixDynamicView(MatrixInterface<MT>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    template<typename MT> MatrixDynamicView(MatrixInterface<MT>& matrix, mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount) :
        data_(matrix.data() + startColumn + startRow * matrix.effectiveColumnCount()),
        rowCount_(rowCount),
        columnCount_(columnCount),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    MatrixDynamicView(T* data, mtl::size_t rowCount, mtl::size_t columnCount, mtl::size_t effectiveColumnCount) :
        data_(data),
        rowCount_(rowCount),
        columnCount_(columnCount),
        effectiveColumnCount_(effectiveColumnCount)
    {}

    mtl::size_t _effectiveColumnCount() const { return effectiveColumnCount_; }
};
#endif

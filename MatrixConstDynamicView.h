#pragma once
#include "MatrixInterface.h"
#include "HeapMatrix.h"
#include "StackMatrix.h"
#include "StaticHeapMatrix.h"

template<typename T> class MatrixConstDynamicView : public MatrixInterface<MatrixConstDynamicView<T>> {
    using Interface = MatrixInterface<MatrixConstDynamicView<T>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = Interface::DynamicValue;
    static const mtl::size_t Cols = Interface::DynamicValue;
    static const mtl::size_t EffCols = Interface::DynamicValue;

private:
    const T* data_;
    mtl::size_t rowCount_;
    mtl::size_t columnCount_;
    mtl::size_t effectiveColumnCount_;

public:
    MatrixConstDynamicView() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0),
        effectiveColumnCount_(0)
    {}
    template<typename MT> MatrixConstDynamicView(const MatrixInterface<MT>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    template<typename MT> MatrixConstDynamicView(const MatrixInterface<MT>& matrix, mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount) :
        data_(matrix.data() + startColumn + startRow * matrix.effectiveColumnCount()),
        rowCount_(rowCount),
        columnCount_(columnCount),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    MatrixConstDynamicView(HeapMatrix<T>&&) = delete;
    template<int RCount, int CCount> MatrixConstDynamicView(StaticHeapMatrix<T, RCount, CCount>&&) = delete;
    template<int RCount, int CCount> MatrixConstDynamicView(StackMatrix<T, RCount, CCount>&&) = delete;

    mtl::size_t _effectiveColumnCount() const { return effectiveColumnCount_; }
};
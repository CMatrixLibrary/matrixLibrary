#pragma once
#include "MatrixInterface.h"
#include "HeapMatrix.h"
#include "StackMatrix.h"
#include "StaticHeapMatrix.h"

template<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_> class MatrixConstStaticSizeView : public MatrixInterface<MatrixConstStaticSizeView<T, RowCount_, ColumnCount_>> {
    using Interface = MatrixInterface<MatrixConstStaticSizeView<T, RowCount_, ColumnCount_>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = RowCount_;
    static const mtl::size_t Cols = ColumnCount_;
    static const mtl::size_t EffCols = Interface::DynamicValue;

private:
    const T* data_;
    mtl::size_t effectiveColumnCount_;

public:
    MatrixConstStaticSizeView() :
        data_(nullptr),
        effectiveColumnCount_(0)
    {}
    template<typename MT> MatrixConstStaticSizeView(const MatrixInterface<MT>& matrix, mtl::size_t startRow=0, mtl::size_t startColumn=0) :
        data_(matrix.data() + startColumn + startRow * matrix.effectiveColumnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    MatrixConstStaticSizeView(HeapMatrix<T>&&) = delete;
    template<int RCount, int CCount> MatrixConstStaticSizeView(StaticHeapMatrix<T, RCount, CCount>&&) = delete;
    template<int RCount, int CCount> MatrixConstStaticSizeView(StackMatrix<T, RCount, CCount>&&) = delete;

    template<typename MT> MatrixConstStaticSizeView& operator=(const MatrixInterface<MT>& matrix) {
        data_ = matrix.data();
        effectiveColumnCount_ = matrix.effectiveColumnCount();
        return *this;
    }
    template<int RCount, int CCount> MatrixConstStaticSizeView& operator=(StaticHeapMatrix<T, RCount, CCount>&&) = delete;
    template<int RCount, int CCount> MatrixConstStaticSizeView& operator=(StackMatrix<T, RCount, CCount>&&) = delete;
    
    mtl::size_t _effectiveColumnCount() const { return effectiveColumnCount_; }
};
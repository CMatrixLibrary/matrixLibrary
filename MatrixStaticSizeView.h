#ifndef MATRIX_STATIC_SIZE_VIEW_H
#define MATRIX_STATIC_SIZE_VIEW_H
#include "MatrixInterface.h"

template<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_> class MatrixStaticSizeView : public MatrixInterface<MatrixStaticSizeView<T, RowCount_, ColumnCount_>> {
    using Interface = MatrixInterface<MatrixStaticSizeView<T, RowCount_, ColumnCount_>>;
    friend Interface;
    
public:
    using ValueType = T;

    static const mtl::size_t Rows = RowCount_;
    static const mtl::size_t Cols = ColumnCount_;
    static const mtl::size_t EffCols = Interface::DynamicValue;

private:
    T* data_;
    mtl::size_t effectiveColumnCount_;

public:
    MatrixStaticSizeView() :
        data_(nullptr),
        effectiveColumnCount_(0)
    {}
    template<typename MT> MatrixStaticSizeView(MatrixInterface<MT>& matrix, mtl::size_t startRow=0, mtl::size_t startColumn=0) :
        data_(matrix.data() + startColumn + startRow * matrix.effectiveColumnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    MatrixStaticSizeView(T* data, mtl::size_t effectiveColumnCount) :
        data_(data),
        effectiveColumnCount_(effectiveColumnCount)
    {}

    template<typename MT> MatrixStaticSizeView& operator=(MatrixInterface<MT>& matrix) {
        data_ = matrix.data();
        effectiveColumnCount_ = matrix.effectiveColumnCount();
        return *this;
    }
    template<int A> MatrixStaticSizeView& operator=(MatrixStaticView<T, RowCount_, ColumnCount_, A> matrix) {
        data_ = matrix.data();
        effectiveColumnCount_ = matrix.effectiveColumnCount();
        return *this;
    }

    mtl::size_t _effectiveColumnCount() const { return effectiveColumnCount_; }
};
#endif

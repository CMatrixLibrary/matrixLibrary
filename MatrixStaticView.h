#ifndef MATRIX_STATIC_VIEW_H
#define MATRIX_STATIC_VIEW_H
#include "MatrixInterface.h"

template<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_, mtl::size_t EffCols_> class MatrixStaticView : public MatrixInterface<MatrixStaticView<T, RowCount_, ColumnCount_, EffCols_>> {
    using Interface = MatrixInterface<MatrixStaticView<T, RowCount_, ColumnCount_, EffCols_>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = RowCount_;
    static const mtl::size_t Cols = ColumnCount_;
    static const mtl::size_t EffCols = EffCols_;

private:
    T* data_;

public:
    MatrixStaticView() :
        data_(nullptr)
    {}
    template<typename MT> MatrixStaticView(MatrixInterface<MT>& matrix, mtl::size_t startRow=0, mtl::size_t startColumn=0) :
        data_(matrix.data() + startColumn + startRow * EffCols_)
    {}
    MatrixStaticView(T* data) :
        data_(data)
    {}
};
#endif

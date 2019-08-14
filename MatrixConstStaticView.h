#ifndef MATRIX_CONST_STATIC_VIEW_H
#define MATRIX_CONST_STATIC_VIEW_H
#include "MatrixInterface.h"
#include "HeapMatrix.h"
#include "StackMatrix.h"
#include "StaticHeapMatrix.h"

template<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_, mtl::size_t EffCols_> class MatrixConstStaticView : public MatrixInterface<MatrixConstStaticView<T, RowCount_, ColumnCount_, EffCols_>> {
    using Interface = MatrixInterface<MatrixConstStaticView<T, RowCount_, ColumnCount_, EffCols_>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = RowCount_;
    static const mtl::size_t Cols = ColumnCount_;
    static const mtl::size_t EffCols = EffCols_;

private:
    const T* data_;

public:
    MatrixConstStaticView() :
        data_(nullptr)
    {}
    template<typename MT> MatrixConstStaticView(const MatrixInterface<MT>& matrix, mtl::size_t startRow=0, mtl::size_t startColumn=0) :
        data_(matrix.data() + startColumn + startRow * EffCols_)
    {}
    MatrixConstStaticView(HeapMatrix<T>&&) = delete;
    template<int RCount, int CCount> MatrixConstStaticView(StaticHeapMatrix<T, RCount, CCount>&&) = delete;
    template<int RCount, int CCount> MatrixConstStaticView(StackMatrix<T, RCount, CCount>&&) = delete;
};
#endif

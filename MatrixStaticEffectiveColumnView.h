#ifndef MATRIX_STATIC_EFFECTIVE_COLUMN_VIEW_H
#define MATRIX_STATIC_EFFECTIVE_COLUMN_VIEW_H
#include "MatrixInterface.h"

template<typename T, mtl::size_t EffCol> class MatrixStaticEffectiveColumnView : public MatrixInterface<MatrixStaticEffectiveColumnView<T, EffCol>> {
    using Interface = MatrixInterface<MatrixStaticEffectiveColumnView<T, EffCol>>;
    friend Interface;
    
public:
    using ValueType = T;

    static const mtl::size_t Rows = Interface::DynamicValue;
    static const mtl::size_t Cols = Interface::DynamicValue;
    static const mtl::size_t EffCols = EffCol;

private:
    T* data_;
    mtl::size_t rowCount_;
    mtl::size_t columnCount_;

public:
    MatrixStaticEffectiveColumnView() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0)
    {}
    template<typename MT> MatrixStaticEffectiveColumnView(MatrixInterface<MT>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount())
    {}
    template<typename MT> MatrixStaticEffectiveColumnView(MatrixInterface<MT>& matrix, mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount) :
        data_(matrix.data() + startColumn + startRow * EffCol),
        rowCount_(rowCount),
        columnCount_(columnCount)
    {}
    MatrixStaticEffectiveColumnView(T* data, mtl::size_t rowCount, mtl::size_t columnCount) :
        data_(data),
        rowCount_(rowCount),
        columnCount_(columnCount)
    {}
};
#endif

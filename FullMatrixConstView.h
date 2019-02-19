#pragma once
#include "FullMatrix.h"
#include "FullMatrixView.h"
#include "VectorConstView.h"
#include "debugAssert.h"

template<typename T> class FullMatrixConstView {
public:
    FullMatrixConstView() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0),
        effectiveColumnCount_(0)
    {}
    FullMatrixConstView(const FullMatrix<T>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.columnCount())
    {}
    FullMatrixConstView(FullMatrix<T>&&) = delete;
    FullMatrixConstView(const FullMatrixView<T>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    FullMatrixConstView(const FullMatrixConstView<T>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    FullMatrixConstView(const FullMatrix<T>& matrix, int startRow, int startColumn, int rowCount, int columnCount) :
        data_(matrix.data() + startColumn + startRow * matrix.columnCount()),
        rowCount_(rowCount),
        columnCount_(columnCount),
        effectiveColumnCount_(matrix.columnCount())
    {
        debugAssertOp(startRow, < , matrix.rowCount());
        debugAssertOp(startColumn, < , matrix.columnCount());
        debugAssert(rowCount <= matrix.rowCount() - startRow, rowCount, " <= ", matrix.rowCount(), " - ", startRow);
        debugAssert(columnCount <= matrix.columnCount() - startColumn, columnCount, " <= ", matrix.columnCount(), " - ", startColumn);
    }
    FullMatrixConstView(const FullMatrixView<T>& matrix, int startRow, int startColumn, int rowCount, int columnCount) :
        data_(matrix.data() + startColumn + startRow * matrix.effectiveColumnCount()),
        rowCount_(rowCount),
        columnCount_(columnCount),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {
        debugAssertOp(startRow, < , matrix.rowCount());
        debugAssertOp(startColumn, < , matrix.columnCount());
        debugAssert(rowCount <= matrix.rowCount() - startRow, rowCount, " <= ", matrix.rowCount(), " - ", startRow);
        debugAssert(columnCount <= matrix.columnCount() - startColumn, columnCount, " <= ", matrix.columnCount(), " - ", startColumn);
    }
    FullMatrixConstView(const FullMatrixConstView<T>& matrix, int startRow, int startColumn, int rowCount, int columnCount) :
        data_(matrix.data() + startColumn + startRow * matrix.effectiveColumnCount()),
        rowCount_(rowCount),
        columnCount_(columnCount),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {
        debugAssertOp(startRow, < , matrix.rowCount());
        debugAssertOp(startColumn, < , matrix.columnCount());
        debugAssert(rowCount <= matrix.rowCount() - startRow, rowCount, " <= ", matrix.rowCount(), " - ", startRow);
        debugAssert(columnCount <= matrix.columnCount() - startColumn, columnCount, " <= ", matrix.columnCount(), " - ", startColumn);
    }

    const T& at(int row, int column) const {
        debugAssertOp(row, < , rowCount_);
        debugAssertOp(column, < , columnCount_);
        return data_[column + row * effectiveColumnCount_];
    }

    VectorConstView<T> operator[](int i) const {
        debugAssertOp(i, < , rowCount_);
        return VectorConstView<T>(data_ + i * effectiveColumnCount_, columnCount_);
    }

    int rowCount() const {
        return rowCount_;
    }
    int columnCount() const {
        return columnCount_;
    }
    int effectiveColumnCount() const {
        return effectiveColumnCount_;
    }

    const T* data() const {
        return data_;
    }

    int size() const {
        return rowCount_ * columnCount_;
    }

    FullMatrixConstView<T> subMatrixView(int startRow, int startCol, int rowsCount, int columnsCount) const {
        return FullMatrixConstView<T>(*this, startRow, startCol, rowsCount, columnsCount);
    }

private:
    const T* data_;
    int rowCount_;
    int columnCount_;
    int effectiveColumnCount_;
};

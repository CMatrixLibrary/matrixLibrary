#pragma once
#include "FullMatrix.h"

template<typename T> class FullMatrixView {
public:
    FullMatrixView() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0),
        effectiveColumnCount_(0)
    {}
    FullMatrixView(FullMatrix<T>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.columnCount())
    {}
    FullMatrixView(FullMatrixView<T>& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    FullMatrixView(FullMatrixView<T>&& matrix) :
        data_(matrix.data()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount()),
        effectiveColumnCount_(matrix.effectiveColumnCount())
    {}
    FullMatrixView(FullMatrix<T>& matrix, int startRow, int startColumn, int rowCount, int columnCount) :
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
    FullMatrixView(FullMatrixView<T>& matrix, int startRow, int startColumn, int rowCount, int columnCount) :
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
    FullMatrixView<T>& operator=(FullMatrixView<T>& matrix) {
        data_ = matrix.data();
        rowCount_ = matrix.rowCount();
        columnCount_ = matrix.columnCount();
        effectiveColumnCount_ = matrix.effectiveColumnCount();
        return *this;
    }
    FullMatrixView<T>& operator=(FullMatrixView<T>&& matrix) {
        data_ = matrix.data();
        rowCount_ = matrix.rowCount();
        columnCount_ = matrix.columnCount();
        effectiveColumnCount_ = matrix.effectiveColumnCount();
        return *this;
    }

    T& at(int column, int row) {
        debugAssertOp(column, < , columnCount_);
        debugAssertOp(row, < , rowCount_);
        return data_[column + row * effectiveColumnCount_];
    }
    const T& at(int column, int row) const {
        debugAssertOp(column, < , columnCount_);
        debugAssertOp(row, < , rowCount_);
        return data_[column + row * effectiveColumnCount_];
    }

    T& operator[](int ind) {
        debugAssertOp(ind, <, size());
        return data_[ind + (ind / columnCount_) * (effectiveColumnCount_ - columnCount_)];
    }
    const T& operator[](int ind) const {
        debugAssertOp(ind, <, size());
        return data_[ind + (ind / columnCount_) * (effectiveColumnCount_ - columnCount_)];
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

    T* data() {
        return data_;
    }
    const T* data() const {
        return data_;
    }

    int size() const {
        return rowCount_ * columnCount_;
    }

    FullMatrixView<T> subMatrixView(int startRow, int startCol, int rowsCount, int columnsCount) {
        return FullMatrixView<T>(*this, startRow, startCol, rowsCount, columnsCount);
    }
    FullMatrixConstView<T> subMatrixView(int startRow, int startCol, int rowsCount, int columnsCount) const {
        return FullMatrixConstView<T>(*this, startRow, startCol, rowsCount, columnsCount);
    }

    template<template<typename> typename Matrix> void copy(const Matrix<T>& matrix) {
        debugAssertOp(rowCount_, >=, matrix.rowCount());
        debugAssertOp(columnCount_, >=, matrix.columnCount());
        for (int row = 0; row < matrix.rowCount(); ++row) {
            for (int column = 0; column < matrix.columnCount(); ++column) {
                at(column, row) = matrix.at(column, row);
            }
        }
    }

private:
    T* data_;
    int rowCount_;
    int columnCount_;
    int effectiveColumnCount_;
};

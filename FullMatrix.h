#pragma once

template<typename T> class FullMatrixView;
template<typename T> class FullMatrixConstView;

template<typename T> class FullMatrix {
public:
    FullMatrix() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0)
    {}
    FullMatrix(int rowCount, int columnCount) :
        data_(new T[rowCount * columnCount]()),
        rowCount_(rowCount),
        columnCount_(columnCount)
    {}
    FullMatrix(const FullMatrix<T>& matrix) :
        data_(new T[matrix.rowCount() * matrix.columnCount()]()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount())
    {
        copy(matrix);
    }
    FullMatrix(const FullMatrixView<T>& matrix) :
        data_(new T[matrix.rowCount() * matrix.columnCount()]()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount())
    {
        copy(matrix);
    }
    FullMatrix(const FullMatrixConstView<T>& matrix) :
        data_(new T[matrix.rowCount() * matrix.columnCount()]()),
        rowCount_(matrix.rowCount()),
        columnCount_(matrix.columnCount())
    {
        copy(matrix);
    }
    ~FullMatrix() {
        delete[] data_;
    }

    template<template<typename> typename Matrix> FullMatrix& operator=(const Matrix<T>& matrix) {
        rowCount_ = matrix.rowCount();
        columnCount_ = matrix.columnCount();

        delete[] data_;
        data_ = new T[rowCount_ * columnCount_]();

        copy(matrix);

        return *this;
    }
    FullMatrix& operator=(FullMatrix<T> matrix) {
        std::swap(rowCount_, matrix.rowCount_);
        std::swap(columnCount_, matrix.columnCount_);
        std::swap(data_, matrix.data_);
        return *this;
    }

    T& at(int column, int row) {
        return data_[column + row * columnCount_];
    }
    const T& at(int column, int row) const {
        return data_[column + row * columnCount_];
    }

    T& operator[](int i) {
        return data_[i];
    }
    const T& operator[](int i) const {
        return data_[i];
    }

    int rowCount() const {
        return rowCount_;
    }
    int columnCount() const {
        return columnCount_;
    }

    const T* data() const {
        return data_;
    }
    T* data() {
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
        for (int row = 0; row < matrix.rowCount(); ++row) {
            for (int column = 0; column < matrix.columnCount(); ++column) {
                at(column, row) = matrix.at(column, row);
            }
        }
    }

    void shrink(int newRowCount, int newColumnCount) {
        auto newData = new T[newRowCount * newColumnCount]();
        for (int row = 0; row < newRowCount; ++row) {
            for (int column = 0; column < newColumnCount; ++column) {
                newData[column + row * newColumnCount] = at(column, row);
            }
        }
        delete[] data_;
        data_ = newData;
        rowCount_ = newRowCount;
        columnCount_ = newColumnCount;
    }

private:
    T* data_;
    int rowCount_;
    int columnCount_;
};

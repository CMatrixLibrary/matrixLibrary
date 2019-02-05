#pragma once
#include <vector>
#include "utilityDetails.h"
#include "naiveOperations.h"

template<typename T> class FullMatrixView;

template<typename T> class FullMatrix {
public:
    FullMatrix(int rowsCount, int colsCount) :
        rowsCount_(rowsCount),
        columnsCount_(colsCount),
        data_(rowsCount * colsCount)
    {}

    T& at(int x, int y) {
        return data_[x + y * columnsCount_];
    }
    const T& at(int x, int y) const {
        return data_[x + y * columnsCount_];
    }

    T& operator[](int i) {
        return data_[i];
    }
    const T& operator[](int i) const {
        return data_[i];
    }

    void insert(const FullMatrix<T>& a) {
        for (int x = 0; x < a.columnsCount(); ++x) {
            for (int y = 0; y < a.rowsCount(); ++y) {
                at(x, y) = a.at(x, y);
            }
        }
    }
    void insert(FullMatrixView<T> a) {
        for (int x = 0; x < a.columnsCount(); ++x) {
            for (int y = 0; y < a.rowsCount(); ++y) {
                at(x, y) = a.at(x, y);
            }
        }
    }
    void shrink(int newRowsCount, int newColumnsCount) {
        std::vector<T> newData(newRowsCount * newColumnsCount);
        for (int x = 0; x < newColumnsCount; ++x) {
            for (int y = 0; y < newRowsCount; ++y) {
                newData[x + y * newColumnsCount] = at(x, y);
            }
        }
        data_ = std::move(newData);
        rowsCount_ = newRowsCount;
        columnsCount_ = newColumnsCount;
    }
    
    int rowsCount() const {
        return rowsCount_;
    }
    int columnsCount() const {
        return columnsCount_;
    }
    const std::vector<T>& data() const {
        return data_;
    }
    std::vector<T>& data() {
        return data_;
    }

private:
    int rowsCount_;
    int columnsCount_;
    std::vector<T> data_;
};
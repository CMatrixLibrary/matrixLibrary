#pragma once
#include "FullMatrix.h"

template<typename T> class FullSubMatrix {
public:
    FullSubMatrix(FullMatrix<T>& ref) :
        ref(ref),
        startRow_(0),
        startCol_(0),
        rowsCount_(ref.rowsCount()),
        columnsCount_(ref.columnsCount())
    {}
    FullSubMatrix(FullMatrix<T>& ref, int startRow, int startCol, int rowsCount, int colsCount) :
        ref(ref),
        startRow_(startRow),
        startCol_(startCol),
        rowsCount_(rowsCount),
        columnsCount_(colsCount)
    {}
    FullSubMatrix(FullSubMatrix<T> ref, int startRow, int startCol, int rowsCount, int colsCount) :
        ref(ref.ref),
        startRow_(ref.startRow_ + startRow),
        startCol_(ref.startCol_ + startCol),
        rowsCount_(rowsCount),
        columnsCount_(colsCount)
    {}

    T& operator[](int i) {
        return ref.data()[startCol_ + startRow_ * ref.columnsCount() + i + (i / columnsCount_) * (ref.columnsCount() - columnsCount_)];
    }
    const T& operator[](int i) const {
        return ref.data()[startCol_ + startRow_ * ref.columnsCount() + i + (i / columnsCount_) * (ref.columnsCount() - columnsCount_)];
    }

    T& at(int x, int y) {
        return ref.data()[startCol_ + startRow_ * ref.columnsCount() + x + y * ref.columnsCount()];
    }
    const T& at(int x, int y) const {
        return ref.data()[startCol_ + startRow_ * ref.columnsCount() + x + y * ref.columnsCount()];
    }

    int rowsCount() const {
        return rowsCount_;
    }
    int columnsCount() const {
        return columnsCount_;
    }
    int size() {
        return rowsCount_ * columnsCount_;
    }
    int startRow() const {
        return startRow_;
    }
    int startCol() const {
        return startCol_;
    }

private:
    FullMatrix<T>& ref;
    int startRow_;
    int startCol_;
    int rowsCount_;
    int columnsCount_;
};
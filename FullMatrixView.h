#pragma once
#include "FullMatrix.h"
#include "FullSubMatrix.h"

template<typename T> class FullMatrixView {
public:
    FullMatrixView(const FullMatrix<T>& ref) :
        ref(ref),
        startRow_(0),
        startCol_(0),
        rowsCount_(ref.rowsCount()),
        columnsCount_(ref.columnsCount())
    {}
    FullMatrixView(const FullMatrix<T>& ref, int startRow, int startCol, int rowsCount, int colsCount) :
        ref(ref),
        startRow_(startRow),
        startCol_(startCol),
        rowsCount_(rowsCount),
        columnsCount_(colsCount)
    {}
    FullMatrixView(const FullSubMatrix<T>& ref) :
        ref(ref.ref),
        startRow_(ref.startRow()),
        startCol_(ref.startCol()),
        rowsCount_(ref.rowsCount()),
        columnsCount_(ref.columnsCount())
    {}
    FullMatrixView(const FullMatrixView<T>& ref, int startRow, int startCol, int rowsCount, int colsCount) :
        ref(ref.ref),
        startRow_(ref.startRow_ + startRow),
        startCol_(ref.startCol_ + startCol),
        rowsCount_(rowsCount),
        columnsCount_(colsCount)
    {}

    const T& operator[](int i) const {
        return ref.data()[startCol_ + startRow_ * ref.columnsCount() + i + (i / columnsCount_) * (ref.columnsCount() - columnsCount_)];
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

private:
    const FullMatrix<T>& ref;
    int startRow_;
    int startCol_;
    int rowsCount_;
    int columnsCount_;
};
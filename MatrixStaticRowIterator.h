#ifndef MATRIX_STATIC_ROW_ITERATOR_H
#define MATRIX_STATIC_ROW_ITERATOR_H
#include "ArrayStaticView.h"

template<typename T, mtl::size_t ColumnCount, mtl::size_t EffectiveColumnCount> class MatrixStaticRowIterator {
public:
    MatrixStaticRowIterator(T* data) :
        view_(data)
    {}

    ArrayStaticView<T, ColumnCount>& operator*() {
        return view_;
    }

    MatrixStaticRowIterator& operator++() {
        view_ = ArrayStaticView<T, ColumnCount>(view_.data() + EffectiveColumnCount);
        return *this;
    }
    MatrixStaticRowIterator& operator--() {
        view_ = ArrayStaticView<T, ColumnCount>(view_.data() - EffectiveColumnCount);
        return *this;
    }

    bool operator==(const MatrixStaticRowIterator other) const {
        return view_.data() == other.view_.data();
    }
    bool operator!=(const MatrixStaticRowIterator other) const {
        return view_.data() != other.view_.data();
    }

private:
    ArrayStaticView<T, ColumnCount> view_;
};

#endif

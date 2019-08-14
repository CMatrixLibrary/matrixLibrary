#ifndef MATRIX_CONST_STATIC_ROW_ITERATOR_H
#define MATRIX_CONST_STATIC_ROW_ITERATOR_H
#include "ArrayConstStaticView.h"

template<typename T, mtl::size_t ColumnCount, mtl::size_t EffectiveColumnCount> class MatrixConstStaticRowIterator {
public:
    MatrixConstStaticRowIterator(const T* data) :
        view_(data)
    {}

    ArrayConstStaticView<T, ColumnCount>& operator*() {
        return view_;
    }

    MatrixConstStaticRowIterator& operator++() {
        view_ = ArrayConstStaticView<T, ColumnCount>(view_.data() + EffectiveColumnCount);
        return *this;
    }
    MatrixConstStaticRowIterator& operator--() {
        view_ = ArrayConstStaticView<T, ColumnCount>(view_.data() - EffectiveColumnCount);
        return *this;
    }

    bool operator==(const MatrixConstStaticRowIterator other) const {
        return view_.data() == other.view_.data();
    }
    bool operator!=(const MatrixConstStaticRowIterator other) const {
        return view_.data() != other.view_.data();
    }

private:
    ArrayConstStaticView<T, ColumnCount> view_;
};

#endif

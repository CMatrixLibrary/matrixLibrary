#pragma once
#include "ArrayView.h"

template<typename T> class MatrixRowIterator {
public:
    MatrixRowIterator(T* data, mtl::size_t columnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(columnCount)
    {}
    MatrixRowIterator(T* data, mtl::size_t columnCount, mtl::size_t effectiveColumnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(effectiveColumnCount)
    {}

    ArrayView<T>& operator*() {
        return view_;
    }

    MatrixRowIterator& operator++() {
        view_ = ArrayView<T>(view_.data() + effectiveColumnCount_, view_.size());
        return *this;
    }
    MatrixRowIterator& operator--() {
        view_ = ArrayView<T>(view_.data() - effectiveColumnCount_, view_.size());
        return *this;
    }

    bool operator==(const MatrixRowIterator& other) const {
        return view_.data() == other.view_.data()
            && view_.size() == other.view_.size()
            && effectiveColumnCount_ == other.effectiveColumnCount_;
    }
    bool operator!=(const MatrixRowIterator& other) const {
        return view_.data() != other.view_.data()
            || view_.size() != other.view_.size()
            || effectiveColumnCount_ != other.effectiveColumnCount_;
    }

private:
    ArrayView<T> view_;
    mtl::size_t effectiveColumnCount_;
};

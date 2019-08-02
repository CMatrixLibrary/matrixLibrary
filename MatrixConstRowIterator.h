#pragma once
#include "ArrayConstView.h"

template<typename T> class MatrixConstRowIterator {
public:
    MatrixConstRowIterator(const T* data, mtl::size_t columnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(columnCount)
    {}
    MatrixConstRowIterator(const T* data, mtl::size_t columnCount, mtl::size_t effectiveColumnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(effectiveColumnCount)
    {}

    ArrayConstView<T>& operator*() {
        return view_;
    }

    MatrixConstRowIterator& operator++() {
        view_ = ArrayConstView<T>(view_.data() + effectiveColumnCount_, view_.size());
        return *this;
    }
    MatrixConstRowIterator& operator--() {
        view_ = ArrayConstView<T>(view_.data() - effectiveColumnCount_, view_.size());
        return *this;
    }

    bool operator==(const MatrixConstRowIterator<T>& other) const {
        return view_.data() == other.view_.data()
            && view_.size() == other.view_.size()
            && effectiveColumnCount_ == other.effectiveColumnCount_;
    }
    bool operator!=(const MatrixConstRowIterator<T>& other) const {
        return view_.data() != other.view_.data()
            || view_.size() != other.view_.size()
            || effectiveColumnCount_ != other.effectiveColumnCount_;
    }

private:
    ArrayConstView<T> view_;
    mtl::size_t effectiveColumnCount_;
};

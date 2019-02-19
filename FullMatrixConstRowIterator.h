#pragma once
#include "VectorConstView.h"

template<typename T> class FullMatrixConstRowIterator {
public:
    FullMatrixConstRowIterator(const T* data, int columnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(columnCount)
    {}
    FullMatrixConstRowIterator(const T* data, int columnCount, int effectiveColumnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(effectiveColumnCount)
    {}

    VectorConstView<T>& operator*() {
        return view_;
    }

    FullMatrixConstRowIterator<T>& operator++() {
        view_ = VectorConstView<T>(view_.data() + effectiveColumnCount_, view_.size());
        return *this;
    }
    FullMatrixConstRowIterator<T>& operator--() {
        view_ = VectorConstView<T>(view_.data() - effectiveColumnCount_, view_.size());
        return *this;
    }

    bool operator==(const FullMatrixConstRowIterator<T>& other) const {
        return view_.data() == other.view_.data()
            && view_.size() == other.view_.size()
            && effectiveColumnCount_ == other.effectiveColumnCount_;
    }
    bool operator!=(const FullMatrixConstRowIterator<T>& other) const {
        return view_.data() != other.view_.data()
            || view_.size() != other.view_.size()
            || effectiveColumnCount_ != other.effectiveColumnCount_;
    }

private:
    VectorConstView<T> view_;
    int effectiveColumnCount_;
};

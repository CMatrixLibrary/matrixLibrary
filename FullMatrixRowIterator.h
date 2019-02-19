#pragma once
#include "VectorView.h"

template<typename T> class FullMatrixRowIterator {
public:
    FullMatrixRowIterator(T* data, int columnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(columnCount)
    {}
    FullMatrixRowIterator(T* data, int columnCount, int effectiveColumnCount) :
        view_(data, columnCount),
        effectiveColumnCount_(effectiveColumnCount)
    {}

    VectorView<T>& operator*() {
        return view_;
    }

    FullMatrixRowIterator<T>& operator++() {
        view_ = VectorView<T>(view_.data() + effectiveColumnCount_, view_.size());
        return *this;
    }
    FullMatrixRowIterator<T>& operator--() {
        view_ = VectorView<T>(view_.data() - effectiveColumnCount_, view_.size());
        return *this;
    }

    bool operator==(const FullMatrixRowIterator<T>& other) const {
        return view_.data() == other.view_.data()
            && view_.size() == other.view_.size()
            && effectiveColumnCount_ == other.effectiveColumnCount_;
    }
    bool operator!=(const FullMatrixRowIterator<T>& other) const {
        return view_.data() != other.view_.data()
            || view_.size() != other.view_.size()
            || effectiveColumnCount_ != other.effectiveColumnCount_;
    }

private:
    VectorView<T> view_;
    int effectiveColumnCount_;
};

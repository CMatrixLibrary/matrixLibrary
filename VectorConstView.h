#pragma once
#include "VectorView.h"

template<typename T> class VectorConstView {
public:
    VectorConstView(const VectorView<T>& view) :
        data_(view.data()),
        size_(view.size())
    {}
    VectorConstView(const T* data, int size) :
        data_(data),
        size_(size)
    {}
    VectorConstView(const T* data, int startIndex, int endIndex) :
        data_(data + startIndex),
        size_(endIndex - startIndex + 1)
    {}

    const T& operator[](int index) const {
        return data_[index];
    }

    const T* data() const {
        return data_;
    }

    int size() const {
        return size_;
    }

    const T* begin() const {
        return data_;
    }
    const T* end() const {
        return data_ + size_;
    }

private:
    const T* data_;
    int size_;
};
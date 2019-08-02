#pragma once
#include "types.h"
#include "debugAssert.h"
#include "ArrayView.h"

template<typename T> class ArrayConstView {
public:
    ArrayConstView(const ArrayView<T>& view) :
        data_(view.data()),
        size_(view.size())
    {}
    ArrayConstView(const T* data, mtl::size_t size) :
        data_(data),
        size_(size)
    {}
    ArrayConstView(const T* data, mtl::size_t startIndex, mtl::size_t endIndex) :
        data_(data + startIndex),
        size_(endIndex - startIndex + 1)
    {}

    const T& operator[](mtl::size_t index) const {
        debugAssertOp(index, < , size());
        return data_[index];
    }

    const T* data() const {
        return data_;
    }

    mtl::size_t size() const {
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
    mtl::size_t size_;
};
#pragma once
#include "types.h"
#include "debugAssert.h"

template<typename T> class ArrayView {
public:
    ArrayView(T* data, mtl::size_t size) :
        data_(data),
        size_(size)
    {}
    ArrayView(T* data, mtl::size_t startIndex, mtl::size_t endIndex) :
        data_(data + startIndex),
        size_(endIndex - startIndex + 1)
    {}

    T& operator[](mtl::size_t index) {
        debugAssertOp(index, < , size());
        return data_[index];
    }
    const T& operator[](mtl::size_t index) const {
        debugAssertOp(index, < , size());
        return data_[index];
    }

    T* data() {
        return data_;
    }
    const T* data() const {
        return data_;
    }

    mtl::size_t size() const {
        return size_;
    }

    T* begin() {
        return data_;
    }
    T* end() {
        return data_ + size_;
    }
    const T* begin() const {
        return data_;
    }
    const T* end() const {
        return data_ + size_;
    }

private:
    T* data_;
    mtl::size_t size_;
};
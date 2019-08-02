#pragma once
#include "types.h"
#include "debugAssert.h"

template<typename T, mtl::size_t Size> class ArrayStaticView {
public:
    ArrayStaticView(T* data) :
        data_(data)
    {}
    ArrayStaticView(T* data, mtl::size_t startIndex) :
        data_(data + startIndex)
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
        return Size;
    }

    T* begin() {
        return data_;
    }
    T* end() {
        return data_ + Size;
    }
    const T* begin() const {
        return data_;
    }
    const T* end() const {
        return data_ + Size;
    }

private:
    T* data_;
};
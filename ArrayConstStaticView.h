#pragma once
#include "types.h"
#include "debugAssert.h"
#include "ArrayStaticView.h"

template<typename T, mtl::size_t Size> class ArrayConstStaticView {
public:
    ArrayConstStaticView(const ArrayStaticView<T, Size> view) :
        data_(view.data())
    {}
    ArrayConstStaticView(const T* data) :
        data_(data)
    {}
    ArrayConstStaticView(const T* data, mtl::size_t startIndex) :
        data_(data + startIndex)
    {}

    const T& operator[](mtl::size_t index) const {
        debugAssertOp(index, < , size());
        return data_[index];
    }

    const T* data() const {
        return data_;
    }

    mtl::size_t size() const {
        return Size;
    }

    const T* begin() const {
        return data_;
    }
    const T* end() const {
        return data_ + Size;
    }

private:
    const T* data_;
};
#pragma once

template<typename T> class VectorView {
public:
    VectorView(T* data, int size) :
        data_(data),
        size_(size)
    {}
    VectorView(T* data, int startIndex, int endIndex) :
        data_(data + startIndex),
        size_(endIndex - startIndex + 1)
    {}

    T& operator[](int index) {
        return data_[index];
    }
    const T& operator[](int index) const {
        return data_[index];
    }

    T* data() {
        return data_;
    }
    const T* data() const {
        return data_;
    }

    int size() const {
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
    int size_;
};
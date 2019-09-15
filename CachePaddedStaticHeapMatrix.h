#ifndef CACHE_PADDED_STATIC_HEAP_MATRIX_H
#define CACHE_PADDED_STATIC_HEAP_MATRIX_H
#include "MatrixInterface.h"
#include <type_traits>

template<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_> 
class CachePaddedStaticHeapMatrix : public MatrixInterface<CachePaddedStaticHeapMatrix<T, RowCount_, ColumnCount_>> {
    using Interface = MatrixInterface<CachePaddedStaticHeapMatrix<T, RowCount_, ColumnCount_>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = RowCount_;
    static const mtl::size_t Cols = ColumnCount_;
    static const mtl::size_t EffCols = ColumnCount_ % 64 != 0 ? ColumnCount_ : ColumnCount_ + 1;

private:
    T* data_;

public:
    CachePaddedStaticHeapMatrix() {
        arrayNew();
    }
    template<typename MT> CachePaddedStaticHeapMatrix(const MatrixInterface<MT>& m) {
        arrayNew();
        this->copy(m);
    }
    CachePaddedStaticHeapMatrix(const CachePaddedStaticHeapMatrix& m) {
        arrayNew();
        this->copy(m);
    }
    CachePaddedStaticHeapMatrix(CachePaddedStaticHeapMatrix&& m) :
        data_(m.data_)
    {
        m.data_ = nullptr;
    }
    ~CachePaddedStaticHeapMatrix() {
        arrayDelete();
    }

    template<typename MT> CachePaddedStaticHeapMatrix& operator=(const MatrixInterface<MT>& m) {
        static_assert(MT::ConstexprRowCount() == Interface::DynamicValue || Rows == MT::ConstexprRowCount(), "mismatch RowCount");
        static_assert(MT::ConstexprColumnCount() == Interface::DynamicValue || Cols == MT::ConstexprColumnCount(), "mismatch ColumnCount");
        
        if (reinterpret_cast<const CachePaddedStaticHeapMatrix*>(&m) == this) {
            return *this;
        }

        arrayDelete();

        arrayNew();
        this->copy(m);

        return *this;
    }
    CachePaddedStaticHeapMatrix& operator=(CachePaddedStaticHeapMatrix&& m) {
        if (&m == this) {
            return *this;
        }

        arrayDelete();

        data_ = m.data_;
        m.data_ = nullptr;

        return *this;
    }

private:
    void arrayNew() {
        data_ = static_cast<T*>(std::malloc(Rows * EffCols * sizeof(T)));
        if constexpr (std::is_class_v<T>) {
            for (std::size_t row = 0; row < Rows; ++row) {
                for (std::size_t column = 0; column < Cols; ++column) {
                    new(&data_[column + row* EffCols]) T();
                }
            }
        }
    }
    void arrayDelete() {
        if constexpr (std::is_class_v<T>) {
            for (std::size_t row = 0; row < Rows; ++row) {
                for (std::size_t column = 0; column < Cols; ++column) {
                    data_[column + row * EffCols].~T();
                }
            }
        }
        std::free(data_);
    }
};
#endif

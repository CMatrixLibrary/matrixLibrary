#pragma once
#include "MatrixInterface.h"
#include <type_traits>

template<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_> 
class StaticHeapMatrix : public MatrixInterface<StaticHeapMatrix<T, RowCount_, ColumnCount_>> {
    using Interface = MatrixInterface<StaticHeapMatrix<T, RowCount_, ColumnCount_>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = RowCount_;
    static const mtl::size_t Cols = ColumnCount_;
    static const mtl::size_t EffCols = ColumnCount_;

private:
    T* data_;

public:
    StaticHeapMatrix() :
        data_(new T[Rows * Cols])
    {}
    explicit StaticHeapMatrix(T* newedPtr) : 
        data_(newedPtr)
    {}
    template<typename MT> StaticHeapMatrix(const MatrixInterface<MT>& m) :
        data_(new T[Rows * Cols])
    {
        this->copy(m);
    }
    StaticHeapMatrix(const StaticHeapMatrix& m) :
        data_(new T[Rows * Cols])
    {
        this->copy(m);
    }
    StaticHeapMatrix(StaticHeapMatrix&& m) :
        data_(m.data_)
    {
        m.data_ = nullptr;
    }
    ~StaticHeapMatrix() {
        delete[] data_;
    }

    template<typename MT> StaticHeapMatrix& operator=(const MatrixInterface<MT>& m) {
        static_assert(MT::ConstexprRowCount() == Interface::DynamicValue || Rows == MT::ConstexprRowCount(), "mismatch RowCount");
        static_assert(MT::ConstexprColumnCount() == Interface::DynamicValue || Cols == MT::ConstexprColumnCount(), "mismatch ColumnCount");
        
        if (reinterpret_cast<const StaticHeapMatrix*>(&m) == this) {
            return *this;
        }

        delete[] data_;

        data_ = new T[m.rowCount() * m.columnCount()];
        this->copy(m);

        return *this;
    }
    StaticHeapMatrix& operator=(StaticHeapMatrix&& m) {
        if (&m == this) {
            return *this;
        }

        delete[] data_;

        data_ = m.data_;
        m.data_ = nullptr;

        return *this;
    }
};
#pragma once
#include "MatrixInterface.h"
#include <array>

template<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_> 
class StackMatrix : public MatrixInterface<StackMatrix<T, RowCount_, ColumnCount_>> {
    using Interface = MatrixInterface<StackMatrix<T, RowCount_, ColumnCount_>>;
    friend class Iterface;

public:
    static const mtl::size_t Rows = RowCount_;
    static const mtl::size_t Cols = ColumnCount_;
    static const mtl::size_t EffCols = ColumnCount_;
    using ValueType = T;

private:
    static constexpr mtl::size_t Size = Rows * Cols;
    std::array<T, Size> data_;

public:
    StackMatrix() : data_(std::array<T, Size>{}) {}
    template<typename MT> StackMatrix(const MatrixInterface<MT>& m) {
        this->copy(m);
    }
    StackMatrix(const StackMatrix& m) {
        this->copy(m);
    }
    StackMatrix(StackMatrix&& m) :
        data_(m.data_)
    {}

    template<typename MT> StackMatrix& assign(const MatrixInterface<MT>& m) {
        static_assert(m.ConstexprRowCount() == Interface::DynamicValue || Rows == m.ConstexprRowCount(), "mismatch RowCount");
        static_assert(m.ConstexprColumnCount() == Interface::DynamicValue || Cols == m.ConstexprColumnCount(), "mismatch ColumnCount");
        
        if (reinterpret_cast<const StackMatrix*>(&m) == this) {
            return *this;
        }
        
        this->copy(m);

        return *this;
    }
    template<typename MT> StackMatrix& operator=(const MatrixInterface<MT>& m) {
        return assign(m);
    }
    StackMatrix& operator=(const StackMatrix& m) {
        return assign(m);
    }
    StackMatrix& operator=(StackMatrix&& m) {
        if (&m == this) {
            return *this;
        }
        
        this->copy(m);

        return *this;
    }

    T* _data()               { return data_.data(); }
    const T* _data()   const { return data_.data(); }
};
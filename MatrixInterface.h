#pragma once

#include "detect.h"
#include "types.h"
#include "debugAssert.h"
#include "ArrayView.h"
#include "ArrayConstView.h"
#include "MatrixRowIterator.h"
#include "MatrixConstRowIterator.h"
#include "MatrixStaticRowIterator.h"
#include "MatrixConstStaticRowIterator.h"

template<typename T> class HeapMatrix;
template<typename T, mtl::size_t, mtl::size_t> class StaticHeapMatrix;
template<typename T, mtl::size_t, mtl::size_t> class StackMatrix;
template<typename T> class MatrixDynamicView;
template<typename T> class MatrixConstDynamicView;
template<typename T, mtl::size_t> class MatrixStaticEffectiveColumnView;
template<typename T, mtl::size_t> class MatrixConstStaticEffectiveColumnView;
template<typename T, mtl::size_t, mtl::size_t> class MatrixStaticSizeView;
template<typename T, mtl::size_t, mtl::size_t> class MatrixConstStaticSizeView;
template<typename T, mtl::size_t, mtl::size_t, mtl::size_t> class MatrixStaticView;
template<typename T, mtl::size_t, mtl::size_t, mtl::size_t> class MatrixConstStaticView;

template<typename T> using createNew_t = decltype(std::declval<const T>()._createNew());
template<typename T> using effectiveColumnCount_t = decltype(std::declval<const T>()._effectiveColumnCount());
template<typename T> using columnCount_t = decltype(std::declval<const T>()._columnCount());
template<typename T> using rowCount_t = decltype(std::declval<const T>()._rowCount());
template<typename T> using data_t = decltype(std::declval<T>()._data());
template<typename T> using constData_t = decltype(std::declval<const T>()._data());
template<typename T> using subMatrixView_t = decltype(std::declval<T>()._subMatrixView(0,0,0,0));
template<typename T> using subMatrixConstView_t = decltype(std::declval<const T>()._subMatrixView(0,0,0,0));
template<typename T> using subMatrixStaticView_t = decltype(std::declval<T>()._subMatrixView<0, 0>(0, 0));
template<typename T> using subMatrixStaticConstView_t = decltype(std::declval<const T>()._subMatrixView<0, 0>(0, 0));

template<typename T> using indexOperator_t = decltype(std::declval<T>()._index(0));
template<typename T> using constindexOperator_t = decltype(std::declval<const T>()._index(0));

template<typename T> using begin_t = decltype(std::declval<T>()._begin());
template<typename T> using end_t = decltype(std::declval<T>()._end());
template<typename T> using constbegin_t = decltype(std::declval<const T>()._begin());
template<typename T> using constend_t = decltype(std::declval<const T>()._end());

template<typename MatrixType> class MatrixInterface {
public:
    static constexpr mtl::size_t ConstexprRowCount() {
        return MatrixType::Rows;
    }
    static constexpr mtl::size_t ConstexprColumnCount() {
        return MatrixType::Cols;
    }
    static constexpr mtl::size_t ConstexprEffectiveColumnCount() {
        return MatrixType::EffCols;
    }
    
    static constexpr mtl::size_t CRow() {
        return ConstexprRowCount();
    }
    static constexpr mtl::size_t CCol() {
        return ConstexprColumnCount();
    }
    static constexpr mtl::size_t CEffCol() {
        return ConstexprEffectiveColumnCount();
    }

    static constexpr mtl::size_t DynamicValue = 0;

    static constexpr bool HasConstexprRowAndColumnCount() {
        return ConstexprColumnCount() != DynamicValue && ConstexprRowCount() != DynamicValue;
    }

    auto createNew() const {
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, createNew_t>{}) {
            return it._createNew();
        } else if constexpr (HasConstexprRowAndColumnCount()) {
            return StaticHeapMatrix<typename MatrixType::ValueType, ConstexprRowCount(), ConstexprColumnCount()>{};
        } else {
            return HeapMatrix<typename MatrixType::ValueType>(rowCount(), columnCount());
        }
    }
    auto createNew(int rowCount, int columnCount) const {
        return HeapMatrix<typename MatrixType::ValueType>(rowCount, columnCount);
    }
    template<mtl::size_t RowCount, mtl::size_t ColumnCount> auto  createNew() const {
        return StaticHeapMatrix<typename MatrixType::ValueType, RowCount, ColumnCount>{};
    }

    auto effectiveColumnCount() const {
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, effectiveColumnCount_t>{}) {
            return it._effectiveColumnCount();
        } else if constexpr (ConstexprEffectiveColumnCount() != DynamicValue) {
            return ConstexprEffectiveColumnCount();
        } else {
            return columnCount();
        }
    }
    auto columnCount() const {
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, columnCount_t>{}) {
            return it._columnCount();
        } else if constexpr (ConstexprColumnCount() != DynamicValue) {
            return ConstexprColumnCount();
        } else {
            return it.columnCount_;
        }
    }
    auto rowCount() const {
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, rowCount_t>{}) {
            return it._rowCount();
        } else if constexpr (ConstexprRowCount() != DynamicValue) {
            return ConstexprRowCount();
        } else {
            return it.rowCount_;
        }
    }
    auto data() { 
        auto& it = static_cast<MatrixType&>(*this);
        if constexpr (detect<MatrixType, data_t>{}) {
            return it._data();
        } else {
            return it.data_;
        }
    }
    auto data() const {
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, constData_t>{}) {
            return it._data();
        } else {
            return it.data_;
        }
    }
    
    auto operator[](mtl::size_t i) {
        debugAssertOp(i, < , rowCount());
        auto& it = static_cast<MatrixType&>(*this);
        if constexpr (detect<MatrixType, indexOperator_t>{}) {
            return it._index(i);
        } else if constexpr (ConstexprColumnCount() != DynamicValue) {
            return ArrayStaticView<typename MatrixType::ValueType, ConstexprColumnCount()>(it.data() + i * it.effectiveColumnCount());
        } else {
            return ArrayView<typename MatrixType::ValueType>(it.data() + i * it.effectiveColumnCount(), it.columnCount());
        }
    }
    auto operator[](mtl::size_t i) const {
        debugAssertOp(i, < , rowCount());
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, constindexOperator_t>{}) {
            return it._index(i);
        } else if constexpr (ConstexprColumnCount() != DynamicValue) {
            return ArrayConstStaticView<typename MatrixType::ValueType, ConstexprColumnCount()>(it.data() + i * it.effectiveColumnCount());
        } else {
            return ArrayConstView<typename MatrixType::ValueType>(it.data() + i * it.effectiveColumnCount(), it.columnCount());
        }
    }
    auto& at(mtl::size_t row, mtl::size_t column) {
        debugAssertOp(row, < , rowCount());
        debugAssertOp(column, < , columnCount());
        auto& it = static_cast<MatrixType&>(*this);
        return it.data()[column + row * it.effectiveColumnCount()];
    }
    const auto& at(mtl::size_t row, mtl::size_t column) const {
        debugAssertOp(row, < , rowCount());
        debugAssertOp(column, < , columnCount());
        auto& it = static_cast<const MatrixType&>(*this);
        return it.data()[column + row * it.effectiveColumnCount()];
    }
    auto size() const {
        auto& it = static_cast<const MatrixType&>(*this);
        return it.rowCount() * it.columnCount();
    }
    template<typename MatrixType2> void copy(const MatrixInterface<MatrixType2>& other) {
        debugAssertOp(rowCount(), >= , other.rowCount());
        debugAssertOp(columnCount(), >= , other.columnCount());
        auto& it = static_cast<MatrixType&>(*this);
        for (mtl::size_t row = 0; row < other.rowCount(); ++row) {
            for (mtl::size_t column = 0; column < other.columnCount(); ++column) {
                it.at(row, column) = other.at(row, column);
            }
        }
    }
    
    auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount) {
        debugAssertOp(startRow, < , this->rowCount());
        debugAssertOp(startColumn, < , this->columnCount());
        debugAssertOp(rowCount, <= , this->rowCount() - startRow);
        debugAssertOp(columnCount, <= , this->columnCount() - startColumn);
        if constexpr (detect<MatrixType, subMatrixView_t>{}) {
            return static_cast<MatrixType*>(this)->_subMatrixView(startRow, startColumn, rowCount, columnCount);
        } else if constexpr(CEffCol() != DynamicValue) {
            return MatrixStaticEffectiveColumnView<typename MatrixType::ValueType, CEffCol()>(*this, startRow, startColumn, rowCount, columnCount);
        } else {
            return MatrixDynamicView<typename MatrixType::ValueType>(*this, startRow, startColumn, rowCount, columnCount);
        }
    }
    auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount) const {
        debugAssertOp(startRow, < , this->rowCount());
        debugAssertOp(startColumn, < , this->columnCount());
        debugAssertOp(rowCount, <= , this->rowCount() - startRow);
        debugAssertOp(columnCount, <= , this->columnCount() - startColumn);
        if constexpr (detect<MatrixType, subMatrixConstView_t>{}) {
            return static_cast<const MatrixType*>(this)->_subMatrixView(startRow, startColumn, rowCount, columnCount);
        }  else if constexpr(CEffCol() != DynamicValue) {
            return MatrixConstStaticEffectiveColumnView<typename MatrixType::ValueType, CEffCol()>(*this, startRow, startColumn, rowCount, columnCount);
        } else {
            return MatrixConstDynamicView<typename MatrixType::ValueType>(*this, startRow, startColumn, rowCount, columnCount);
        }
    }

    template<mtl::size_t RowCount, mtl::size_t ColumnCount> auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn) {
        debugAssertOp(startRow, < , this->rowCount());
        debugAssertOp(startColumn, < , this->columnCount());
        debugAssertOp(RowCount, <= , this->rowCount() - startRow);
        debugAssertOp(ColumnCount, <= , this->columnCount() - startColumn);
        if constexpr (CEffCol() != DynamicValue) {
            return MatrixStaticView<typename MatrixType::ValueType, RowCount, ColumnCount, CEffCol()>(*this, startRow, startColumn);
        } else {
            return MatrixStaticSizeView<typename MatrixType::ValueType, RowCount, ColumnCount>(*this, startRow, startColumn);
        }
    }
    template<mtl::size_t RowCount, mtl::size_t ColumnCount> auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn) const {
        debugAssertOp(startRow, < , this->rowCount());
        debugAssertOp(startColumn, < , this->columnCount());
        debugAssertOp(RowCount, <= , this->rowCount() - startRow);
        debugAssertOp(ColumnCount, <= , this->columnCount() - startColumn);
        if constexpr (CEffCol() != DynamicValue) {
            return MatrixConstStaticView<typename MatrixType::ValueType, RowCount, ColumnCount, CEffCol()>(*this, startRow, startColumn);
        } else {
            return MatrixConstStaticSizeView<typename MatrixType::ValueType, RowCount, ColumnCount>(*this, startRow, startColumn);
        }
    }
    auto begin() {
        auto& it = static_cast<MatrixType&>(*this);
        if constexpr (detect<MatrixType, begin_t>{}) {
            return it._begin();
        } else if constexpr (CCol() != DynamicValue && CEffCol() != DynamicValue) {
            return MatrixStaticRowIterator<typename MatrixType::ValueType, CCol(), CEffCol()>(it.data());
        } else {
            return MatrixRowIterator<typename MatrixType::ValueType>(it.data(), it.columnCount(), it.effectiveColumnCount());
        }
    }
    auto end() {
        auto& it = static_cast<MatrixType&>(*this);
        if constexpr (detect<MatrixType, begin_t>{}) {
            return it._end();
        } else if constexpr (CCol() != DynamicValue && CEffCol() != DynamicValue) {
            return MatrixStaticRowIterator<typename MatrixType::ValueType, CCol(), CEffCol()>(it.data() + CEffCol() * it.rowCount());
        } else {
            return MatrixRowIterator<typename MatrixType::ValueType>(it.data() + it.effectiveColumnCount() * it.rowCount(), it.columnCount(), it.effectiveColumnCount());
        }
    }
    auto begin() const {
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, constbegin_t>{}) {
            return it._begin();
        }  else if constexpr (CCol() != DynamicValue && CEffCol() != DynamicValue) {
            return MatrixConstStaticRowIterator<typename MatrixType::ValueType, CCol(), CEffCol()>(it.data());
        } else {
            return MatrixConstRowIterator<typename MatrixType::ValueType>(it.data(), it.columnCount(), it.effectiveColumnCount());
        }
    }
    auto end() const {
        auto& it = static_cast<const MatrixType&>(*this);
        if constexpr (detect<MatrixType, constend_t>{}) {
            return it._end();
        } else if constexpr (CCol() != DynamicValue && CEffCol() != DynamicValue) {
            return MatrixConstStaticRowIterator<typename MatrixType::ValueType, CCol(), CEffCol()>(it.data() + CEffCol() * it.rowCount());
        } else {
            return MatrixConstRowIterator<typename MatrixType::ValueType>(it.data() + it.effectiveColumnCount() * it.rowCount(), it.columnCount(), it.effectiveColumnCount());
        }
    }
};
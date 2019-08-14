#ifndef MATRIX_EXTENDED_FUNCTIONS_H
#define MATRIX_EXTENDED_FUNCTIONS_H
#include "MatrixInterface.h"
#include <array>

template<typename M> auto shrink(const MatrixInterface<M>& matrix, int newRowCount, int newColumnCount) {
    auto result = matrix.createNew(newRowCount, newColumnCount);
    for (mtl::size_t row = 0; row < result.rowCount(); ++row) {
        for (mtl::size_t column = 0; column < result.columnCount(); ++column) {
            result.at(row, column) = matrix.at(row, column);
        }
    }
    return result;
}
template<int NewRowCount, int NewColumnCount, typename M> 
auto shrink(const MatrixInterface<M>& matrix) {
    auto result = matrix.template createNew<NewRowCount, NewColumnCount>();
    for (mtl::size_t row = 0; row < result.rowCount(); ++row) {
        for (mtl::size_t column = 0; column < result.columnCount(); ++column) {
            result.at(row, column) = matrix.at(row, column);
        }
    }
    return result;
}

template<mtl::size_t Rows, mtl::size_t Columns, typename M>
auto matrixDivideView(const MatrixInterface<M>& m) {
    std::array<std::array<decltype(std::declval<const M>().subMatrixView(0, 0, 0, 0)), Columns>, Rows> result;
    auto rowSize = m.rowCount() / Rows;
    auto columnSize = m.columnCount() / Columns;
    for (mtl::size_t y = 0; y < Rows; ++y) {
        for (mtl::size_t x = 0; x < Columns; ++x) {
            result[y][x] = m.subMatrixView(y*rowSize, x*columnSize, rowSize, columnSize);
        }
    }
    return result;
}
template<mtl::size_t Rows, mtl::size_t Columns, typename M>
auto matrixDivideView(MatrixInterface<M>& m) {
    if constexpr (M::HasConstexprRowAndColumnCount()) {
        constexpr auto RowSize = M::CRow() / Rows;
        constexpr auto ColumnSize = M::CCol() / Columns;
        std::array<std::array<decltype(std::declval<M>().template subMatrixView<RowSize, ColumnSize>(0, 0)), Columns>, Rows> result;
        for (mtl::size_t y = 0; y < Rows; ++y) {
            for (mtl::size_t x = 0; x < Columns; ++x) {
                result[y][x] = m.template subMatrixView<RowSize, ColumnSize>(y*RowSize, x*ColumnSize);
            }
        }
        return result;
    } else {
        std::array<std::array<decltype(std::declval<M>().subMatrixView(0, 0, 0, 0)), Columns>, Rows> result;
        auto rowSize = m.rowCount() / Rows;
        auto columnSize = m.columnCount() / Columns;
        for (mtl::size_t y = 0; y < Rows; ++y) {
            for (mtl::size_t x = 0; x < Columns; ++x) {
                result[y][x] = m.subMatrixView(y*rowSize, x*columnSize, rowSize, columnSize);
            }
        }
        return result;
    }
}
template<mtl::size_t Rows, mtl::size_t Columns, typename M>
auto matrixDivide(const MatrixInterface<M>& m) {
    if constexpr (M::HasConstexprRowAndColumnCount()) {
        constexpr auto RowSize = M::CRow() / Rows;
        constexpr auto ColumnSize = M::CCol() / Columns;
        std::array<std::array<decltype(std::declval<M>().template createNew<RowSize, ColumnSize>()), Columns>, Rows> result;
        for (mtl::size_t y = 0; y < Rows; ++y) {
            for (mtl::size_t x = 0; x < Columns; ++x) {
                result[y][x] = m.template subMatrixView<RowSize, ColumnSize>(y*RowSize, x*ColumnSize);
            }
        }
        return result;
    } else {
        std::array<std::array<decltype(std::declval<M>().createNew()), Columns>, Rows> result;
        auto rowSize = m.rowCount() / Rows;
        auto columnSize = m.columnCount() / Columns;
        for (mtl::size_t y = 0; y < Rows; ++y) {
            for (mtl::size_t x = 0; x < Columns; ++x) {
                result[y][x] = m.subMatrixView(y*rowSize, x*columnSize, rowSize, columnSize);
            }
        }
        return result;
    }
}

#endif

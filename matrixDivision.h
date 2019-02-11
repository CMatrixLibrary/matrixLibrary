#pragma once
#include <tuple>
#include <array>
#include <cassert>
#include "FullMatrixView.h"
#include "FullMatrixConstView.h"

/*
    Evenly divide matrix and save as std::array[Rows][Columns] of FullMatrixView
    It's required that the matrix size can be divided to given Rows and Columns
*/
template<int Rows, int Columns, typename T, template<typename> typename Matrix>
std::array<std::array<FullMatrixConstView<T>, Columns>, Rows> matrixDivide(const Matrix<T>& m) {
    assert(m.rowCount() % Rows == 0);
    assert(m.columnCount() % Columns == 0);

    std::array<std::array<FullMatrixConstView<T>, Columns>, Rows> result;
    auto rowSize = m.rowCount() / Rows;
    auto columnSize = m.columnCount() / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[y][x] = m.subMatrixView(y*rowSize, x*columnSize, rowSize, columnSize);
        }
    }

    return result;
}
template<int Rows, int Columns, typename T>
std::array<std::array<FullMatrixView<T>, Columns>, Rows> matrixDivide(FullMatrix<T>& m) {
    assert(m.rowCount() % Rows == 0);
    assert(m.columnCount() % Columns == 0);

    std::array<std::array<FullMatrixView<T>, Columns>, Rows> result;
    auto rowSize = m.rowCount() / Rows;
    auto columnSize = m.columnCount() / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[y][x] = m.subMatrixView(y*rowSize, x*columnSize, rowSize, columnSize);
        }
    }

    return result;
}
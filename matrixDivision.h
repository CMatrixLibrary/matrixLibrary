#pragma once
#include <tuple>
#include <array>
#include <cassert>
#include "FullMatrixView.h"

namespace details {
    template <typename T, std::size_t...Is>
    std::array<T, sizeof...(Is)> make_array(const T& value, std::index_sequence<Is...>) {
        return { {(static_cast<void>(Is), value)...} };
    }

    template <std::size_t N, typename T>
    std::array<T, N> make_array(const T& value) {
        return make_array(value, std::make_index_sequence<N>());
    }
}

/*
    Evenly divide matrix and save as std::array[Rows][Columns] of FullMatrixView
    It's required that the matrix size can be divided to given Rows and Columns
*/
template<int Rows, int Columns, typename T, template<typename> typename Matrix>
std::array<std::array<FullMatrixView<T>, Columns>, Rows> matrixDivide(const Matrix<T>& m) {
    assert(m.rowsCount() % Rows == 0);
    assert(m.columnsCount() % Columns == 0);

    auto result = details::make_array<Rows, std::array<FullMatrixView<T>, Columns>>(
        details::make_array<Columns>(FullMatrixView(m))
    );
    auto rowSize = m.rowsCount() / Rows;
    auto columnSize = m.columnsCount() / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[y][x] = FullMatrixView(m, y*rowSize, x*columnSize, rowSize, columnSize);
        }
    }
    return result;
}
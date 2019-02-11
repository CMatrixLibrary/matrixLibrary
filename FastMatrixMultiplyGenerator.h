#pragma once
#include <array>
#include <tuple>
#include "FullMatrix.h"
#include "FullMatrixView.h"
#include "naiveOperations.h"
#include "operators.h"
#include "matrixDivision.h"

template <int... Ints> struct Values {};
template <typename... Ts> struct Types {};

template<int Rows, int Columns, typename T, template<typename> typename Matrix>
std::array<FullMatrixConstView<T>, Columns * Rows> matrixDivideToConstView(const Matrix<T>& m) {
    assert(m.rowCount() % Rows == 0);
    assert(m.columnCount() % Columns == 0);

    std::array<FullMatrixConstView<T>, Columns * Rows> result;
    auto rowSize = m.rowCount() / Rows;
    auto columnSize = m.columnCount() / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[x + y * Columns] = FullMatrixConstView<T>(m, y*rowSize, x*columnSize, rowSize, columnSize);
        }
    }
    return result;
}
template<int Rows, int Columns, typename T, template<typename> typename Matrix>
std::array<FullMatrixView<T>, Columns * Rows> matrixDivideToView(Matrix<T>& m) {
    assert(m.rowCount() % Rows == 0);
    assert(m.columnCount() % Columns == 0);

    std::array<FullMatrixView<T>, Columns * Rows> result;
    auto rowSize = m.rowCount() / Rows;
    auto columnSize = m.columnCount() / Columns;
    for (int y = 0; y < Rows; ++y) {
        for (int x = 0; x < Columns; ++x) {
            result[x + y * Columns] = FullMatrixView<T>(m, y*rowSize, x*columnSize, rowSize, columnSize);
        }
    }
    return result;
}

template<int Index, int SizeX, int SizeY, typename T, int coff, int... CoffInts>
FullMatrix<T> apply(
    const std::array<FullMatrixConstView<T>, SizeX * SizeY>& dA,
    Values<coff, CoffInts...> coffs
) {
    if constexpr (Index == SizeX * SizeY - 1) {
        return coff * dA[Index];
    } else if constexpr (coff == 0) {
        return apply<Index + 1, SizeX, SizeY>(dA, Values<CoffInts...>());
    } else if constexpr (coff == 1) {
        return dA[Index] + apply<Index + 1, SizeX, SizeY>(dA, Values<CoffInts...>());
    } else {
        return coff * dA[Index] + apply<Index + 1, SizeX, SizeY>(dA, Values<CoffInts...>());
    }
}

template<
    int SizeX, int SizeY, int SizeZ,
    int MulCount,
    typename T,
    template<typename> typename MatrixA,
    template<typename> typename MatrixB,
    int ACoffInt,
    int BCoffInt,
    int CCoffInt,
    int... ACoffInts,
    int... BCoffInts,
    int... CCoffInts,
    typename... ACoffs,
    typename... BCoffs,
    typename... CCoffs
>
FullMatrix<T> multiplyRecursive(
    const MatrixA<T>& a,
    const MatrixB<T>& b,
    Types<Values<ACoffInt, ACoffInts...>, ACoffs...> aCoff,
    Types<Values<BCoffInt, BCoffInts...>, BCoffs...> bCoff,
    Types<Values<CCoffInt, CCoffInts...>, CCoffs...> cCoff
);

template<
    int Index,
    int SizeX, int SizeY, int SizeZ,
    int MulCount,
    typename T,
    int... ACoffInts,
    int... BCoffInts,
    typename... ACoffs,
    typename... BCoffs,
    typename... ACoffTypes,
    typename... BCoffTypes,
    typename... CCoffTypes
>
void createMArray(
    std::array<FullMatrix<T>, MulCount>& matrix,
    std::array<FullMatrixConstView<T>, SizeX * SizeY>& dA,
    std::array<FullMatrixConstView<T>, SizeY * SizeZ>& dB,
    Types<Values<ACoffInts...>, ACoffs...> aCoffIter,
    Types<Values<BCoffInts...>, BCoffs...> bCoffIter,
    Types<ACoffTypes...> aCoff,
    Types<BCoffTypes...> bCoff,
    Types<CCoffTypes...> cCoff
) {
    if constexpr (Index < MulCount) {
        matrix[Index] = multiplyRecursive<SizeX, SizeY, SizeZ, MulCount, T>(
            apply<0, SizeX, SizeY>(dA, Values<ACoffInts...>()),
            apply<0, SizeY, SizeZ>(dB, Values<BCoffInts...>()),
            aCoff, bCoff, cCoff
        );
        createMArray<Index + 1, SizeX, SizeY, SizeZ>(matrix, dA, dB, Types<ACoffs...>(), Types<BCoffs...>(), aCoff, bCoff, cCoff);
    }
}


template<int Index, int MulCount, typename T, int coff, int... CoffInts>
void applyC(
    FullMatrixView<T>& dC,
    const std::array<FullMatrix<T>, MulCount>& m,
    Values<coff, CoffInts...> coffs
) {
    if constexpr (Index < MulCount) {
        if constexpr (coff == 1) {
            dC += m[Index];
        } else if constexpr (coff != 0) {
            dC += coff * m[Index];
        }
        applyC<Index + 1>(dC, m, Values<CoffInts...>());
    }
}
template<
    int Index,
    int SizeX, int SizeY,
    int MulCount,
    typename T,
    int... CCoffInts,
    typename... Coffs
>
void createCArray(
    std::array<FullMatrixView<T>, SizeX * SizeY>& dC,
    const std::array<FullMatrix<T>, MulCount>& m,
    Types<Values<CCoffInts...>, Coffs...> cCoff
) {
    if constexpr (Index < SizeX * SizeY) {
        applyC<0>(dC[Index], m, Values<CCoffInts...>());
        createCArray<Index + 1, SizeX, SizeY>(dC, m, Types<Coffs...>());
    }
}


template<
    int SizeX, int SizeY, int SizeZ, 
    int MulCount,
    typename T,
    template<typename> typename MatrixA,
    template<typename> typename MatrixB,
    int ACoffInt,
    int BCoffInt,
    int CCoffInt,
    int... ACoffInts,
    int... BCoffInts,
    int... CCoffInts,
    typename... ACoffs,
    typename... BCoffs,
    typename... CCoffs
>
FullMatrix<T> multiplyRecursive(
    const MatrixA<T>& a,
    const MatrixB<T>& b,
    Types<Values<ACoffInt, ACoffInts...>, ACoffs...> aCoff,
    Types<Values<BCoffInt, BCoffInts...>, BCoffs...> bCoff,
    Types<Values<CCoffInt, CCoffInts...>, CCoffs...> cCoff
) {
    if (a.rowCount() <= 32) {
        return naiveMul(a, b);
    }
    auto dA = matrixDivideToConstView<SizeX, SizeY>(a);
    auto dB = matrixDivideToConstView<SizeY, SizeZ>(b);

    std::array<FullMatrix<T>, MulCount> m;
    createMArray<0, SizeX, SizeY, SizeZ, MulCount, T>(m, dA, dB, aCoff, bCoff, aCoff, bCoff, cCoff);

    FullMatrix<T> c(a.rowCount(), b.columnCount());
    auto dC = matrixDivideToView<SizeX, SizeZ>(c);
    createCArray<0, SizeX, SizeY>(dC, m, cCoff);

    return c;
}


template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> generatedStrassen(const MatrixA<T>& a, const MatrixB<T>& b) {
    Types <
        Values<1, 0, 0, 1>,
        Values<0, 0, 1, 1>,
        Values<1, 0, 0, 0>,
        Values<0, 0, 0, 1>,
        Values<1, 1, 0, 0>,
        Values<-1, 0, 1, 0>,
        Values<0, 1, 0, -1>,
        Values<>
    > aCoff;
    Types <
        Values<1, 0, 0, 1>,
        Values<1, 0, 0, 0>,
        Values<0, 1, 0, -1>,
        Values<-1, 0, 1, 0>,
        Values<0, 0, 0, 1>,
        Values<1, 1, 0, 0>,
        Values<0, 0, 1, 1>,
        Values<>
    > bCoff;
    Types <
        Values<1, 0, 0, 1, -1, 0, 1, 1234>,
        Values<0, 0, 1, 0, 1, 0, 0, 1234>,
        Values<0, 1, 0, 1, 0, 0, 0, 1234>,
        Values<1, -1, 1, 0, 0, 1, 0, 1234>,
        Values<>
    > cCoff;
    return multiplyRecursive<2, 2, 2, 7>(a, b, aCoff, bCoff, cCoff);
}

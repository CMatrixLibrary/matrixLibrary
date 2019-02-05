#pragma once
#include <iostream>
#include <iomanip>
#include "naiveOperations.h"
#include "strassenMultiply.h"

template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> operator+(const MatrixA<T>& a, const MatrixB<T>& b) {
    return naiveAdd(a, b);
}
template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> operator-(const MatrixA<T>& a, const MatrixB<T>& b) {
    return naiveSub(a, b);
}
template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
FullMatrix<T> operator*(const MatrixA<T>& a, const MatrixB<T>& b) {
    return strassenMul(a, b);
}

template<typename T, template<typename> typename Matrix> 
std::ostream& operator<<(std::ostream& out, Matrix<T>& m) {
    for (int i = 0; i < m.rowsCount(); ++i) {
        for (int j = 0; j < m.columnsCount(); ++j) {
            out << std::setw(5) << m.at(j, i);
        }
        out << '\n';
    }
    return out;
}
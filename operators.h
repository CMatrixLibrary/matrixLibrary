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
FullMatrix<T> operator-(const Matrix<T>& m) {
    return naiveNeg(m);
}
template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
void operator+=(MatrixA<T>& a, const MatrixB<T>& b) {
    return naiveAddAssign(a, b);
}
template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>
void operator-=(MatrixA<T>& a, const MatrixB<T>& b) {
    return naiveSubAssign(a, b);
}

template<typename T, template<typename> typename Matrix>
FullMatrix<T> operator*(const Matrix<T>& m, const T& scalar) {
    return naiveMul(m, scalar);
}
template<typename T, template<typename> typename Matrix>
FullMatrix<T> operator*(const T& scalar, const Matrix<T>& m) {
    return naiveMul(m, scalar);
}

template<typename T, template<typename> typename Matrix> 
std::ostream& operator<<(std::ostream& out, Matrix<T>& m) {
    for (int i = 0; i < m.rowCount(); ++i) {
        for (int j = 0; j < m.columnCount(); ++j) {
            out << std::setw(5) << m.at(i, j);
        }
        out << '\n';
    }
    return out;
}
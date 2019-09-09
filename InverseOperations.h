#ifndef INVERSE_OPERATIONS_H
#define INVERSE_OPERATIONS_H
#include "MatrixInterface.h"
#include "Matrix.h"

template<typename M>
void zeroLowerTriangularColumnWithGauss(MatrixInterface<M>& m, int columnToZero, MatrixInterface<M>& eye) {
	int n = m.rowCount();
	auto zeroingElement = m.at(columnToZero, columnToZero);
	for (mtl::size_t i = columnToZero + 1; i < n; ++i) {
		auto multiplier = m.at(i, columnToZero) / zeroingElement;
		for (mtl::size_t j = 0; j < n; ++j) {
			m.at(i, j) -= multiplier * m.at(columnToZero, j);
			eye.at(i, j) -= multiplier * eye.at(columnToZero, j);
		}
	}
}

template<typename M>
void zeroUpperTriangularColumnWithGauss(MatrixInterface<M>& m, int columnToZero, MatrixInterface<M>& eye) {
	int n = m.rowCount();
	auto zeroingElement = m.at(columnToZero, columnToZero);
	for (mtl::size_t i = columnToZero - 1; i >= 0; --i) {
		auto multiplier = m.at(i, columnToZero) / zeroingElement;
		for (mtl::size_t j = n - 1; j >= 0; --j) {
			m.at(i, j) -= multiplier * m.at(columnToZero, j);
			eye.at(i, j) -= multiplier * eye.at(columnToZero, j);
		}
	}
}

template<typename M>
auto& gaussInverse(MatrixInterface<M>& matrixToInverse, MatrixInterface<M>& eye) {
	int n = matrixToInverse.columnCount();
	for (mtl::size_t i = 0; i < n - 1; ++i) {
		zeroLowerTriangularColumnWithGauss(matrixToInverse, i, eye);
	}
	for (mtl::size_t i = n - 1; i > 0; --i) {
		zeroUpperTriangularColumnWithGauss(matrixToInverse, i, eye);
	}
	for (mtl::size_t i = 0; i < n; ++i) {
		for (mtl::size_t j = 0; j < n; ++j) {
			eye.at(i, j) /= matrixToInverse.at(i, i);
		}
	}

	std::cout << matrixToInverse << std::endl;
	std::cout << eye << std::endl;
	return matrixToInverse;
}

#endif

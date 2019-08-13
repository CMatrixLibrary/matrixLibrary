#pragma once
#include "MatrixInterface.h"
#include "Matrix.h"

template<typename M> auto transpose(const MatrixInterface<M>& matrix) {
    auto result = matrix.createNew();
    for (int i = 0; i < matrix.rowCount(); ++i) {
        for (int j = 0; j < matrix.columnCount(); ++j) {
            result.at(j, i) = matrix.at(i, j);
        }
    }
    return result;
}

template<typename MR, typename M1, typename M2>
void naiveMul(MatrixInterface<MR>& r, const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    for (mtl::size_t i = 0; i < a.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < b.columnCount(); ++j) {
            typename MR::ValueType sum{};
            for (mtl::size_t k = 0; k < a.columnCount(); ++k) {
                sum += a.at(i, k) * b.at(k, j);
            }
            r.at(i, j) = sum;
        }
    }
}

template<typename M1, typename M2>
auto naiveMul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
        naiveMul(result, a, b);
        return result;
    } else {
        Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
        naiveMul(result, a, b);
        return result;
    }
}

template<typename M1, typename M2>
auto naiveAdd(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    auto result = a.createNew();
    for (mtl::size_t i = 0; i < result.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < result.columnCount(); ++j) {
            result.at(i, j) = a.at(i, j) + b.at(i, j);
        }
    }
    return result;
}

template<typename M1, typename M2>
auto& naiveAddAssign(MatrixInterface<M1>& result, const MatrixInterface<M2>& m) {
    for (mtl::size_t i = 0; i < result.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < result.columnCount(); ++j) {
            result.at(i, j) += m.at(i, j);
        }
    }
    return result;
}

template<typename M1, typename M2>
auto naiveSub(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
    auto result = a.createNew();
    for (mtl::size_t i = 0; i < result.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < result.columnCount(); ++j) {
            result.at(i, j) = a.at(i, j) - b.at(i, j);
        }
    }
    return result;
}

template<typename M1, typename M2>
auto& naiveSubAssign(MatrixInterface<M1>& result, const MatrixInterface<M2>& m) {
    for (mtl::size_t i = 0; i < result.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < result.columnCount(); ++j) {
            result.at(i, j) -= m.at(i, j);
        }
    }
    return result;
}

template<typename M>
auto naiveNeg(const MatrixInterface<M>& m) {
    auto result = m.createNew();
    for (mtl::size_t i = 0; i < result.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < result.columnCount(); ++j) {
            result.at(i, j) -= m.at(i, j);
        }
    }
    return result;
}

template<typename M>
auto naiveScalarMul(const MatrixInterface<M>& m, const typename M::ValueType& scalar) {
    auto result = m.createNew();
    for (mtl::size_t i = 0; i < result.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < result.columnCount(); ++j) {
            result.at(i, j) = m.at(i, j) * scalar;
        }
    }
    return result;
}
template<typename M>
auto& naiveScalarMulAssign(MatrixInterface<M>& result, const typename M::ValueType& scalar) {
    for (mtl::size_t i = 0; i < result.rowCount(); ++i) {
        for (mtl::size_t j = 0; j < result.columnCount(); ++j) {
            result.at(i, j) *= scalar;
        }
    }
    return result;
}

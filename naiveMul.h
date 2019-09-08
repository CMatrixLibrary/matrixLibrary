#ifndef NAIVE_MUL_H
#define NAIVE_MUL_H

#include "MatrixInterface.h"
#include "baseMulUtility.h"

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

CREATE_BASE_MUL_WRAPPERS(naiveMul)

#endif

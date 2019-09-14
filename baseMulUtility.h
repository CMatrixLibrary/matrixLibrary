#ifndef BASE_MUL_UTILITY_H
#define BASE_MUL_UTILITY_H

#include "MatrixInterface.h"
#include "Matrix.h"
#include "MatrixView.h"

// macro creating high and low level wrappers which call base multiplication function
#define CREATE_BASE_MUL_WRAPPERS(FUNCTION_NAME) \
template<typename M1, typename M2>\
auto FUNCTION_NAME(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {\
    if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {\
        Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;\
        FUNCTION_NAME(result, a, b);\
        return result;\
    } else {\
        auto result = a.createNew(a.rowCount(), b.columnCount()); \
        FUNCTION_NAME(result, a, b);\
        return result;\
    }\
}\
template<typename T> void FUNCTION_NAME(T* c, const T* a, const T* b, int n, int m, int p, int effC, int effA, int effB) {\
    MatrixView<T> cc(c, n, p, effC);\
    FUNCTION_NAME(cc, MatrixConstView<T>(a, n, m, effA), MatrixConstView<T>(b, m, p, effB));\
}\
template<typename T> void FUNCTION_NAME(T* c, const T* a, const T* b, int n, int m, int p) {\
    FUNCTION_NAME(c, a, b, n, m, p, p, m, p);\
}\
template<int n, int m, int p, int effC, int effA, int effB, typename T> void FUNCTION_NAME(T* c, const  T* a, const T* b) {\
    MatrixView<T, n, p, effC> cc(c);\
    FUNCTION_NAME(cc, MatrixConstView<T, n, m, effA>(a), MatrixConstView<T, m, p, effB>(b));\
}\
template<int n, int m, int p, typename T> void FUNCTION_NAME(T* c, const T* a, const T* b) {\
    return FUNCTION_NAME<n, m, p, p, m, p>(c, a, b);\
}

#endif
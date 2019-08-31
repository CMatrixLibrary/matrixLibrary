#ifndef BLAS_MUL_H
#define BLAS_MUL_H

#include <complex>
#include <type_traits>
#include "always_false.h"
#include "Matrix.h"

#define Blas_StaticAssertMessage "blas is not available. Maybe you're missing \"#define USE_BLAS\"?"

namespace blas {
    template<typename T> constexpr bool IsCompatible =
        std::is_same_v<T, float> || std::is_same_v<T, double> ||
        std::is_same_v<T, std::complex<float>> || std::is_same_v<T, std::complex<double>>;

#ifdef USE_BLAS
    #define BLAS_IS_AVAILABLE
    constexpr bool IsAvailable = true;
#else
    constexpr bool IsAvailable = false;
#endif
}

namespace blas::detail {
    template<typename T> void mul(T* c, const T* a, const T* b, int n, int m, int q, int effC, int effA, int effB) {
        if constexpr (IsAvailable) {
            static_assert(always_false_v<T>, "blas function is not available for this type");
        } else {
            static_assert(always_false_v<T>, Blas_StaticAssertMessage);
        }
    }

    template<typename T> void setToZero(T* a, int n, int m, int effM) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                a[j + i * effM] = T{};
            }
        }
    }

#ifdef USE_BLAS
    template<> void mul<float>(float* c, const float* a, const float* b, int n, int m, int q, int effC, int effA, int effB) {
        setToZero(c, n, m, effC);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, 1, a, effA, b, effB, 1, c, effC);
    }
    template<> void mul<double>(double* c, const double* a, const double* b, int n, int m, int q, int effC, int effA, int effB) {
        setToZero(c, n, m, effC);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, 1, a, effA, b, effB, 1, c, effC);
    }
    template<> void mul<std::complex<float>>(std::complex<float>* c, const std::complex<float>* a, const std::complex<float>* b, int n, int m, int q, int effC, int effA, int effB) {
        setToZero(c, n, m, effC);
        auto one = std::complex<float>(1);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, &one, a, effA, b, effB, &one, c, effC);
    }
    template<> void mul<std::complex<double>>(std::complex<double>* c, const std::complex<double>* a, const std::complex<double>* b, int n, int m, int q, int effC, int effA, int effB) {
        setToZero(c, n, m, effC);
        auto one = std::complex<double>(1);
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, &one, a, effA, b, effB, &one, c, effC);
    }
#endif
}

namespace blas {

    template<typename T> void mul(T* c, const T* a, const T* b, int n, int m, int q, int effC, int effA, int effB) {
        detail::mul(c, a, b, n, m, q, effC, effA, effB);
    }
    template<typename T> void mul(T* c, const T* a, const T* b, int n, int m, int q) {
        detail::mul(c, a, b, n, m, q, q, m, q);
    }

    template<typename M1, typename M2>
    auto mul(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b) {
        if constexpr (IsAvailable) {
            if constexpr (M1::HasConstexprRowAndColumnCount() && M2::HasConstexprRowAndColumnCount()) {
                Matrix<typename M1::ValueType, M1::CRow(), M2::CCol()> result;
                for (int i = 0; i < a.size(); ++i) result.data()[i] = typename M1::ValueType{};
                mul(result.data(), a.data(), b.data(), M1::CRow(), M1::CCol(), M2::CCol(), result.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
                return result;
            } else {
                Matrix<typename M1::ValueType> result(a.rowCount(), b.columnCount());
                for (int i = 0; i < a.size(); ++i) result.data()[i] = typename M1::ValueType{};
                mul(result.data(), a.data(), b.data(), a.rowCount(), a.columnCount(), b.columnCount(), result.effectiveColumnCount(), a.effectiveColumnCount(), b.effectiveColumnCount());
                return result;
            }
        } else {
            static_assert(always_false_v<M1>, "blas is not available");
        }
    }
}

#endif

#pragma once

#include <complex>
#include <type_traits>
#include "always_false.h"

namespace blas {
#ifdef USE_BLAS
    #define BLAS_IS_AVAILABLE
    constexpr bool IsAvailable = true;
#else
    constexpr bool IsAvailable = false;
#endif

    template<typename T> void mul(T* c, const T* a, const T* b, int n, int m, int q, int effC, int effA, int effB) {
        if constexpr (IsAvailable) {
            static_assert(always_false_v<T>, "blas function is not available for this type");
        } else {
            static_assert(always_false_v<T>, "blas functions are not available, missing \"#define USE_BLAS\"?");
        }
    }

#ifdef USE_BLAS
    template<> void mul<float>(float* c, const float* a, const float* b, int n, int m, int q, int effC, int effA, int effB) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, 1, a, effA, b, effB, 1, c, effC);
    }
    template<> void mul<double>(double* c, const double* a, const double* b, int n, int m, int q, int effC, int effA, int effB) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, 1, a, effA, b, effB, 1, c, effC);
    }
    template<> void mul<std::complex<float>>(std::complex<float>* c, const std::complex<float>* a, const std::complex<float>* b, int n, int m, int q, int effC, int effA, int effB) {
        auto one = std::complex<float>(1);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, &one, a, effA, b, effB, &one, c, effC);
    }
    template<> void mul<std::complex<double>>(std::complex<double>* c, const std::complex<double>* a, const std::complex<double>* b, int n, int m, int q, int effC, int effA, int effB) {
        auto one = std::complex<double>(1);
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, m, &one, a, effA, b, effB, &one, c, effC);
    }
#endif

    template<typename T> void mul(T* c, const T* a, const T* b, int n, int m, int q) {
        mul(c, a, b, n, m, q, q, m, q);
    }
}

#pragma once
#include <immintrin.h>
#include <type_traits>
#include <new>
#include <cstdint>
#include "alignedAllocation.h"

template<typename ValueType> class AVX256Type;

template<> class AVX256Type<int32_t> {
    __m256i value;
public:
    AVX256Type(__m256i value) : value(value) {}
    operator __m256i() { return value; }
};
template<> class AVX256Type<float> {
    __m256 value;
public:
    AVX256Type(__m256 value) : value(value) {}
    operator __m256() { return value; }
};
template<> class AVX256Type<double> {
    __m256d value;
public:
    AVX256Type(__m256d value) : value(value) {}
    operator __m256d() { return value; }
};

namespace AVX256 {
    template<typename T> constexpr int packedCount() {
        return 256/8 / sizeof(T);
    }

    std::size_t Alignment = 32;

    template<typename T> T* alignedAlloc() {
        return alignedMalloc<T>(Alignment);
    }
    template<typename T> T* alignedArrayAlloc(std::size_t size) {
        return alignedMalloc<T>(size, Alignment);
    }
    void alignedDealloc(void* ptr) {
        alignedFree(ptr);
    }
    void alignedArrayDealloc(void* ptr) {
        alignedFree(ptr);
    }

    AVX256Type<int32_t> loadAligned(const int32_t* ptr) { return _mm256_load_si256((__m256i*)ptr); }
    AVX256Type<float> loadAligned(const float* ptr) { return _mm256_load_ps(ptr); }
    AVX256Type<double> loadAligned(const double* ptr) { return _mm256_load_pd(ptr); }

    AVX256Type<int32_t> loadUnaligned(const int32_t* ptr) { return _mm256_loadu_si256((__m256i*)ptr); }
    AVX256Type<float> loadUnaligned(const float* ptr) { return _mm256_loadu_ps(ptr); }
    AVX256Type<double> loadUnaligned(const double* ptr) { return _mm256_loadu_pd(ptr); }

    void storeAligned(int32_t* dst, AVX256Type<int32_t> src) { return _mm256_store_si256((__m256i*)dst, src); }
    void storeAligned(float* dst, AVX256Type<float> src) { return _mm256_store_ps(dst, src); }
    void storeAligned(double* dst, AVX256Type<double> src) { return _mm256_store_pd(dst, src); }

    void storeUnaligned(int32_t* dst, AVX256Type<int32_t> src) { return _mm256_storeu_si256((__m256i*)dst, src); }
    void storeUnaligned(float* dst, AVX256Type<float> src) { return _mm256_storeu_ps(dst, src); }
    void storeUnaligned(double* dst, AVX256Type<double> src) { return _mm256_storeu_pd(dst, src); }

    AVX256Type<int32_t> add(AVX256Type<int32_t> a, AVX256Type<int32_t> b) { return _mm256_add_epi32(a, b); }
    AVX256Type<float> add(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_add_ps(a, b); }
    AVX256Type<double> add(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_add_pd(a, b); }

    AVX256Type<int32_t> sub(AVX256Type<int32_t> a, AVX256Type<int32_t> b) { return _mm256_sub_epi32(a, b); }
    AVX256Type<float> sub(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_sub_ps(a, b); }
    AVX256Type<double> sub(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_sub_pd(a, b); }

    AVX256Type<int32_t> mul(AVX256Type<int32_t> a, AVX256Type<int32_t> b) { return _mm256_mul_epi32(a, b); }
    AVX256Type<float> mul(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_mul_ps(a, b); }
    AVX256Type<double> mul(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_mul_pd(a, b); }

    template<typename T> AVX256Type<T> zero() { return AVX256Type<T>(0); }
    template<> AVX256Type<int32_t> zero() { return _mm256_setzero_si256(); }
    template<> AVX256Type<float> zero() { return _mm256_setzero_ps(); }
    template<> AVX256Type<double> zero() { return _mm256_setzero_pd(); }

    AVX256Type<int32_t> setAllElements(int32_t value) { return _mm256_set1_epi32(value); }
    AVX256Type<float> setAllElements(float value) { return _mm256_set1_ps(value); }
    AVX256Type<double> setAllElements(double value) { return _mm256_set1_pd(value); }

    AVX256Type<int32_t> fma(AVX256Type<int32_t> a, AVX256Type<int32_t> b, AVX256Type<int32_t> c) {
        return add(mul(a, b), c);
    }
    AVX256Type<float> fma(AVX256Type<float> a, AVX256Type<float> b, AVX256Type<float> c) {
        return _mm256_fmadd_ps(a, b, c);
    }
    AVX256Type<double> fma(AVX256Type<double> a, AVX256Type<double> b, AVX256Type<double> c) {
        return _mm256_fmadd_pd(a, b, c);
    }
};

template<typename T> AVX256Type<T> operator+(AVX256Type<T> a, AVX256Type<T> b) {
    return AVX256::add(a, b);
}
template<typename T> AVX256Type<T> operator-(AVX256Type<T> a, AVX256Type<T> b) {
    return AVX256::sub(a, b);
}
template<typename T> AVX256Type<T> operator*(AVX256Type<T> a, AVX256Type<T> b) {
    return AVX256::mul(a, b);
}
template<typename T> AVX256Type<T> operator+=(AVX256Type<T>& a, AVX256Type<T> b) {
    return a = a + b;
}
template<typename T> AVX256Type<T> operator-=(AVX256Type<T>& a, AVX256Type<T> b) {
    return a = a - b;
}
template<typename T> AVX256Type<T> operator*=(AVX256Type<T>& a, AVX256Type<T> b) {
    return a = a * b;
}
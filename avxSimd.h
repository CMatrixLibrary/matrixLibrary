#pragma once
#include <immintrin.h>
#include <type_traits>
#include <new>

template<typename ValueType> class AVX256Type;

template<> class AVX256Type<int> {
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
        return sizeof(T) / 256;
    }

    auto Alignment = static_cast<std::align_val_t>(32);

    template<typename T> T* AlignedAlloc() {
        return reinterpret_cast<T*>(operator new(sizeof(T), Alignment));
    }
    template<typename T> T* AlignedArrayAlloc(std::size_t size) {
        return reinterpret_cast<T*>(operator new[](sizeof(T)*size, Alignment));
    }
    void AlignedDealloc(void* ptr) {
        operator delete(ptr, Alignment);
    }
    void AlignedArrayDealloc(void* ptr) {
        operator delete[](ptr, Alignment);
    }
    

    AVX256Type<int> loadAligned(const int* ptr) { return _mm256_load_si256((__m256i*)ptr); }
    AVX256Type<float> loadAligned(const float* ptr) { return _mm256_load_ps(ptr); }
    AVX256Type<double> loadAligned(const double* ptr) { return _mm256_load_pd(ptr); }

    AVX256Type<int> loadUnaligned(const int* ptr) { return _mm256_loadu_si256((__m256i*)ptr); }
    AVX256Type<float> loadUnaligned(const float* ptr) { return _mm256_loadu_ps(ptr); }
    AVX256Type<double> loadUnaligned(const double* ptr) { return _mm256_loadu_pd(ptr); }

    void storeAligned(int* dst, AVX256Type<int> src) { return _mm256_store_si256((__m256i*)dst, src); }
    void storeAligned(float* dst, AVX256Type<float> src) { return _mm256_store_ps(dst, src); }
    void storeAligned(double* dst, AVX256Type<double> src) { return _mm256_store_pd(dst, src); }

    void storeUnaligned(int* dst, AVX256Type<int> src) { return _mm256_storeu_si256((__m256i*)dst, src); }
    void storeUnaligned(float* dst, AVX256Type<float> src) { return _mm256_storeu_ps(dst, src); }
    void storeUnaligned(double* dst, AVX256Type<double> src) { return _mm256_storeu_pd(dst, src); }

    AVX256Type<int> add(AVX256Type<int> a, AVX256Type<int> b) { return _mm256_add_epi32(a, b); }
    AVX256Type<float> add(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_add_ps(a, b); }
    AVX256Type<double> add(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_add_pd(a, b); }

    AVX256Type<int> sub(AVX256Type<int> a, AVX256Type<int> b) { return _mm256_sub_epi32(a, b); }
    AVX256Type<float> sub(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_sub_ps(a, b); }
    AVX256Type<double> sub(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_sub_pd(a, b); }

    AVX256Type<int> mul(AVX256Type<int> a, AVX256Type<int> b) { return _mm256_mullo_epi32(a, b); }
    AVX256Type<float> mul(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_mul_ps(a, b); }
    AVX256Type<double> mul(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_mul_pd(a, b); }

    template<typename T> AVX256Type<T> zero() { return AVX256Type<T>(0); }
    template<> AVX256Type<int> zero() { return _mm256_setzero_si256(); }
    template<> AVX256Type<float> zero() { return _mm256_setzero_ps(); }
    template<> AVX256Type<double> zero() { return _mm256_setzero_pd(); }

    AVX256Type<int> setAllElements(int value) { return _mm256_set1_epi32(value); }
    AVX256Type<float> setAllElements(float value) { return _mm256_set1_ps(value); }
    AVX256Type<double> setAllElements(double value) { return _mm256_set1_pd(value); }
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
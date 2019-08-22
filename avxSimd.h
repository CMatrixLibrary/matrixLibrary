#ifndef AVX_SIMD_H
#define AVX_SIMD_H
#include <immintrin.h>
#include <type_traits>
#include <new>
#include <cstdint>
#include <limits>
#include "alignedAllocation.h"
#include "always_false.h"
#include "compilerMacros.h"

#define AVX256_StaticAssertMessage "avx2 is not available. Maybe you're missing a compilation flag?"

namespace AVX256 {
#ifdef __AVX2__
    #define AVX2_IS_AVAILABLE
    #if (defined(COMPILER_MSVC) || defined(COMPILER_INTEL)) && !defined(__FMA__)
        #define __FMA__
    #endif

    constexpr bool IsAvailable = true;        
#else
    constexpr bool IsAvailable = false;
#endif

#ifdef __FMA__
    #define FMA_IS_AVAILABLE
#endif
}

template<typename ValueType, typename Enable = void> class AVX256Type;

#ifdef AVX2_IS_AVAILABLE
template<typename ValueType> class AVX256Type<ValueType, typename std::enable_if_t<std::is_integral_v<ValueType>>> {
    __m256i value;
public:
    AVX256Type(__m256i value) : value(value) {}
    operator __m256i() { return value; }
};
template<> class AVX256Type<float, typename std::enable_if_t<sizeof(float) == 4 && std::numeric_limits<float>::is_iec559>> {
    __m256 value;
public:
    AVX256Type(__m256 value) : value(value) {}
    operator __m256() { return value; }
};
template<> class AVX256Type<double, typename std::enable_if_t<sizeof(double) == 8 && std::numeric_limits<double>::is_iec559>> {
    __m256d value;
public:
    AVX256Type(__m256d value) : value(value) {}
    operator __m256d() { return value; }
};
#else
template<> class AVX256Type<float> {};
template<> class AVX256Type<double> {};
#endif

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

#ifdef AVX2_IS_AVAILABLE

    template<typename T> AVX256Type<T> loadAligned(const T* ptr) { 
        return _mm256_load_si256((__m256i*)ptr); 
    }
    AVX256Type<float> loadAligned(const float* ptr) { return _mm256_load_ps(ptr); }
    AVX256Type<double> loadAligned(const double* ptr) { return _mm256_load_pd(ptr); }

    template<typename T> AVX256Type<T> loadUnaligned(const T* ptr) { return _mm256_loadu_si256((__m256i*)ptr); }
    AVX256Type<float> loadUnaligned(const float* ptr) { return _mm256_loadu_ps(ptr); }
    AVX256Type<double> loadUnaligned(const double* ptr) { return _mm256_loadu_pd(ptr); }

    template<typename T> void storeAligned(T* dst, AVX256Type<T> src) { return _mm256_store_si256((__m256i*)dst, src); }
    void storeAligned(float* dst, AVX256Type<float> src) { return _mm256_store_ps(dst, src); }
    void storeAligned(double* dst, AVX256Type<double> src) { return _mm256_store_pd(dst, src); }

    template<typename T> void storeUnaligned(T* dst, AVX256Type<T> src) { return _mm256_storeu_si256((__m256i*)dst, src); }
    void storeUnaligned(float* dst, AVX256Type<float> src) { return _mm256_storeu_ps(dst, src); }
    void storeUnaligned(double* dst, AVX256Type<double> src) { return _mm256_storeu_pd(dst, src); }

    template<typename T> AVX256Type<T> add(AVX256Type<T> a, AVX256Type<T> b) { 
        if constexpr (sizeof(T) == 1) return _mm256_add_epi8(a, b);
        if constexpr (sizeof(T) == 2) return _mm256_add_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_add_epi32(a, b);
        if constexpr (sizeof(T) == 8) return _mm256_add_epi64(a, b);
    }
    AVX256Type<float> add(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_add_ps(a, b); }
    AVX256Type<double> add(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_add_pd(a, b); }

    template<typename T> AVX256Type<T> sub(AVX256Type<T> a, AVX256Type<T> b) {
        if constexpr (sizeof(T) == 1) return _mm256_sub_epi8(a, b);
        if constexpr (sizeof(T) == 2) return _mm256_sub_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_sub_epi32(a, b);
        if constexpr (sizeof(T) == 8) return _mm256_sub_epi64(a, b);
    }
    AVX256Type<float> sub(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_sub_ps(a, b); }
    AVX256Type<double> sub(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_sub_pd(a, b); }

    template<typename T> AVX256Type<T> mul(AVX256Type<T> a, AVX256Type<T> b) {
        if constexpr (sizeof(T) == 1) {
            auto even = _mm256_mullo_epi16(a, b);
            auto odd = _mm256_mullo_epi16(_mm256_srli_epi16(a, 8), _mm256_srli_epi16(b, 8));
            return _mm256_or_si256(_mm256_slli_epi16(odd, 8), _mm256_and_si256(even, _mm256_set1_epi16(0xFF)));
        }
        if constexpr (sizeof(T) == 2) return _mm256_mullo_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_mullo_epi32(a, b);
        if constexpr (sizeof(T) == 8) {
            // taken from https://stackoverflow.com/a/37320416
            __m256i bswap = _mm256_shuffle_epi32(b, 0xB1);
            __m256i prodlh = _mm256_mullo_epi32(a, bswap);

            __m256i prodlh2 = _mm256_srli_epi64(prodlh, 32);
            __m256i prodlh3 = _mm256_add_epi32(prodlh2, prodlh);
            __m256i prodlh4 = _mm256_and_si256(prodlh3, _mm256_set1_epi64x(0x00000000FFFFFFFF));

            __m256i prodll = _mm256_mul_epu32(a, b);
            __m256i prod = _mm256_add_epi64(prodll, prodlh4);
            return  prod;
        }
    }
    AVX256Type<float> mul(AVX256Type<float> a, AVX256Type<float> b) { return _mm256_mul_ps(a, b); }
    AVX256Type<double> mul(AVX256Type<double> a, AVX256Type<double> b) { return _mm256_mul_pd(a, b); }

    template<typename T> AVX256Type<T> zero() { return _mm256_setzero_si256(); }
    template<> AVX256Type<int8_t> zero() { return _mm256_setzero_si256(); }
    template<> AVX256Type<float> zero() { return _mm256_setzero_ps(); }
    template<> AVX256Type<double> zero() { return _mm256_setzero_pd(); }

    template<typename T> AVX256Type<T> setAllElements(T value) {
        if constexpr (sizeof(T) == 1) return _mm256_set1_epi8(value);
        if constexpr (sizeof(T) == 2) return _mm256_set1_epi16(value);
        if constexpr (sizeof(T) == 4) return _mm256_set1_epi32(value);
        if constexpr (sizeof(T) == 8) return _mm256_set1_epi64x(value);
    }
    AVX256Type<float> setAllElements(float value) { return _mm256_set1_ps(value); }
    AVX256Type<double> setAllElements(double value) { return _mm256_set1_pd(value); }

    template<typename T> AVX256Type<T> fma(AVX256Type<T> a, AVX256Type<T> b, AVX256Type<T> c) {
        return add(mul(a, b), c);
    }

    AVX256Type<float> fma(AVX256Type<float> a, AVX256Type<float> b, AVX256Type<float> c) {
        #ifdef FMA_IS_AVAILABLE
        return _mm256_fmadd_ps(a, b, c);
        #else
        return add(mul(a, b), c);
        #endif
    }
    AVX256Type<double> fma(AVX256Type<double> a, AVX256Type<double> b, AVX256Type<double> c) {
        #ifdef FMA_IS_AVAILABLE
        return _mm256_fmadd_pd(a, b, c);
        #else
        return add(mul(a, b), c);
        #endif
    }
#else
    template<typename T> AVX256Type<T> loadAligned(const T* ptr) { return AVX256Type<T>{}; }
    template<typename T> AVX256Type<T> loadUnaligned(const T* ptr) { return AVX256Type<T>{}; }
    template<typename T> void storeAligned(T* dst, AVX256Type<T> src) {}
    template<typename T> void storeUnaligned(T* dst, AVX256Type<T> src) {}
    template<typename T> AVX256Type<T> add(AVX256Type<T> a, AVX256Type<T> b) { return AVX256Type<T>{}; }
    template<typename T> AVX256Type<T> sub(AVX256Type<T> a, AVX256Type<T> b) { return AVX256Type<T>{}; }
    template<typename T> AVX256Type<T> mul(AVX256Type<T> a, AVX256Type<T> b) { return AVX256Type<T>{}; }
    template<typename T> AVX256Type<T> zero() { return AVX256Type<T>{}; }
    template<typename T> AVX256Type<T> setAllElements(T value) { return AVX256Type<T>{}; }
    template<typename T> AVX256Type<T> fma(AVX256Type<T> a, AVX256Type<T> b, AVX256Type<T> c) { return AVX256Type<T>{}; }
#endif
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

#endif

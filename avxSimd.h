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

#define avx_StaticAssertMessage "avx2 is not available. Maybe you're missing a compilation flag?"

namespace avx {
    template<typename T> constexpr bool IsCompatible = std::is_integral_v<T> || std::numeric_limits<T>::is_iec559;

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

template<typename ValueType, typename Enable = void> class AvxType;

#ifdef AVX2_IS_AVAILABLE
template<typename ValueType> class AvxType<ValueType, typename std::enable_if_t<std::is_integral_v<ValueType>>> {
    __m256i value;
public:
    AvxType(__m256i value) : value(value) {}
    operator __m256i() { return value; }
};
template<> class AvxType<float, typename std::enable_if_t<sizeof(float) == 4 && std::numeric_limits<float>::is_iec559>> {
    __m256 value;
public:
    AvxType(__m256 value) : value(value) {}
    operator __m256() { return value; }
};
template<> class AvxType<double, typename std::enable_if_t<sizeof(double) == 8 && std::numeric_limits<double>::is_iec559>> {
    __m256d value;
public:
    AvxType(__m256d value) : value(value) {}
    operator __m256d() { return value; }
};
#else
template<> class AvxType<float> {};
template<> class AvxType<double> {};
#endif

namespace avx {
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

    template<typename T> AvxType<T> loadAligned(const T* ptr) { 
        return _mm256_load_si256((__m256i*)ptr); 
    }
    AvxType<float> loadAligned(const float* ptr) { return _mm256_load_ps(ptr); }
    AvxType<double> loadAligned(const double* ptr) { return _mm256_load_pd(ptr); }

    template<typename T> AvxType<T> loadUnaligned(const T* ptr) { return _mm256_loadu_si256((__m256i*)ptr); }
    AvxType<float> loadUnaligned(const float* ptr) { return _mm256_loadu_ps(ptr); }
    AvxType<double> loadUnaligned(const double* ptr) { return _mm256_loadu_pd(ptr); }

    template<typename T> void storeAligned(T* dst, AvxType<T> src) { return _mm256_store_si256((__m256i*)dst, src); }
    void storeAligned(float* dst, AvxType<float> src) { return _mm256_store_ps(dst, src); }
    void storeAligned(double* dst, AvxType<double> src) { return _mm256_store_pd(dst, src); }

    template<typename T> void storeUnaligned(T* dst, AvxType<T> src) { return _mm256_storeu_si256((__m256i*)dst, src); }
    void storeUnaligned(float* dst, AvxType<float> src) { return _mm256_storeu_ps(dst, src); }
    void storeUnaligned(double* dst, AvxType<double> src) { return _mm256_storeu_pd(dst, src); }

    template<typename T> AvxType<T> add(AvxType<T> a, AvxType<T> b) { 
        if constexpr (sizeof(T) == 1) return _mm256_add_epi8(a, b);
        if constexpr (sizeof(T) == 2) return _mm256_add_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_add_epi32(a, b);
        if constexpr (sizeof(T) == 8) return _mm256_add_epi64(a, b);
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    AvxType<float> add(AvxType<float> a, AvxType<float> b) { return _mm256_add_ps(a, b); }
    AvxType<double> add(AvxType<double> a, AvxType<double> b) { return _mm256_add_pd(a, b); }

    template<typename T> AvxType<T> sub(AvxType<T> a, AvxType<T> b) {
        if constexpr (sizeof(T) == 1) return _mm256_sub_epi8(a, b);
        if constexpr (sizeof(T) == 2) return _mm256_sub_epi16(a, b);
        if constexpr (sizeof(T) == 4) return _mm256_sub_epi32(a, b);
        if constexpr (sizeof(T) == 8) return _mm256_sub_epi64(a, b);
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    AvxType<float> sub(AvxType<float> a, AvxType<float> b) { return _mm256_sub_ps(a, b); }
    AvxType<double> sub(AvxType<double> a, AvxType<double> b) { return _mm256_sub_pd(a, b); }

    template<typename T> AvxType<T> mul(AvxType<T> a, AvxType<T> b) {
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
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    AvxType<float> mul(AvxType<float> a, AvxType<float> b) { return _mm256_mul_ps(a, b); }
    AvxType<double> mul(AvxType<double> a, AvxType<double> b) { return _mm256_mul_pd(a, b); }

    template<typename T> AvxType<T> zero() { return _mm256_setzero_si256(); }
    template<> AvxType<int8_t> zero() { return _mm256_setzero_si256(); }
    template<> AvxType<float> zero() { return _mm256_setzero_ps(); }
    template<> AvxType<double> zero() { return _mm256_setzero_pd(); }

    template<typename T> AvxType<T> setAllElements(T value) {
        if constexpr (sizeof(T) == 1) return _mm256_set1_epi8(value);
        if constexpr (sizeof(T) == 2) return _mm256_set1_epi16(value);
        if constexpr (sizeof(T) == 4) return _mm256_set1_epi32(value);
        if constexpr (sizeof(T) == 8) return _mm256_set1_epi64x(value);
        return __m256i{}; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    }
    AvxType<float> setAllElements(float value) { return _mm256_set1_ps(value); }
    AvxType<double> setAllElements(double value) { return _mm256_set1_pd(value); }

    template<typename T> AvxType<T> fma(AvxType<T> a, AvxType<T> b, AvxType<T> c) {
        return add(mul(a, b), c);
    }

    AvxType<float> fma(AvxType<float> a, AvxType<float> b, AvxType<float> c) {
        #ifdef FMA_IS_AVAILABLE
        return _mm256_fmadd_ps(a, b, c);
        #else
        return add(mul(a, b), c);
        #endif
    }
    AvxType<double> fma(AvxType<double> a, AvxType<double> b, AvxType<double> c) {
        #ifdef FMA_IS_AVAILABLE
        return _mm256_fmadd_pd(a, b, c);
        #else
        return add(mul(a, b), c);
        #endif
    }
#else
    template<typename T> AvxType<T> loadAligned(const T* ptr) { return AvxType<T>{}; }
    template<typename T> AvxType<T> loadUnaligned(const T* ptr) { return AvxType<T>{}; }
    template<typename T> void storeAligned(T* dst, AvxType<T> src) {}
    template<typename T> void storeUnaligned(T* dst, AvxType<T> src) {}
    template<typename T> AvxType<T> add(AvxType<T> a, AvxType<T> b) { return AvxType<T>{}; }
    template<typename T> AvxType<T> sub(AvxType<T> a, AvxType<T> b) { return AvxType<T>{}; }
    template<typename T> AvxType<T> mul(AvxType<T> a, AvxType<T> b) { return AvxType<T>{}; }
    template<typename T> AvxType<T> zero() { return AvxType<T>{}; }
    template<typename T> AvxType<T> setAllElements(T value) { return AvxType<T>{}; }
    template<typename T> AvxType<T> fma(AvxType<T> a, AvxType<T> b, AvxType<T> c) { return AvxType<T>{}; }
#endif
};

template<typename T> AvxType<T> operator+(AvxType<T> a, AvxType<T> b) {
    return avx::add(a, b);
}
template<typename T> AvxType<T> operator-(AvxType<T> a, AvxType<T> b) {
    return avx::sub(a, b);
}
template<typename T> AvxType<T> operator*(AvxType<T> a, AvxType<T> b) {
    return avx::mul(a, b);
}
template<typename T> AvxType<T> operator+=(AvxType<T>& a, AvxType<T> b) {
    return a = a + b;
}
template<typename T> AvxType<T> operator-=(AvxType<T>& a, AvxType<T> b) {
    return a = a - b;
}
template<typename T> AvxType<T> operator*=(AvxType<T>& a, AvxType<T> b) {
    return a = a * b;
}

#endif

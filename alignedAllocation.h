#ifndef ALIGNED_ALLOCATION_H
#define ALIGNED_ALLOCATION_H

// Aligned versions of new, new[], delete, delete[] operators.
// new and new[] have additional "Construct" versions. Unlike standard versions they force initialization.
// This is because operators new and new[] do not initialize non-class objects, unless you add empty parentheses,
// which we simulate here by adding "Construct" to function name instead of parentheses.
//
// new vs alignedNew example comparison for non-class type (where 'Align' is chosen alignment):
// new int       == alignedNew<int>(Allign)
// new int()     == alignedNew<int>(Allign, 0)  {OR} alignedConstructNew<int>(Allign)
// new int(123)  == alignedNew<int>(Align, 123) {OR} alignedConstructNew<int>(Allign, 123)
// new int[10]   == alignedArrayNew<int>(10, Align)
// new int[10]() == alignedArrayNew<int>(10, Align, 0) {OR} alignedArrayConstructNew<int>(10, Align)
//
// In case of allocating class type the "Construct" and non-"Construct" variants are equivalent

#include <new>
#include <cstdlib>
#include <utility>
#include <type_traits>
#include <algorithm>
#include "compilerMacros.h"

// General aligned malloc/free working with any C++17 compliant compilator + Intel, Visual, g++ and clang compilers.
void* alignedMalloc(std::size_t size, std::size_t alignment) {
#if defined(COMPILER_INTEL) || defined(COMPILER_MSVC) || defined(COMPILER_GCC) || defined(COMPILER_CLANG)
    return _aligned_malloc(size, alignment);
#else
    return std::aligned_alloc(alignment, size);
#endif
}
void alignedFree(void* ptr) {
#if defined(COMPILER_INTEL) || defined(COMPILER_MSVC) || defined(COMPILER_GCC) || defined(COMPILER_CLANG)
    _aligned_free(ptr);
#else
    std::free(memory);
#endif
}


// Typed malloc
template<typename T> T* alignedMalloc(std::size_t alignment) {
    return reinterpret_cast<T*>(alignedMalloc(sizeof(T), alignment));
}
template<typename T> T* alignedMalloc(std::size_t size, std::size_t alignment) {
    return reinterpret_cast<T*>(alignedMalloc(sizeof(T)*size, alignment));
}


// array allocation with its size and starting address saved in the same block of memory
namespace detail {
    template<typename T> T* alignedSavedSizeMalloc(std::size_t size, std::size_t alignment) {
        auto spaceForSizeAndPtr = std::max(sizeof(std::size_t) + sizeof(void*), alignment);
        auto memory = alignedMalloc(sizeof(T)*size + spaceForSizeAndPtr, alignment);
        *reinterpret_cast<std::size_t*>(reinterpret_cast<char*>(memory) + spaceForSizeAndPtr - sizeof(std::size_t)) = size;
        *reinterpret_cast<void**>(reinterpret_cast<char*>(memory) + spaceForSizeAndPtr - sizeof(std::size_t) - sizeof(void*)) = memory;
        return reinterpret_cast<T*>(reinterpret_cast<char*>(memory) + spaceForSizeAndPtr);
    }
    void alignedSavedSizeFree(void* ptr) {
        alignedFree(*reinterpret_cast<void**>(reinterpret_cast<char*>(ptr) - sizeof(std::size_t) - sizeof(void*)));
    }
    std::size_t alignedSavedSizeGetSize(void* ptr) {
        return *reinterpret_cast<std::size_t*>(reinterpret_cast<char*>(ptr) - sizeof(std::size_t));
    }
}

// allocate + construct
template<typename T, typename... Args> T* alignedConstructNew(std::size_t alignment, Args&&... args) {
    return new(alignedMalloc<T>(alignment)) T(std::forward<Args>(args)...);
}
template<typename T, typename... Args> T* alignedArrayConstructNew(std::size_t size, std::size_t alignment, Args&&... args) {
    auto memory = detail::alignedSavedSizeMalloc<T>(size, alignment);
    for (std::size_t i = 0; i < size; ++i) {
        new(&memory[i]) T(std::forward<Args>(args)...);
    }
    return memory;
}

// "operator new" counterpart
template<typename T, typename... Args> T* alignedNew(std::size_t alignment, Args&&... args) {
    if constexpr (std::is_class_v<T> || sizeof...(Args) != 0) {
        return alignedConstructNew<T>(alignment, std::forward<Args>(args)...);
    } else {
        return alignedMalloc<T>(alignment);
    }
}

// "operator new[]" counterpart
template<typename T, typename... Args> T* alignedArrayNew(std::size_t size, std::size_t alignment, Args&&... args) {
    if constexpr (std::is_class_v<T> || sizeof...(Args) != 0) {
        return alignedArrayConstructNew<T>(size, alignment, std::forward<Args>(args)...);
    } else {
        return detail::alignedSavedSizeMalloc<T>(size, alignment);
    }
}

// "operator delete" counterpart
template<typename T> void alignedDelete(T* ptr) {
    if constexpr (std::is_class_v<T>) {
        ptr->~T();
    }
    alignedFree(ptr);
}

// "operator delete[]" counterpart
template<typename T> void alignedArrayDelete(T* ptr) {
    if constexpr (std::is_class_v<T>) {
        auto size = detail::alignedSavedSizeGetSize(ptr);
        for (std::size_t i = 0; i < size; ++i) {
            ptr[i].~T();
        }
    }
    detail::alignedSavedSizeFree(ptr);
}

#endif

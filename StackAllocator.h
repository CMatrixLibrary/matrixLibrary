#pragma once
#include "alignedAllocation.h"
#include <utility>

template<typename T> class StackAllocator {
    char* memory;

public:
    static const int Allignment = 32;
    static int Allign(int numberOfObjects) {
        return ((numberOfObjects*sizeof(T) + Allignment - 1) / Allignment) * Allignment;
    }

    StackAllocator(std::size_t sizeInBytes) {
        memory = reinterpret_cast<char*>(alignedMalloc(sizeInBytes, Allignment));
    }
    ~StackAllocator() {
        alignedFree(memory);
    }
    StackAllocator(const StackAllocator&) = delete;
    StackAllocator& operator=(const StackAllocator&) = delete;

    T* alloc(int objectCount) {
        T* allocatedArrayPtr = reinterpret_cast<T*>(memory);
        if constexpr (std::is_class_v<T>) {
            for (int i = 0; i < objectCount; ++i) {
                new (&allocatedArrayPtr[i]) T();
            }
        }
        memory += Allign(objectCount);
        return allocatedArrayPtr;
    }
    void dealloc(T* ptr, int objectCount) {
        if constexpr (std::is_class_v<T>) {
            for (int i = 0; i < objectCount; ++i) {
                ptr[i].~T();
            }
        }
        memory = reinterpret_cast<char*>(ptr);
    }
};
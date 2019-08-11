#pragma once
#include "alignedAllocation.h"

template<typename T> class StackAllocator {
    T* memory;

public:
    static const int Allignment = 32;
    static int Allign(int valueToAllign) {
        return ((valueToAllign + Allignment - 1) / Allignment) * Allignment;
    }

    StackAllocator(std::size_t size) {
        memory = alignedArrayNew<T>(size, Allignment);
    }
    ~StackAllocator() {
        alignedArrayDelete(memory);
    }
    StackAllocator(const StackAllocator&) = delete;
    StackAllocator& operator=(const StackAllocator&) = delete;

    T* alloc(int size) {
        T* oldMemory = memory;
        memory += Allign(size);
        return oldMemory;
    }
    void dealloc(T* ptr) {
        memory = ptr;
    }
};
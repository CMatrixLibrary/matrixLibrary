#pragma once
#include <cstdlib>
#include <new>

class StackAllocator {
    int* memory;

public:
    static const int Allignment = 32;
    static int Allign(int valueToAllign) {
        return ((valueToAllign + Allignment - 1) / Allignment) * Allignment;
    }

    StackAllocator(int size) {
        memory = reinterpret_cast<int*>(operator new[](sizeof(int)*size, static_cast<std::align_val_t>(Allignment)));
    }
    ~StackAllocator() {
        operator delete[](memory, static_cast<std::align_val_t>(Allignment));
    }
    StackAllocator(const StackAllocator&) = delete;
    StackAllocator& operator=(const StackAllocator&) = delete;

    int* alloc(int size) {
        int* oldMemory = memory;
        memory += Allign(size);
        return oldMemory;
    }
    void dealloc(int* ptr) {
        memory = ptr;
    }
};
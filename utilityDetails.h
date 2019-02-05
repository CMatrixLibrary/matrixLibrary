#pragma once
#include <cstdint>

namespace details {
    int nextPowerOf2(uint32_t n) {
        n -= 1;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        return n + 1;
    }

    bool isPowerOf2(int n) {
        return n && !(n & (n - 1));
    }
}
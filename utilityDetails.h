#pragma once
#include <cstdint>
#include <cmath>

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

    int nextPowerOf(int number, int power) {
        int nextPower = 1;
        while (nextPower < number) {
            nextPower *= power;
        }
        return nextPower;
    }
    bool isPowerOf(int number, int power) {
        double value = std::log(power) / std::log(number);
        return trunc(value) == value;
    }
}
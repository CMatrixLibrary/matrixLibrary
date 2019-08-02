#pragma once

template<typename T> class Range {
    class Iterator {
    public:
        Iterator(T value) : value(value) {}

        Iterator& operator++() {
            ++value;
            return *this;
        }
        bool operator==(const Iterator& other) {
            return value == other.value;
        }
        bool operator!=(const Iterator& other) {
            return value != other.value;
        }
        T operator*() {
            return value;
        }

    private:
        T value;
    };

public:
    Range(T start, T end) :
        start_(start),
        end_(end)
    {}
    Range(T end) :
        start_(0),
        end_(end)
    {}

    Iterator begin() {
        return Iterator(start_);
    }
    Iterator end() {
        return Iterator(end_);
    }

private:
    T start_;
    T end_;
};

template<typename Container> auto indicies(const Container& container) -> Range<decltype(container.size())> {
    return Range<decltype(container.size())>(0, container.size());
}

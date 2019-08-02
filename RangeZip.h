#pragma once
#include <tuple>
#include <cstddef>
#include <utility>

template<typename... Ts> class RangeZip {
    template <typename T> using IteratorType = decltype(((T*)nullptr)->begin());
    using TupleType = std::tuple<IteratorType<Ts>...>;

    template<typename T> using DereferenceType = decltype(**((T*)nullptr));
    using ReturnRefType = std::tuple<DereferenceType<IteratorType<Ts>>...>;

    class Iterator {
    public:
        Iterator(TupleType tuple) : values(tuple) {}
        
        void operator++() {
            increment<sizeof...(Ts)>();
        }
        bool operator!=(const Iterator& other) {
            // here better would be to check only the shortest of iterators, or assume they
            // are of the same length and only check first iterator.
            return notEqual<sizeof...(Ts)>(other);
        }
        ReturnRefType operator*() {
            return dereference(std::index_sequence_for<Ts...>());
        }

    private:
        template<std::size_t N> void increment() {
            ++std::get<N - 1>(values);
            if constexpr (N > 1) increment<N - 1>();
        }

        template<std::size_t N> bool notEqual(const Iterator& other) {
            auto result = std::get<N - 1>(values) != std::get<N - 1>(other.values);
            if constexpr (N > 1) return result && notEqual<N - 1>(other);
            else return result;
        }

        template <std::size_t... N> ReturnRefType dereference(std::index_sequence<N...>) {
            return std::tie(*std::get<N>(values)...);
        }

        TupleType values;
    };

public:
    RangeZip(Ts&... values) :
        beginTuple(createBeginTuple(values..., std::index_sequence_for<Ts...>())),
        endTuple(createEndTuple(values..., std::index_sequence_for<Ts...>()))
    {}

    Iterator begin() {
        return Iterator(beginTuple);
    }
    Iterator end() {
        return Iterator(endTuple);
    }

private:
    template <std::size_t... N> TupleType createBeginTuple(Ts&... values, std::index_sequence<N...>) {
        return std::make_tuple(values.begin()...);
    }
    template <std::size_t... N> TupleType createEndTuple(Ts&... values, std::index_sequence<N...>) {
        return std::make_tuple(values.end()...);
    }

    TupleType beginTuple;
    TupleType endTuple;
};


template<typename... Ts> RangeZip<Ts...> rangeZip(Ts&... values) {
    return RangeZip(values...);
}

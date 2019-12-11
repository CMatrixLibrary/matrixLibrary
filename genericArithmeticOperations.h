#ifndef GENERIC_ARITHMETIC_OPERATIONS_H
#define GENERIC_ARITHMETIC_OPERATIONS_H

#include <cstdlib>

/*
    meta-programming utility for creating long arithmetic expressions.
    Examples:

    calculate<Add>(a)                            >> a              // just return value
    calculate<Sub>(a)                            >> -a             // negation
    calculate<Add, Add>(a, b)                    >> a + b
    calculate<Assign, Add, Add>(a, b, c)         >> a = b + c      
    calculate<Assign, Add, Sub, Add>(a, b, c, d) >> a = b - c + d

    The seemingly unnecessery "Add"s are required to make the utility complete.
    For example you need the ability to distinguish between say "a = b" and "a = -b"
    by having calculate<Assign, Add>(a, b) and calculate<Assign, Sub>(a, b).

    But there are shortcuts added, so that "Assign" or "Add" operation are automaticly added.
    here are some examples of shorthand notation:

    calculate<Assign>(a, b) >> a = b
    calculate<Add>(a, b, c) >> a = b + c
    calculate<Sub>(a, b, c) >> a = b - c


    This utility is useless by itself as you could just type the expression instead of
    calling the "calculate" function.
    It does however allow to use it in larger functions.
    For example if you want to apply some arithmetic expression to all elements of 3 vectors
    you could create a function:

    template<OpType op1, OpType op2, OpType op3> operateOnVectors(int* a, int* b, int* C, int size) {
        for (int i = 0; i < size; ++i) {
            calculate<op1, op2, op3>(a[i], b[i], c[i]);
        }
    }

    And now you can do any permutation of operations on those 3 vectors. if you want to add vectors
    you can call operateOnVectors<Assign, Add, Add>(...).
    Or if you maybe want to copy 'c' vector to both 'a' and 'b'? (so a[i] = b[i] = c[i])
    then you call operateOnVectors<Assign, Assign, Add>(...).
    And so on.
*/

#include "always_false.h"

namespace ArithmeticOperation {
    enum OpType {
        Add,
        Sub,
        Assign,
        AddAssign,
        SubAssign
    };
}

template <ArithmeticOperation::OpType...> struct ArithmeticOperations {};

template<ArithmeticOperation::OpType op, typename T>
inline auto calculate(T&& arg) {
    if constexpr (op == ArithmeticOperation::OpType::Add) return arg;
    else if constexpr (op == ArithmeticOperation::OpType::Sub) return -arg;
    else static_assert(always_false_v<T>, "Incorrect operation");
    return arg; // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
}

template<ArithmeticOperation::OpType op, ArithmeticOperation::OpType... ops, typename T, typename... Ts>
inline auto calculate(T&& arg, Ts&&... args) {
    if constexpr (sizeof...(Ts) > sizeof...(ops)) {
        if constexpr (op == ArithmeticOperation::Add || op == ArithmeticOperation::Sub) {
            return arg = calculate<ArithmeticOperation::Add, op, ops...>(args...);
        }
        else if constexpr (op == ArithmeticOperation::Assign) return arg = calculate<ArithmeticOperation::Add, ops...>(args...);
        else if constexpr (op == ArithmeticOperation::AddAssign) return arg += calculate<ArithmeticOperation::Add, ops...>(args...);
        else return arg -= calculate<ArithmeticOperation::Add, ops...>(args...);
    }
    else {
        if constexpr (op == ArithmeticOperation::Assign) return arg = calculate<ops...>(args...);
        else if constexpr (op == ArithmeticOperation::AddAssign) return arg += calculate<ops...>(args...);
        else if constexpr (op == ArithmeticOperation::SubAssign) return arg -= calculate<ops...>(args...);
        else if constexpr (op == ArithmeticOperation::Add) return calculate<ops...>(args...) + arg;
        else return calculate<ops...>(args...) - arg;
    }
    std::abort(); // impossible to hit, but some compilerers can't see that so needed to avoid warnings 
    //return std::result_of<decltype(calculate<op, ops..., T, Ts...>(arg, args...))>::type{}; 
}

template<ArithmeticOperation::OpType... ops, typename... Ts> void operation(int n, int m, Ts&&... args) {
    for (int i = 0; i < n*m; ++i) {
        calculate<ops...>(args[i]...);
    }
}

template<ArithmeticOperation::OpType... ops, typename T, typename... Ts> void operationEffRow(int m, T&& arg, Ts&&... args) {
    for (int j = 0; j < m; ++j) {
        calculate<ops...>(arg[j], args[j]...);
    }
}
template<ArithmeticOperation::OpType... ops, typename T, typename... Ts> void operationEff(int n, int m, int effFirst, int effRest, T&& arg, Ts&&... args) {
    for (int i = 0; i < n; ++i) {
        operationEffRow<ops...>(m, &arg[i * effFirst], &args[i * effRest]...);
    }
}
template<ArithmeticOperation::OpType... ops, typename... Ts> void operationEff(int n, int m, Ts&&... args) {
    for (int i = 0; i < n; ++i) {
        operationEffRow<ops...>(m, &args.first[i * args.second]...);
    }
}
template<int n, int m, ArithmeticOperation::OpType... ops, typename... Ts> void operation(Ts&&... args) {
    for (int i = 0; i < n*m; ++i) {
        calculate<ops...>(args[i]...);
    }
}
template<int n, int m, int effFirst, int effRest, ArithmeticOperation::OpType... ops, typename T, typename... Ts> void operationEff(T&& arg, Ts&&... args) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            calculate<ops...>(arg[j + i * effFirst], args[j + i * effRest]...);
        }
    }
}
template<int n, int m, int... effs, ArithmeticOperation::OpType... ops, typename... Ts> void operationEff(ArithmeticOperations<ops...> o, Ts&&... args) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            calculate<ops...>(args[j + i * effs]...);
        }
    }
}

#endif

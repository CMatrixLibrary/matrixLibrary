#ifndef ALWAYS_FALSE_H
#define ALWAYS_FALSE_H

template<typename T> struct always_false { const static bool value = false; };
template<typename T> inline constexpr bool always_false_v = always_false<T>::value;

#endif

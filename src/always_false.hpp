#ifndef EXACT_HUBBARD_ALWAYS_FALSE_HPP
#define EXACT_HUBBARD_ALWAYS_FALSE_HPP


#include <type_traits>


template <typename>
struct alwaysFalse : std::false_type {};

template <typename T>
constexpr static bool alwaysFalse_v = alwaysFalse<T>::value;

#endif //EXACT_HUBBARD_ALWAYS_FALSE_HPP

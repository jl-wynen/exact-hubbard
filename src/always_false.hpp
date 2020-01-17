#ifndef EXACT_HUBBARD_ALWAYS_FALSE_HPP
#define EXACT_HUBBARD_ALWAYS_FALSE_HPP


#include <type_traits>


/**
 * A type that always represents `false` but defers evaluation using a template.
 * \attention Do not specialise the template!
 */
template <typename>
struct alwaysFalse : std::false_type {};

/// Helper variable for alwaysFalse.
template <typename T>
constexpr static bool alwaysFalse_v = alwaysFalse<T>::value;

#endif //EXACT_HUBBARD_ALWAYS_FALSE_HPP

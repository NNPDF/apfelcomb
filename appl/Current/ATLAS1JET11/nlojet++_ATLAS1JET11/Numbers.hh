#ifndef __NUMBERS_HH__
#define __NUMBERS_HH__

#include <limits>

// return true if value is NaN 
template<typename T>
inline bool isnan(T value) {
  return value != value;
}


// return true if value is infinite (seems to work
// with + and - infinity)
template<typename T>
inline bool isinf(T value) {
  return std::numeric_limits<T>::has_infinity &&
    value == std::numeric_limits<T>::infinity();
}

// return true if value is either NaN or infinite
template<typename T>
inline bool isbad(T value) {
  return isnan(value) || isinf(value);
}


#endif // __NUMBERS_HH__

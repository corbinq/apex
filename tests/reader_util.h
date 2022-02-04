#ifndef APEX_READER_UTIL_H
#define APEX_READER_UTIL_H

#include <string>
#include <limits>
#include <exception>
#include <cmath>
#include <fstream>
#include <limits>

/**
 * Check if a filepath exists.
 * @param name Path to file.
 * @return True if path exists.
 */
inline bool filepath_exists(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

/**
 * Extract a floating point value of type T from a given string.
 * @tparam T Type of floating point value (float, double, long double, etc.)
 * @param s String of the value
 * @return Parsed value as type T
 */
template <typename T> T extract_fp(const std::string &s) {
  T d;
  try {
    if (std::is_same<T, long double>::value) {
      d = stold(s);
    }
    else if (std::is_same<T, double>::value) {
      d = stod(s);
    }
    else if (std::is_same<T, float>::value) {
      d = stof(s);
    }
    else {
      throw std::invalid_argument("Invalid return type when extracting floating point type number from string");
    }
  }
  catch (const std::exception &e) {
    d = std::numeric_limits<T>::quiet_NaN();
  }
  return d;
}

/**
 * Generic function for testing approximate equality of two floating point values.
 * @tparam T Type of floating point value (float, double, long double, etc.)
 * @param x Value 1
 * @param y Value 2
 * @param max_rel_diff Maximum relative difference between the two values. Set by default to machine epsilon.
 * @return True if values are within acceptable difference of each other, false otherwise.
 */
template <typename T> bool approx_equal(T x, T y, T max_rel_diff = std::numeric_limits<T>::epsilon()) {
  T diff = fabs(x - y);
  x = fabs(x);
  y = fabs(y);
  T largest = (y > x) ? y : x;
  if (diff <= largest * max_rel_diff) return true;
  return false;
}

/**
 * Similar to approx_equal, but also handles the case where one or both values may be NAN.
 */
//template <typename T> bool approx_handle_nan(const T &a, const T &b) {
//  if (std::isnan(a) && std::isnan(b)) {
//    return true;
//  }
//  else if (std::isnan(a) ^ std::isnan(b)) {
//    return false;
//  }
//  else {
//    return approx_equal(a, b);
//  }
//}

#endif //APEX_READER_UTIL_H

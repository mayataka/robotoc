#ifndef ROBOTOC_PYBIND11_MACROS_HPP_
#define ROBOTOC_PYBIND11_MACROS_HPP_

#include <iostream>
#include <sstream>

namespace robotoc {

#define DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(CLASS) \
.def("clone", [](const CLASS& self) { \
  auto copy = self; \
  return copy; \
}) 

#define DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(CLASS) \
.def("__str__", [](const CLASS& self) { \
  std::stringstream ss; \
  ss << self; \
  return ss.str(); \
})

} // namespace robotoc 

#endif // ROBOTOC_PYBIND11_MACROS_HPP_
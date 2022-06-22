#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "Eigen/Core"

#include "robotoc/utils/rotation.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(rotation, m) {
  m.def("omp_get_thread_num", &omp_get_thread_num);
}

} // namespace python
} // namespace robotoc
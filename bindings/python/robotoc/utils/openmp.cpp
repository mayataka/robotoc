#include <pybind11/pybind11.h>

#include <omp.h>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(openmp, m) {
  m.def("omp_get_thread_num", &omp_get_thread_num);
}

} // namespace python
} // namespace robotoc
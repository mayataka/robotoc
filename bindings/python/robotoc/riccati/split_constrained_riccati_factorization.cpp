#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_constrained_riccati_factorization, m) {
  py::class_<SplitConstrainedRiccatiFactorization>(m, "SplitConstrainedRiccatiFactorization")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    .def("__str__", [](const SplitConstrainedRiccatiFactorization& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc
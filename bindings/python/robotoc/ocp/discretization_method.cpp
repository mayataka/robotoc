#include <pybind11/pybind11.h>

#include "robotoc/ocp/discretization_method.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(discretization_method, m) {
  py::enum_<DiscretizationMethod>(m, "DiscretizationMethod", py::arithmetic())
    .value("GridBased",  DiscretizationMethod::GridBased)
    .value("PhaseBased", DiscretizationMethod::PhaseBased)
    .export_values();
}

} // namespace python
} // namespace robotoc
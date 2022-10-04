#include <pybind11/pybind11.h>

#include "robotoc/solver/interpolation_order.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(interpolation_order, m) {
  py::enum_<InterpolationOrder>(m, "InterpolationOrder", py::arithmetic())
    .value("Linear",  InterpolationOrder::Linear)
    .value("Zero", InterpolationOrder::Zero)
    .export_values();
}

} // namespace python
} // namespace robotoc
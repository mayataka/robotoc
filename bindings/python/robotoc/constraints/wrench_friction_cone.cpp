#include <pybind11/pybind11.h>

#include "robotoc/constraints/wrench_friction_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(wrench_friction_cone, m) {
  py::class_<WrenchFrictionCone, ConstraintComponentBase, 
             std::shared_ptr<WrenchFrictionCone>>(m, "WrenchFrictionCone")
    .def(py::init<const Robot&, const double, const double>(),
          py::arg("robot"), py::arg("X"), py::arg("Y"))
    .def("set_rectangular", &WrenchFrictionCone::setRectangular,
          py::arg("X"), py::arg("Y"));
}

} // namespace python
} // namespace robotoc
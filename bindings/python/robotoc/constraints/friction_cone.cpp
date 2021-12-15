#include <pybind11/pybind11.h>

#include "robotoc/constraints/friction_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(friction_cone, m) {
  py::class_<FrictionCone, ConstraintComponentBase, 
             std::shared_ptr<FrictionCone>>(m, "FrictionCone")
    .def(py::init<const Robot&, const double>(),
          py::arg("robot"), py::arg("mu"))
    .def("set_friction_coefficient", &FrictionCone::setFrictionCoefficient,
          py::arg("mu"));
}

} // namespace python
} // namespace robotoc
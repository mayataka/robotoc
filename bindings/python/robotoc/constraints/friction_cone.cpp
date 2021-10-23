#include <pybind11/pybind11.h>

#include "robotoc/constraints/friction_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(friction_cone, m) {
  py::class_<FrictionCone, ConstraintComponentBase, 
             std::shared_ptr<FrictionCone>>(m, "FrictionCone")
    .def(py::init<const Robot&, const double, const double, const double>(),
         py::arg("robot"), py::arg("mu"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rule")=0.995)
    .def("set_friction_coefficient", &FrictionCone::setFrictionCoefficient);
}

} // namespace python
} // namespace robotoc
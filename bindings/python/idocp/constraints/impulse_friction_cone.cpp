#include <pybind11/pybind11.h>

#include "idocp/constraints/impulse_friction_cone.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_friction_cone, m) {
  py::class_<ImpulseFrictionCone, ImpulseConstraintComponentBase, 
             std::shared_ptr<ImpulseFrictionCone>>(m, "ImpulseFrictionCone")
    .def(py::init<const Robot&, const double, const double, const double>(),
         py::arg("robot"), py::arg("mu"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rule")=0.995)
    .def("set_friction_coefficient", &ImpulseFrictionCone::setFrictionCoefficient);
}

} // namespace python
} // namespace idocp
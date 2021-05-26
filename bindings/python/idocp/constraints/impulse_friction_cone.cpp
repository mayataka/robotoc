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
         py::arg("fraction_to_boundary_rate")=0.995)
    .def("set_friction_coefficient", &ImpulseFrictionCone::setFrictionCoefficient);

  m.def("create_impulse_friction_cone", [](const Robot& robot, const double mu) {
    return std::make_shared<ImpulseFrictionCone>(robot, mu);
  });
}

} // namespace python
} // namespace idocp
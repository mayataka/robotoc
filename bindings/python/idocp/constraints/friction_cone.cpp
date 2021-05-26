#include <pybind11/pybind11.h>

#include "idocp/constraints/friction_cone.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(friction_cone, m) {
  py::class_<FrictionCone, ConstraintComponentBase, 
             std::shared_ptr<FrictionCone>>(m, "FrictionCone")
    .def(py::init<const Robot&, const double, const double, const double>(),
         py::arg("robot"), py::arg("mu"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rate")=0.995)
    .def("set_friction_coefficient", &FrictionCone::setFrictionCoefficient);

  m.def("create_friction_cone", [](const Robot& robot, const double mu) {
    return std::make_shared<FrictionCone>(robot, mu);
  });
}

} // namespace python
} // namespace idocp
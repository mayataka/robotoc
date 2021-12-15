#include <pybind11/pybind11.h>

#include "robotoc/constraints/impulse_friction_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_friction_cone, m) {
  py::class_<ImpulseFrictionCone, ImpulseConstraintComponentBase, 
             std::shared_ptr<ImpulseFrictionCone>>(m, "ImpulseFrictionCone")
    .def(py::init<const Robot&, const double>(),
          py::arg("robot"), py::arg("mu"))
    .def("set_friction_coefficient", &ImpulseFrictionCone::setFrictionCoefficient,
          py::arg("mu"));
}

} // namespace python
} // namespace robotoc
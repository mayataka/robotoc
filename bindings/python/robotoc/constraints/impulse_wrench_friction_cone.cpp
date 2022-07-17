#include <pybind11/pybind11.h>

#include "robotoc/constraints/impulse_wrench_friction_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_wrench_friction_cone, m) {
  py::class_<ImpulseWrenchFrictionCone, ImpulseConstraintComponentBase, 
             std::shared_ptr<ImpulseWrenchFrictionCone>>(m, "ImpulseWrenchFrictionCone")
    .def(py::init<const Robot&, const double, const double, const double>(),
          py::arg("robot"), py::arg("mu"), py::arg("X"), py::arg("Y"))
    .def("clone", [](const ImpulseWrenchFrictionCone& self) {
       auto other = self;
       return other;
     })
    .def("clone", &ImpulseWrenchFrictionCone::clone)
    .def("set_friction_coefficient", &ImpulseWrenchFrictionCone::setFrictionCoefficient,
          py::arg("mu"))
    .def("set_rectangular", &ImpulseWrenchFrictionCone::setRectangular,
          py::arg("X"), py::arg("Y"));
}

} // namespace python
} // namespace robotoc
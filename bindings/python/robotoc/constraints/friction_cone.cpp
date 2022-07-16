#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/constraints/friction_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(friction_cone, m) {
  py::class_<FrictionCone, ConstraintComponentBase, 
             std::shared_ptr<FrictionCone>>(m, "FrictionCone")
    .def(py::init<const Robot&, const std::vector<double>&>(),
          py::arg("robot"), py::arg("mu"))
    .def(py::init<const Robot&, const double>(),
          py::arg("robot"), py::arg("mu"))
    .def("clone", [](const FrictionCone& self) {
       auto other = self;
       return other;
     })
    .def("set_friction_coefficient", static_cast<void (FrictionCone::*)(const std::vector<double>&)>(&FrictionCone::setFrictionCoefficient),
          py::arg("mu"))
    .def("set_friction_coefficient", static_cast<void (FrictionCone::*)(const double)>(&FrictionCone::setFrictionCoefficient),
          py::arg("mu"));
}

} // namespace python
} // namespace robotoc
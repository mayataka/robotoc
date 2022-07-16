#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/robot_properties.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(robot_properties, m) {
  py::class_<RobotProperties>(m, "RobotProperties")
    .def(py::init<>())
    .def("clone", [](const RobotProperties& self) {
       auto other = self;
       return other;
     })
    .def_readwrite("generalized_momentum_bias", &RobotProperties::generalized_momentum_bias);
}

} // namespace python
} // namespace robotoc
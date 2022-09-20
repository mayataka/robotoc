#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/robot_properties.hpp"
#include "robotoc/utils/pybind11_macros.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(robot_properties, m) {
  py::class_<RobotProperties>(m, "RobotProperties")
    .def(py::init<>())
    .def_readwrite("generalized_momentum_bias", &RobotProperties::generalized_momentum_bias)
    .def_readwrite("has_generalized_momentum_bias", &RobotProperties::has_generalized_momentum_bias)
     DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(RobotProperties);
}

} // namespace python
} // namespace robotoc
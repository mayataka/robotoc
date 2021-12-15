#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_position_lower_limit.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_position_lower_limit, m) {
  py::class_<JointPositionLowerLimit, ConstraintComponentBase, 
             std::shared_ptr<JointPositionLowerLimit>>(m, "JointPositionLowerLimit")
    .def(py::init<const Robot&>(),
         py::arg("robot"));
}

} // namespace python
} // namespace robotoc
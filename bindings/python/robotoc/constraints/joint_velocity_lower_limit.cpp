#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_velocity_lower_limit.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_velocity_lower_limit, m) {
  py::class_<JointVelocityLowerLimit, ConstraintComponentBase, 
             std::shared_ptr<JointVelocityLowerLimit>>(m, "JointVelocityLowerLimit")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    .def("clone", [](const JointVelocityLowerLimit& self) {
       auto other = self;
       return other;
     });
}

} // namespace python
} // namespace robotoc
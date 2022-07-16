#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_velocity_upper_limit.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_velocity_upper_limit, m) {
  py::class_<JointVelocityUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointVelocityUpperLimit>>(m, "JointVelocityUpperLimit")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    .def("clone", [](const JointVelocityUpperLimit& self) {
       auto other = self;
       return other;
     });
}

} // namespace python
} // namespace robotoc
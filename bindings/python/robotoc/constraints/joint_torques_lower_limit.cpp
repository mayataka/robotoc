#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_torques_lower_limit.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_torques_lower_limit, m) {
  py::class_<JointTorquesLowerLimit, ConstraintComponentBase, 
             std::shared_ptr<JointTorquesLowerLimit>>(m, "JointTorquesLowerLimit")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    .def("clone", &JointTorquesLowerLimit::clone);
}

} // namespace python
} // namespace robotoc
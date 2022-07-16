#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_acceleration_lower_limit.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_acceleration_lower_limit, m) {
  py::class_<JointAccelerationLowerLimit, ConstraintComponentBase, 
             std::shared_ptr<JointAccelerationLowerLimit>>(m, "JointAccelerationLowerLimit")
    .def(py::init<const Robot&, const Eigen::VectorXd&>(),
         py::arg("robot"), py::arg("amin"))
    .def("clone", [](const JointAccelerationLowerLimit& self) {
       auto other = self;
       return other;
     });
}

} // namespace python
} // namespace robotoc
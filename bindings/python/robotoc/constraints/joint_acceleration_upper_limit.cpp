#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_acceleration_upper_limit.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_acceleration_upper_limit, m) {
  py::class_<JointAccelerationUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointAccelerationUpperLimit>>(m, "JointAccelerationUpperLimit")
    .def(py::init<const Robot&, const Eigen::VectorXd&>(),
         py::arg("robot"), py::arg("amax"))
    .def("clone", [](const JointAccelerationUpperLimit& self) {
       auto other = self;
       return other;
     });
}

} // namespace python
} // namespace robotoc
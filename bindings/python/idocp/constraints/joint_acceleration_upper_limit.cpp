#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/constraints/joint_acceleration_upper_limit.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_acceleration_upper_limit, m) {
  py::class_<JointAccelerationUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointAccelerationUpperLimit>>(m, "JointAccelerationUpperLimit")
    .def(py::init<const Robot&, const Eigen::VectorXd&, const double, const double>(),
         py::arg("robot"), py::arg("amax"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rate")=0.995);

  m.def("create_joint_acceleration_upper_limit", [](const Robot& robot, 
                                                    const Eigen::VectorXd& amax) {
    return std::make_shared<JointAccelerationUpperLimit>(robot, amax);
  });
}

} // namespace python
} // namespace idocp
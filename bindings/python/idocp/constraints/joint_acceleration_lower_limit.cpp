#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/constraints/joint_acceleration_lower_limit.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_acceleration_lower_limit, m) {
  py::class_<JointAccelerationLowerLimit, ConstraintComponentBase, 
             std::shared_ptr<JointAccelerationLowerLimit>>(m, "JointAccelerationLowerLimit")
    .def(py::init<const Robot&, const Eigen::VectorXd&, const double, const double>(),
         py::arg("robot"), py::arg("amin"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rule")=0.995);

  m.def("create_joint_acceleration_lower_limit", [](const Robot& robot, 
                                                    const Eigen::VectorXd& amin) {
    return std::make_shared<JointAccelerationLowerLimit>(robot, amin);
  });
}

} // namespace python
} // namespace idocp
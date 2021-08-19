#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/constraints/joint_velocity_upper_limit.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_velocity_upper_limit, m) {
  py::class_<JointVelocityUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointVelocityUpperLimit>>(m, "JointVelocityUpperLimit")
    .def(py::init<const Robot&, const double, const double>(),
         py::arg("robot"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rule")=0.995);

  m.def("create_joint_velocity_upper_limit", [](const Robot& robot) {
    return std::make_shared<JointVelocityUpperLimit>(robot);
  });
}

} // namespace python
} // namespace idocp
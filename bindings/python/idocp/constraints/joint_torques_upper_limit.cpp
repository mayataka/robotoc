#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/constraints/joint_torques_upper_limit.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_torques_upper_limit, m) {
  py::class_<JointTorquesUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointTorquesUpperLimit>>(m, "JointTorquesUpperLimit")
    .def(py::init<const Robot&, const double, const double>(),
         py::arg("robot"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rule")=0.995);

  m.def("create_joint_torques_upper_limit", [](const Robot& robot) {
    return std::make_shared<JointTorquesUpperLimit>(robot);
  });
}

} // namespace python
} // namespace idocp
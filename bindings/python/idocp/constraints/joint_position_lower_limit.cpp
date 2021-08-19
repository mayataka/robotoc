#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/constraints/joint_position_lower_limit.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_position_lower_limit, m) {
  py::class_<JointPositionLowerLimit, ConstraintComponentBase, 
             std::shared_ptr<JointPositionLowerLimit>>(m, "JointPositionLowerLimit")
    .def(py::init<const Robot&, const double, const double>(),
         py::arg("robot"), py::arg("barrier")=1.0e-04,
         py::arg("fraction_to_boundary_rule")=0.995);

  m.def("create_joint_position_lower_limit", [](const Robot& robot) {
    return std::make_shared<JointPositionLowerLimit>(robot);
  });
}

} // namespace python
} // namespace idocp
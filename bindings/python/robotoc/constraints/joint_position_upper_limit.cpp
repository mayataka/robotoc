#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_position_upper_limit, m) {
  py::class_<JointPositionUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointPositionUpperLimit>>(m, "JointPositionUpperLimit")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(JointPositionUpperLimit);
}

} // namespace python
} // namespace robotoc
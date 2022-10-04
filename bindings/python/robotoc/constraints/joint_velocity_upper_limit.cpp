#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_velocity_upper_limit, m) {
  py::class_<JointVelocityUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointVelocityUpperLimit>>(m, "JointVelocityUpperLimit")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(JointVelocityUpperLimit);
}

} // namespace python
} // namespace robotoc
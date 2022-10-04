#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(joint_torques_upper_limit, m) {
  py::class_<JointTorquesUpperLimit, ConstraintComponentBase, 
             std::shared_ptr<JointTorquesUpperLimit>>(m, "JointTorquesUpperLimit")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(JointTorquesUpperLimit);
}

} // namespace python
} // namespace robotoc
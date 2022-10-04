#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/task_space_3d_cost.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(task_space_3d_cost, m) {
  py::class_<TaskSpace3DCost, CostFunctionComponentBase,
             std::shared_ptr<TaskSpace3DCost>>(m, "TaskSpace3DCost")
    .def(py::init<const Robot&, const int>(),
          py::arg("robot"), py::arg("frame_id"))
    .def(py::init<const Robot&, const std::string&>(),
          py::arg("robot"), py::arg("frame_name"))
    .def(py::init<const Robot&, const int, const std::shared_ptr<TaskSpace3DRefBase>&>(),
          py::arg("robot"), py::arg("frame_id"), py::arg("ref"))
    .def(py::init<const Robot&, const int, const Eigen::Vector3d&>(),
          py::arg("robot"), py::arg("frame_id"), py::arg("const_ref"))
    .def(py::init<const Robot&, const std::string&, const std::shared_ptr<TaskSpace3DRefBase>&>(),
          py::arg("robot"), py::arg("frame_name"), py::arg("ref"))
    .def(py::init<const Robot&, const std::string&, const Eigen::Vector3d&>(),
          py::arg("robot"), py::arg("frame_name"), py::arg("const_ref"))
    .def(py::init<>())
    .def("set_ref", &TaskSpace3DCost::set_ref,
          py::arg("ref"))
    .def("set_const_ref", &TaskSpace3DCost::set_const_ref,
          py::arg("const_ref"))
    .def("set_weight", &TaskSpace3DCost::set_weight,
          py::arg("weight"))
    .def("set_weight_terminal", &TaskSpace3DCost::set_weight_terminal,
          py::arg("weight_terminal"))
    .def("set_weight_impact", &TaskSpace3DCost::set_weight_impact,
          py::arg("weight_impact"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(TaskSpace3DCost);
}

} // namespace python
} // namespace robotoc
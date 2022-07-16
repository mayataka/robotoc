#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/task_space_6d_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(task_space_6d_cost, m) {
  py::class_<TaskSpace6DCost, CostFunctionComponentBase,
             std::shared_ptr<TaskSpace6DCost>>(m, "TaskSpace6DCost")
    .def(py::init<const Robot&, const int>(),
          py::arg("robot"), py::arg("frame_id"))
    .def(py::init<const Robot&, const std::string&>(),
          py::arg("robot"), py::arg("frame_name"))
    .def(py::init<const Robot&, const int, const std::shared_ptr<TaskSpace6DRefBase>&>(),
          py::arg("robot"), py::arg("frame_id"), py::arg("ref"))
    .def(py::init<const Robot&, const int, const SE3&>(),
          py::arg("robot"), py::arg("frame_id"), py::arg("const_ref"))
    .def(py::init<const Robot&, const int, const Eigen::Vector3d&, const Eigen::Matrix3d&>(),
          py::arg("robot"), py::arg("frame_id"), py::arg("const_position_ref"), 
          py::arg("const_rotation_ref"))
    .def(py::init<const Robot&, const std::string&, const SE3&>(),
          py::arg("robot"), py::arg("frame_name"), py::arg("const_ref"))
    .def(py::init<const Robot&, const std::string&, const Eigen::Vector3d&, const Eigen::Matrix3d&>(),
          py::arg("robot"), py::arg("frame_name"), py::arg("const_position_ref"), 
          py::arg("const_rotation_ref"))
    .def(py::init<>())
    .def("clone", [](const TaskSpace6DCost& self) {
       auto other = self;
       return other;
     })
    .def("set_ref", &TaskSpace6DCost::set_ref,
          py::arg("ref"))
    .def("set_const_ref", static_cast<void (TaskSpace6DCost::*)(const SE3&)>(&TaskSpace6DCost::set_const_ref),
          py::arg("const_ref"))
    .def("set_const_ref", static_cast<void (TaskSpace6DCost::*)(const Eigen::Vector3d&, const Eigen::Matrix3d&)>(&TaskSpace6DCost::set_const_ref),
          py::arg("const_position_ref"), py::arg("const_rotation_ref"))
    .def("set_weight", &TaskSpace6DCost::set_weight,
          py::arg("weight_position"), py::arg("weight_rotation"))
    .def("set_weight_terminal", &TaskSpace6DCost::set_weight_terminal,
          py::arg("weight_position_terminal"), 
          py::arg("weight_rotation_terminal"))
    .def("set_weight_impulse", &TaskSpace6DCost::set_weight_impulse,
          py::arg("weight_position_impulse"), 
          py::arg("weight_rotation_impulse"));
}

} // namespace python
} // namespace robotoc
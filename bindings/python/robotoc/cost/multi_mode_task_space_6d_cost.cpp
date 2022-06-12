#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/cost/multi_mode_task_space_6d_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(multi_mode_task_space_6d_cost, m) {
  py::class_<MultiModeTaskSpace6DCost, CostFunctionComponentBase,
             std::shared_ptr<MultiModeTaskSpace6DCost>>(m, "MultiModeTaskSpace6DCost")
    .def(py::init<const Robot&, const int>(),
          py::arg("robot"), py::arg("frame_id"))
    .def(py::init<const Robot&, const std::string&>(),
          py::arg("robot"), py::arg("frame_name"))
    .def("set_x6d_ref", 
          static_cast<void (MultiModeTaskSpace6DCost::*)(const Eigen::Vector3d&, const Eigen::Matrix3d&, const int)>(&MultiModeTaskSpace6DCost::set_x6d_ref),
          py::arg("trans_ref"), py::arg("rot_ref"), py::arg("contact_mode_id")=0)
    .def("set_x6d_ref", 
          static_cast<void (MultiModeTaskSpace6DCost::*)(const Eigen::Vector3d&, const Eigen::Matrix3d&, const std::vector<int>&)>(&MultiModeTaskSpace6DCost::set_x6d_ref),
          py::arg("trans_ref"), py::arg("rot_ref"), py::arg("contact_mode_ids"))
    .def("set_x6d_weight", 
          static_cast<void (MultiModeTaskSpace6DCost::*)(const Eigen::Vector3d&, const Eigen::Vector3d&, const int)>(&MultiModeTaskSpace6DCost::set_x6d_weight),
          py::arg("trans_weight"), py::arg("rot_weight"), py::arg("contact_mode_id")=0)
    .def("set_x6d_weight", 
          static_cast<void (MultiModeTaskSpace6DCost::*)(const Eigen::Vector3d&, const Eigen::Vector3d&, const std::vector<int>&)>(&MultiModeTaskSpace6DCost::set_x6d_weight),
          py::arg("trans_weight"), py::arg("rot_weight"), py::arg("contact_mode_ids"))
    .def("set_x6df_weight", &MultiModeTaskSpace6DCost::set_x6df_weight, 
          py::arg("trans_weight"), py::arg("rot_weight"))
    .def("set_x6di_weight", &MultiModeTaskSpace6DCost::set_x6di_weight, 
          py::arg("trans_weight"), py::arg("rot_weight"));
}

} // namespace python
} // namespace robotoc
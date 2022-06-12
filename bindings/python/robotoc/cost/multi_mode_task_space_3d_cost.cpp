#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/cost/multi_mode_task_space_3d_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(multi_mode_task_space_3d_cost, m) {
  py::class_<MultiModeTaskSpace3DCost, CostFunctionComponentBase,
             std::shared_ptr<MultiModeTaskSpace3DCost>>(m, "MultiModeTaskSpace3DCost")
    .def(py::init<const Robot&, const int>(),
          py::arg("robot"), py::arg("frame_id"))
    .def(py::init<const Robot&, const std::string&>(),
          py::arg("robot"), py::arg("frame_name"))
    .def("set_x3d_ref", 
          static_cast<void (MultiModeTaskSpace3DCost::*)(const Eigen::Vector3d&, const int)>(&MultiModeTaskSpace3DCost::set_x3d_ref),
          py::arg("x3d_ref"), py::arg("contact_mode_id")=0)
    .def("set_x3d_ref", 
          static_cast<void (MultiModeTaskSpace3DCost::*)(const Eigen::Vector3d&, const std::vector<int>&)>(&MultiModeTaskSpace3DCost::set_x3d_ref),
          py::arg("x3d_ref"), py::arg("contact_mode_ids"))
    .def("set_x3d_weight", 
          static_cast<void (MultiModeTaskSpace3DCost::*)(const Eigen::Vector3d&, const int)>(&MultiModeTaskSpace3DCost::set_x3d_weight),
          py::arg("x3d_weight"), py::arg("contact_mode_id")=0)
    .def("set_x3d_weight", 
          static_cast<void (MultiModeTaskSpace3DCost::*)(const Eigen::Vector3d&, const std::vector<int>&)>(&MultiModeTaskSpace3DCost::set_x3d_weight),
          py::arg("x3d_weight"), py::arg("contact_mode_ids"))
    .def("set_x3df_weight", &MultiModeTaskSpace3DCost::set_x3df_weight, 
          py::arg("x3df_weight"))
    .def("set_x3di_weight", &MultiModeTaskSpace3DCost::set_x3di_weight, 
          py::arg("x3di_weight"));
}

} // namespace python
} // namespace robotoc
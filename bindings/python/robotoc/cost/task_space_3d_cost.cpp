#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/task_space_3d_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(task_space_3d_cost, m) {
  py::class_<TaskSpace3DCost, CostFunctionComponentBase,
             std::shared_ptr<TaskSpace3DCost>>(m, "TaskSpace3DCost")
    .def(py::init<const Robot&, const int>(),
          py::arg("robot"), py::arg("frame_id"))
    .def("set_x3d_ref", &TaskSpace3DCost::set_x3d_ref,
          py::arg("x3d_ref"))
    .def("set_x3d_weight", &TaskSpace3DCost::set_x3d_weight, 
          py::arg("x3d_weight"))
    .def("set_x3df_weight", &TaskSpace3DCost::set_x3df_weight, 
          py::arg("x3df_weight"))
    .def("set_x3di_weight", &TaskSpace3DCost::set_x3di_weight, 
          py::arg("x3di_weight"));
}

} // namespace python
} // namespace robotoc
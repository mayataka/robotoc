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
    .def("set_x6d_ref", &TaskSpace6DCost::set_x6d_ref,
          py::arg("trans_ref"), py::arg("rot_ref"))
    .def("set_x6d_weight", &TaskSpace6DCost::set_x6d_weight,
          py::arg("trans_weight"), py::arg("rot_weight"))
    .def("set_x6df_weight", &TaskSpace6DCost::set_x6df_weight,
          py::arg("trans_weight"), py::arg("rot_weight"))
    .def("set_x6di_weight", &TaskSpace6DCost::set_x6di_weight,
          py::arg("trans_weight"), py::arg("rot_weight"));
}

} // namespace python
} // namespace robotoc
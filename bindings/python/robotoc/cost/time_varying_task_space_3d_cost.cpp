#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_varying_task_space_3d_cost, m) {
  py::class_<TimeVaryingTaskSpace3DCost, CostFunctionComponentBase,
             std::shared_ptr<TimeVaryingTaskSpace3DCost>>(m, "TimeVaryingTaskSpace3DCost")
    .def(py::init<const Robot&, const int, 
                  const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>&>(),
          py::arg("robot"), py::arg("frame_id"), py::arg("x3d_ref"))
    .def("set_x3d_ref", &TimeVaryingTaskSpace3DCost::set_x3d_ref,
          py::arg("x3d_ref"))
    .def("set_x3d_weight", &TimeVaryingTaskSpace3DCost::set_x3d_weight,
          py::arg("x3d_weight"))
    .def("set_x3df_weight", &TimeVaryingTaskSpace3DCost::set_x3df_weight,
          py::arg("x3df_weight"))
    .def("set_x3di_weight", &TimeVaryingTaskSpace3DCost::set_x3di_weight,
          py::arg("x3di_weight"));
}

} // namespace python
} // namespace robotoc
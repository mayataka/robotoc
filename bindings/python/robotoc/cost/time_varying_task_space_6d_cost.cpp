#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/time_varying_task_space_6d_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_varying_task_space_6d_cost, m) {
  py::class_<TimeVaryingTaskSpace6DCost, CostFunctionComponentBase,
             std::shared_ptr<TimeVaryingTaskSpace6DCost>>(m, "TimeVaryingTaskSpace6DCost")
    .def(py::init<const Robot&, const int, 
                  const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>&>(),
          py::arg("robot"), py::arg("frame_id"), py::arg("x6d_ref"))
    .def(py::init<const Robot&, const std::string&, 
                  const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>&>(),
          py::arg("robot"), py::arg("frame_name"), py::arg("x6d_ref"))
    .def("set_x6d_ref", &TimeVaryingTaskSpace6DCost::set_x6d_ref,
          py::arg("x6d_ref"))
    .def("set_x6d_weight", &TimeVaryingTaskSpace6DCost::set_x6d_weight,
          py::arg("trans_weight"), py::arg("rot_weight"))
    .def("set_x6df_weight", &TimeVaryingTaskSpace6DCost::set_x6df_weight,
          py::arg("trans_weight"), py::arg("rot_weight"))
    .def("set_x6di_weight", &TimeVaryingTaskSpace6DCost::set_x6di_weight,
          py::arg("trans_weight"), py::arg("rot_weight"));
}

} // namespace python
} // namespace robotoc
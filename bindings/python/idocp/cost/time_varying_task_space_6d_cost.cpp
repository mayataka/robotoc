#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/cost/time_varying_task_space_6d_cost.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_varying_task_space_6d_cost, m) {
  py::class_<TimeVaryingTaskSpace6DCost, CostFunctionComponentBase,
             std::shared_ptr<TimeVaryingTaskSpace6DCost>>(m, "TimeVaryingTaskSpace6DCost")
    .def(py::init<const Robot&, const int, 
                  const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>&>())
    .def("set_ref", &TimeVaryingTaskSpace6DCost::set_ref)
    .def("set_q_weight", &TimeVaryingTaskSpace6DCost::set_q_weight)
    .def("set_qf_weight", &TimeVaryingTaskSpace6DCost::set_qf_weight)
    .def("set_qi_weight", &TimeVaryingTaskSpace6DCost::set_qi_weight);

  m.def("create_time_varying_task_space_6d_cost", [](const Robot& robot, 
                                                     const int frame_id,
                                                     const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& ref) {
    return std::make_shared<TimeVaryingTaskSpace6DCost>(robot, frame_id, ref);
  });
}

} // namespace python
} // namespace idocp
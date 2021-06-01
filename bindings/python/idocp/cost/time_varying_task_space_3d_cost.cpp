#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/cost/time_varying_task_space_3d_cost.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_varying_task_space_3d_cost, m) {
  py::class_<TimeVaryingTaskSpace3DCost, CostFunctionComponentBase,
             std::shared_ptr<TimeVaryingTaskSpace3DCost>>(m, "TimeVaryingTaskSpace3DCost")
    .def(py::init<const Robot&, const int, 
                  const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>&>())
    .def("set_ref", &TimeVaryingTaskSpace3DCost::set_ref)
    .def("set_q_weight", &TimeVaryingTaskSpace3DCost::set_q_weight)
    .def("set_qf_weight", &TimeVaryingTaskSpace3DCost::set_qf_weight)
    .def("set_qi_weight", &TimeVaryingTaskSpace3DCost::set_qi_weight);

  m.def("create_time_varying_task_space_3d_cost", [](const Robot& robot, 
                                                     const int frame_id,
                                                     const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>& ref) {
    return std::make_shared<TimeVaryingTaskSpace3DCost>(robot, frame_id, ref);
  });
}

} // namespace python
} // namespace idocp
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/cost/task_space_6d_cost.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(task_space_6d_cost, m) {
  py::class_<TaskSpace6DCost, CostFunctionComponentBase,
             std::shared_ptr<TaskSpace6DCost>>(m, "TaskSpace6DCost")
    .def(py::init<const Robot&, const int>())
    .def("set_q_6d_ref", &TaskSpace6DCost::set_q_6d_ref)
    .def("set_q_weight", &TaskSpace6DCost::set_q_weight)
    .def("set_qf_weight", &TaskSpace6DCost::set_qf_weight)
    .def("set_qi_weight", &TaskSpace6DCost::set_qi_weight);
}

} // namespace python
} // namespace idocp
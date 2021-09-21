#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/cost/time_varying_com_cost.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_varying_com_cost, m) {
  py::class_<TimeVaryingCoMCost, 
             CostFunctionComponentBase,
             std::shared_ptr<TimeVaryingCoMCost>>(m, "TimeVaryingCoMCost")
    .def(py::init<const Robot&, 
                  const std::shared_ptr<TimeVaryingCoMRefBase>&>())
    .def("set_ref", &TimeVaryingCoMCost::set_ref)
    .def("set_q_weight", &TimeVaryingCoMCost::set_q_weight)
    .def("set_qf_weight", &TimeVaryingCoMCost::set_qf_weight)
    .def("set_qi_weight", &TimeVaryingCoMCost::set_qi_weight);
}

} // namespace python
} // namespace idocp
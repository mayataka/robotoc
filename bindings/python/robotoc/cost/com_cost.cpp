#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/com_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(com_cost, m) {
  py::class_<CoMCost, CostFunctionComponentBase,
             std::shared_ptr<CoMCost>>(m, "CoMCost")
    .def(py::init<const Robot&>())
    .def("set_com_ref", &CoMCost::set_CoM_ref)
    .def("set_q_weight", &CoMCost::set_q_weight)
    .def("set_qf_weight", &CoMCost::set_qf_weight)
    .def("set_qi_weight", &CoMCost::set_qi_weight);
}

} // namespace python
} // namespace robotoc
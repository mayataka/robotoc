#include <pybind11/pybind11.h>

#include "robotoc/cost/cost_function.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(cost_function, m) {
  py::class_<CostFunction, std::shared_ptr<CostFunction>>(m, "CostFunction")
    .def(py::init<const double, const double>(),
          py::arg("discount_factor"), py::arg("discount_time_step"))
    .def(py::init<>())
    .def("clone", [](const CostFunction& self) {
       auto other = self;
       return other;
     })
    .def("set_discount_factor", &CostFunction::setDiscountFactor,
          py::arg("discount_factor"), py::arg("discount_time_step"))
    .def("discount_factor", &CostFunction::discountFactor)
    .def("discount_time_step", &CostFunction::discountTimeStep)
    .def("push_back", &CostFunction::push_back)
    .def("clear", &CostFunction::clear)
    .def("create_cost_function_data", &CostFunction::createCostFunctionData,
          py::arg("robot"))
    .def("eval_stage_cost", &CostFunction::evalStageCost,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), 
          py::arg("grid_info"), py::arg("s"))
    .def("linearize_stage_cost", &CostFunction::linearizeStageCost,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), 
          py::arg("grid_info"), py::arg("s"), py::arg("kkt_residual"))
    .def("quadratize_stage_cost", &CostFunction::quadratizeStageCost,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), 
          py::arg("grid_info"), py::arg("s"), py::arg("kkt_residual"), 
          py::arg("kkt_matrix"))
    .def("eval_terminal_cost", &CostFunction::evalTerminalCost,
          py::arg("robot"), py::arg("data"), py::arg("grid_info"), py::arg("s"))
    .def("linearize_terminal_cost", &CostFunction::linearizeTerminalCost,
          py::arg("robot"), py::arg("data"), py::arg("grid_info"), py::arg("s"),
          py::arg("kkt_residual"))
    .def("quadratize_terminal_cost", &CostFunction::quadratizeTerminalCost,
          py::arg("robot"), py::arg("data"), py::arg("grid_info"), py::arg("s"),
          py::arg("kkt_residual"), py::arg("kkt_hessian"))
    .def("eval_impulse_cost", &CostFunction::evalImpulseCost,
          py::arg("robot"), py::arg("impulse_status"), py::arg("data"), 
          py::arg("grid_info"), py::arg("s"))
    .def("linearize_impulse_cost", &CostFunction::linearizeImpulseCost,
          py::arg("robot"), py::arg("impulse_status"), py::arg("data"), 
          py::arg("grid_info"), py::arg("s"), py::arg("kkt_residual"))
    .def("quadratize_impulse_cost", &CostFunction::quadratizeImpulseCost,
          py::arg("robot"), py::arg("impulse_status"), py::arg("data"), 
          py::arg("grid_info"), py::arg("s"), py::arg("kkt_residual"), 
          py::arg("kkt_matrix"));
}

} // namespace python
} // namespace robotoc
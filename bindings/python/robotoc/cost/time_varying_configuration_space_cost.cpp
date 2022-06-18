#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/time_varying_configuration_space_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_varying_configuration_space_cost, m) {
  py::class_<TimeVaryingConfigurationSpaceCost, 
             CostFunctionComponentBase,
             std::shared_ptr<TimeVaryingConfigurationSpaceCost>>(m, "TimeVaryingConfigurationSpaceCost")
    .def(py::init<const Robot&, 
                  const std::shared_ptr<TimeVaryingConfigurationRefBase>&>())
    .def("set_q_ref", &TimeVaryingConfigurationSpaceCost::set_q_ref,
          py::arg("q_ref"))
    .def("set_q_weight", &TimeVaryingConfigurationSpaceCost::set_q_weight,
          py::arg("q_weight"))
    .def("set_q_weight_terminal", &TimeVaryingConfigurationSpaceCost::set_q_weight_terminal,
          py::arg("q_weight_terminal"))
    .def("set_q_weight_impulse", &TimeVaryingConfigurationSpaceCost::set_q_weight_impulse,
          py::arg("q_weight_impulse"));
}

} // namespace python
} // namespace robotoc
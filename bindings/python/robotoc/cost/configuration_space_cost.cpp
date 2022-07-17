#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/configuration_space_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(configuration_space_cost, m) {
  py::class_<ConfigurationSpaceCost, CostFunctionComponentBase,
             std::shared_ptr<ConfigurationSpaceCost>>(m, "ConfigurationSpaceCost")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    .def(py::init<const Robot&, const std::shared_ptr<ConfigurationSpaceRefBase>&>(),
          py::arg("robot"), py::arg("ref"))
    .def("clone", &ConfigurationSpaceCost::clone)
    .def("set_ref", &ConfigurationSpaceCost::set_ref,
          py::arg("ref"))
    .def("set_q_ref", &ConfigurationSpaceCost::set_q_ref,
          py::arg("q_ref"))
    .def("set_v_ref", &ConfigurationSpaceCost::set_v_ref,
          py::arg("v_ref"))
    .def("set_u_ref", &ConfigurationSpaceCost::set_u_ref,
          py::arg("u_ref"))
    .def("set_q_weight", &ConfigurationSpaceCost::set_q_weight,
          py::arg("q_weight"))
    .def("set_v_weight", &ConfigurationSpaceCost::set_v_weight,
          py::arg("v_weight"))
    .def("set_a_weight", &ConfigurationSpaceCost::set_a_weight,
          py::arg("a_weight"))
    .def("set_u_weight", &ConfigurationSpaceCost::set_u_weight,
          py::arg("u_weight"))
    .def("set_q_weight_terminal", &ConfigurationSpaceCost::set_q_weight_terminal,
          py::arg("q_weight_terminal"))
    .def("set_v_weight_terminal", &ConfigurationSpaceCost::set_v_weight_terminal,
          py::arg("v_weight_terminal"))
    .def("set_q_weight_impulse", &ConfigurationSpaceCost::set_q_weight_impulse,
          py::arg("q_weight_impulse"))
    .def("set_v_weight_impulse", &ConfigurationSpaceCost::set_v_weight_impulse,
          py::arg("v_weight_impulse"))
    .def("set_dv_weight_impulse", &ConfigurationSpaceCost::set_dv_weight_impulse,
          py::arg("dv_weight_impulse"));
}

} // namespace python
} // namespace robotoc
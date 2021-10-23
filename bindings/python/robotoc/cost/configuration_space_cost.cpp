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
    .def(py::init<const Robot&>())
    .def("set_q_ref", &ConfigurationSpaceCost::set_q_ref)
    .def("set_v_ref", &ConfigurationSpaceCost::set_v_ref)
    .def("set_u_ref", &ConfigurationSpaceCost::set_u_ref)
    .def("set_q_weight", &ConfigurationSpaceCost::set_q_weight)
    .def("set_v_weight", &ConfigurationSpaceCost::set_v_weight)
    .def("set_a_weight", &ConfigurationSpaceCost::set_a_weight)
    .def("set_u_weight", &ConfigurationSpaceCost::set_u_weight)
    .def("set_qf_weight", &ConfigurationSpaceCost::set_qf_weight)
    .def("set_vf_weight", &ConfigurationSpaceCost::set_vf_weight)
    .def("set_qi_weight", &ConfigurationSpaceCost::set_qi_weight)
    .def("set_vi_weight", &ConfigurationSpaceCost::set_vi_weight)
    .def("set_dvi_weight", &ConfigurationSpaceCost::set_dvi_weight);
}

} // namespace python
} // namespace robotoc
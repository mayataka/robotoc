#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/cost/multi_mode_configuration_space_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(multi_mode_configuration_space_cost, m) {
  py::class_<MultiModeConfigurationSpaceCost, CostFunctionComponentBase,
             std::shared_ptr<MultiModeConfigurationSpaceCost>>(m, "MultiModeConfigurationSpaceCost")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    .def("set_q_ref", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const int)>(&MultiModeConfigurationSpaceCost::set_q_ref),
          py::arg("q_ref"), py::arg("contact_mode_id")=0)
    .def("set_q_ref", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const std::vector<int>&)>(&MultiModeConfigurationSpaceCost::set_q_ref),
          py::arg("q_ref"), py::arg("contact_mode_ids"))
    .def("set_v_ref", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const int)>(&MultiModeConfigurationSpaceCost::set_v_ref),
          py::arg("v_ref"), py::arg("contact_mode_id")=0)
    .def("set_v_ref", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const std::vector<int>&)>(&MultiModeConfigurationSpaceCost::set_v_ref),
          py::arg("v_ref"), py::arg("contact_mode_ids"))
    .def("set_u_ref", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const int)>(&MultiModeConfigurationSpaceCost::set_u_ref),
          py::arg("u_ref"), py::arg("contact_mode_id")=0)
    .def("set_u_ref", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const std::vector<int>&)>(&MultiModeConfigurationSpaceCost::set_u_ref),
          py::arg("u_ref"), py::arg("contact_mode_ids"))
    .def("set_q_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const int)>(&MultiModeConfigurationSpaceCost::set_q_weight),
          py::arg("q_weight"), py::arg("contact_mode_id")=0)
    .def("set_q_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const std::vector<int>&)>(&MultiModeConfigurationSpaceCost::set_q_weight),
          py::arg("set_q_weight"), py::arg("contact_mode_ids"))
    .def("set_v_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const int)>(&MultiModeConfigurationSpaceCost::set_v_weight),
          py::arg("v_weight"), py::arg("contact_mode_id")=0)
    .def("set_v_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const std::vector<int>&)>(&MultiModeConfigurationSpaceCost::set_v_weight),
          py::arg("set_v_weight"), py::arg("contact_mode_ids"))
    .def("set_a_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const int)>(&MultiModeConfigurationSpaceCost::set_a_weight),
          py::arg("a_weight"), py::arg("contact_mode_id")=0)
    .def("set_a_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const std::vector<int>&)>(&MultiModeConfigurationSpaceCost::set_a_weight),
          py::arg("set_a_weight"), py::arg("contact_mode_ids"))
    .def("set_u_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const int)>(&MultiModeConfigurationSpaceCost::set_u_weight),
          py::arg("u_weight"), py::arg("contact_mode_id")=0)
    .def("set_u_weight", 
          static_cast<void (MultiModeConfigurationSpaceCost::*)(const Eigen::VectorXd&, const std::vector<int>&)>(&MultiModeConfigurationSpaceCost::set_u_weight),
          py::arg("set_u_weight"), py::arg("contact_mode_ids"))
    .def("set_q_weight_terminal", &MultiModeConfigurationSpaceCost::set_q_weight_terminal,
          py::arg("q_weight_terminal"))
    .def("set_v_weight_terminal", &MultiModeConfigurationSpaceCost::set_v_weight_terminal,
          py::arg("v_weight_terminal"))
    .def("set_q_weight_impulse", &MultiModeConfigurationSpaceCost::set_q_weight_impulse,
          py::arg("q_weight_impulse"))
    .def("set_v_weight_impulse", &MultiModeConfigurationSpaceCost::set_v_weight_impulse,
          py::arg("v_weight_impulse"))
    .def("set_dv_weight_impulse", &MultiModeConfigurationSpaceCost::set_dv_weight_impulse,
          py::arg("dv_weight_impulse"));
}

} // namespace python
} // namespace robotoc
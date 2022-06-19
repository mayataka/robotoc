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
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    .def(py::init<const Robot&, const std::shared_ptr<CoMRefBase>&>(),
          py::arg("robot"), py::arg("ref"))
    .def(py::init<const Robot&, const Eigen::Vector3d&>(),
          py::arg("robot"), py::arg("const_ref"))
    .def(py::init<>())
    .def("set_ref", &CoMCost::set_ref,
          py::arg("ref"))
    .def("set_const_ref", &CoMCost::set_const_ref,
          py::arg("const_ref"))
    .def("set_weight", &CoMCost::set_weight,
          py::arg("weight"))
    .def("set_weight_terminal", &CoMCost::set_weight_terminal,
          py::arg("weight_terminal"))
    .def("set_weight_impulse", &CoMCost::set_weight_impulse,
          py::arg("weight_impulse"));
}

} // namespace python
} // namespace robotoc
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
    .def("set_com_ref", &CoMCost::set_com_ref,
          py::arg("com_ref"))
    .def("set_com_weight", &CoMCost::set_com_weight,
          py::arg("com_weight"))
    .def("set_comf_weight", &CoMCost::set_comf_weight,
          py::arg("comf_weight"))
    .def("set_comi_weight", &CoMCost::set_comi_weight,
          py::arg("comi_weight"));
}

} // namespace python
} // namespace robotoc
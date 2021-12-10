#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/time_varying_com_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_varying_com_cost, m) {
  py::class_<TimeVaryingCoMCost, 
             CostFunctionComponentBase,
             std::shared_ptr<TimeVaryingCoMCost>>(m, "TimeVaryingCoMCost")
    .def(py::init<const Robot&, 
                  const std::shared_ptr<TimeVaryingCoMRefBase>&>(),
          py::arg("robot"), py::arg("com_ref"))
    .def("set_com_ref", &TimeVaryingCoMCost::set_com_ref,
          py::arg("com_ref"))
    .def("set_com_weight", &TimeVaryingCoMCost::set_com_weight,
          py::arg("com_weight"))
    .def("set_comf_weight", &TimeVaryingCoMCost::set_comf_weight,
          py::arg("comf_weight"))
    .def("set_comi_weight", &TimeVaryingCoMCost::set_comi_weight,
          py::arg("comi_weight"));
}

} // namespace python
} // namespace robotoc
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/sto/sto_cost_function.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(sto_cost_function, m) {
  py::class_<STOCostFunction, std::shared_ptr<STOCostFunction>>(m, "STOCostFunction")
    .def(py::init<>())
    .def("exist", &STOCostFunction::exist,
          py::arg("name"))
    .def("add", &STOCostFunction::add,
          py::arg("name"), py::arg("cost"))
    .def("erase", &STOCostFunction::erase,
          py::arg("name"))
    .def("get", &STOCostFunction::get,
          py::arg("name"))
    .def("clear", &STOCostFunction::clear)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(STOCostFunction);
}

} // namespace python
} // namespace robotoc
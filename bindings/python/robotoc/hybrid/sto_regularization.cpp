#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/hybrid/sto_regularization.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(sto_regularization, m) {
  py::enum_<STORegularizationType>(m, "STORegularizationType", py::arithmetic())
    .value("Const",  STORegularizationType::Const)
    .value("Abs", STORegularizationType::Abs)
    .value("Quad", STORegularizationType::Quad)
    .value("Exp", STORegularizationType::Exp)
    .value("Exp2", STORegularizationType::Exp2)
    .value("None", STORegularizationType::None)
    .export_values();

  py::class_<STORegularization>(m, "STORegularization")
    .def(py::init<const STORegularizationType&, const double>(),
         py::arg("reg_type"), py::arg("w"))
    .def(py::init(&STORegularization::defaultSTORegularization))
    .def("set_regularization", &STORegularization::setRegularization);
}

} // namespace python
} // namespace robotoc 
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/riccati/riccati_factorization.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(riccati_factorization, m) {
  py::class_<RiccatiFactorization>(m, "RiccatiFactorization")
    .def(py::init<const Robot&, const int, const int>(),
          py::arg("robot"), py::arg("N"), py::arg("max_num_each_discrete_events"))
    .def(py::init<>())
    .def("__getitem__", [](const RiccatiFactorization& self, const int i) {
        return self.data[i];
      })
    .def("__setitem__", [](RiccatiFactorization& self, const int i, 
                           const SplitRiccatiFactorization& riccati) {
        self.data[i] = riccati;
      })
    .def_readwrite("aux", &RiccatiFactorization::aux)
    .def_readwrite("lift", &RiccatiFactorization::lift)
    .def_readwrite("impulse", &RiccatiFactorization::impulse)
    .def("__str__", [](const RiccatiFactorization& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc
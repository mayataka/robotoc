#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/riccati/split_riccati_factorization.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_riccati_factorization, m) {
  py::class_<SplitRiccatiFactorization>(m, "SplitRiccatiFactorization")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    .def_readwrite("P", &SplitRiccatiFactorization::P)
    .def_readwrite("s", &SplitRiccatiFactorization::s)
    .def_readwrite("psi_x", &SplitRiccatiFactorization::psi_x)
    .def_readwrite("psi_u", &SplitRiccatiFactorization::psi_u)
    .def_readwrite("Psi", &SplitRiccatiFactorization::Psi)
    .def_readwrite("phi_x", &SplitRiccatiFactorization::phi_x)
    .def_readwrite("phi_u", &SplitRiccatiFactorization::phi_u)
    .def_readwrite("Phi", &SplitRiccatiFactorization::Phi)
    .def_readwrite("xi", &SplitRiccatiFactorization::xi)
    .def_readwrite("chi", &SplitRiccatiFactorization::chi)
    .def_readwrite("rho", &SplitRiccatiFactorization::rho)
    .def_readwrite("eta", &SplitRiccatiFactorization::eta)
    .def_readwrite("iota", &SplitRiccatiFactorization::iota);
}

} // namespace python
} // namespace robotoc
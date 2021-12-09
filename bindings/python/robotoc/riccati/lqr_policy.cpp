#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/riccati/lqr_policy.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(lqr_policy, m) {
  py::class_<LQRPolicy>(m, "LQRPolicy")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    .def_readwrite("K", &LQRPolicy::K)
    .def_readwrite("k", &LQRPolicy::k)
    .def_readwrite("T", &LQRPolicy::T)
    .def_readwrite("W", &LQRPolicy::W)
    .def("Kq", &LQRPolicy::Kq)
    .def("Kv", &LQRPolicy::Kv);
}

} // namespace python
} // namespace robotoc
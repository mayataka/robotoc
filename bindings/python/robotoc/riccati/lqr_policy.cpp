#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/hybrid/hybrid_container.hpp"


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

  py::class_<hybrid_container<LQRPolicy>>(m, "LQRPolicies")
    .def(py::init<const Robot&, const int, const int>(),
          py::arg("robot"), py::arg("N"), py::arg("max_num_each_discrete_events"))
    .def(py::init<>())
    .def("__getitem__", [](const hybrid_container<LQRPolicy>& self, const int i) {
        return self.data[i];
      })
    .def("__setitem__", [](hybrid_container<LQRPolicy>& self, const int i, 
                           const LQRPolicy& lqr_policy) {
        self.data[i] = lqr_policy;
      })
    .def_readwrite("aux", &hybrid_container<LQRPolicy>::aux)
    .def_readwrite("lift", &hybrid_container<LQRPolicy>::lift);
}

} // namespace python
} // namespace robotoc
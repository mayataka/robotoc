#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/trotting_swing_foot_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(trotting_swing_foot_ref, m) {
  py::class_<TrottingSwingFootRef, SwingFootRefBase,
             std::shared_ptr<TrottingSwingFootRef>>(m, "TrottingSwingFootRef")
    .def(py::init<const int, const int, const double, const double>())
    .def("update_q_3d_ref", &TrottingSwingFootRef::update_q_3d_ref);
}

} // namespace python
} // namespace robotoc
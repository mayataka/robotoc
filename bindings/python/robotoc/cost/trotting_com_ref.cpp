#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/trotting_com_ref.hpp" 


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(trotting_com_ref, m) {
  py::class_<TrottingCoMRef, TimeVaryingCoMRefBase,
             std::shared_ptr<TrottingCoMRef>>(m, "TrottingCoMRef")
    .def(py::init<const Eigen::Vector3d&, const double, const int>(),
          py::arg("com0"), py::arg("step_length"), py::arg("max_steps")=100)
    .def("next_step", &TrottingCoMRef::nextStep)
    .def("set_contact_sequence", &TrottingCoMRef::setContactSequence)
    .def("set_initial_phase_rate", &TrottingCoMRef::setInitialPhaseRate)
    .def("set_last_phase_rate", &TrottingCoMRef::setLastPhaseRate)
    .def("update_com_ref", &TrottingCoMRef::update_com_ref)
    .def("is_active", &TrottingCoMRef::isActive);
}

} // namespace python
} // namespace robotoc
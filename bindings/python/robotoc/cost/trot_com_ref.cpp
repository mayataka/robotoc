#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/trot_com_ref.hpp" 


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(trot_com_ref, m) {
  py::class_<TrotCoMRef, CoMRefBase,
             std::shared_ptr<TrotCoMRef>>(m, "TrotCoMRef")
    .def(py::init<const Eigen::Vector3d&, const double, const int>(),
          py::arg("com0"), py::arg("step_length"), py::arg("max_steps")=100)
    .def("next_step", &TrotCoMRef::nextStep)
    .def("set_contact_sequence", &TrotCoMRef::setContactSequence)
    .def("set_initial_phase_rate", &TrotCoMRef::setInitialPhaseRate)
    .def("set_last_phase_rate", &TrotCoMRef::setLastPhaseRate)
    .def("updateRef", &TrotCoMRef::updateRef)
    .def("is_active", &TrotCoMRef::isActive);
}

} // namespace python
} // namespace robotoc
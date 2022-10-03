#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/ocp/discretization_method.hpp"
#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_discretization, m) {
  py::enum_<DiscretizationMethod>(m, "DiscretizationMethod", py::arithmetic())
    .value("GridBased",  DiscretizationMethod::GridBased)
    .value("PhaseBased", DiscretizationMethod::PhaseBased)
    .export_values();

  py::class_<TimeDiscretization>(m, "TimeDiscretization")
    .def(py::init<const double, const int, const int>(), 
          py::arg("T"), py::arg("N"), py::arg("reserved_num_discrete_events")=0)
    .def("N", &TimeDiscretization::N)
    .def("size", &TimeDiscretization::size) 
    .def("grid", &TimeDiscretization::grid,
          py::arg("i")) 
    .def("front", &TimeDiscretization::front)
    .def("back", &TimeDiscretization::back)
    .def("__getitem__", [](const TimeDiscretization& self, const int i) {
        return self[i];
     })
    .def("max_time_step", &TimeDiscretization::maxTimeStep)
    .def("discretize", &TimeDiscretization::discretize,
          py::arg("contact_sequence"), py::arg("t")) 
    .def("correct_time_steps", &TimeDiscretization::correctTimeSteps,
          py::arg("contact_sequence"), py::arg("t")) 
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(TimeDiscretization)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(TimeDiscretization);
}

} // namespace python
} // namespace robotoc
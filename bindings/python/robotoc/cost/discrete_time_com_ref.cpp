#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/cost/discrete_time_com_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(discrete_time_com_ref, m) {
  py::class_<DiscreteTimeCoMRef, TimeVaryingCoMRefBase,
             std::shared_ptr<DiscreteTimeCoMRef>>(m, "DiscreteTimeCoMRef")
    .def(py::init<const std::vector<Eigen::Vector3d>&>(), 
          py::arg("com_position_to_contact_position"))
    .def("set_com_ref", 
          static_cast<void (DiscreteTimeCoMRef::*)(const std::shared_ptr<ContactSequence>&)>(&DiscreteTimeCoMRef::setCoMRef),
          py::arg("contact_sequence"))
    .def("set_com_ref", 
          static_cast<void (DiscreteTimeCoMRef::*)(const std::shared_ptr<ContactSequence>&,
                                                   const Eigen::Vector3d&, const Eigen::Vector3d&, 
                                                   const double, const double)>(&DiscreteTimeCoMRef::setCoMRef),
          py::arg("contact_sequence"), py::arg("first_com_ref"), py::arg("last_com_ref"),
          py::arg("first_rate")=0, py::arg("last_com_ref")=0)
    .def("update_com_ref", &DiscreteTimeCoMRef::update_com_ref,
          py::arg("grid_info"), py::arg("com_ref"))
    .def("is_active", &DiscreteTimeCoMRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
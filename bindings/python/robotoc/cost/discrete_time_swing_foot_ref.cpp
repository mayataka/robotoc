#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/discrete_time_swing_foot_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(discrete_time_swing_foot_ref, m) {
  py::class_<DiscreteTimeSwingFootRef, TimeVaryingTaskSpace3DRefBase,
             std::shared_ptr<DiscreteTimeSwingFootRef>>(m, "DiscreteTimeSwingFootRef")
    .def(py::init<const int, const double>(), 
          py::arg("contact_index"), py::arg("swing_height"))
    .def("set_swing_foot_ref", 
          static_cast<void (DiscreteTimeSwingFootRef::*)(const std::shared_ptr<ContactSequence>&)>(&DiscreteTimeSwingFootRef::setSwingFootRef),
          py::arg("contact_sequence"))
    .def("set_swing_foot_ref", 
          static_cast<void (DiscreteTimeSwingFootRef::*)(const std::shared_ptr<ContactSequence>&,
                                                         const Eigen::Vector3d&, const Eigen::Vector3d&, 
                                                         const double, const double)>(&DiscreteTimeSwingFootRef::setSwingFootRef),
          py::arg("contact_sequence"), py::arg("first_contact_position"), py::arg("last_contact_position"),
          py::arg("first_rate")=0, py::arg("last_com_ref")=0)
    .def("update_x3d_ref", &DiscreteTimeSwingFootRef::update_x3d_ref,
          py::arg("grid_info"), py::arg("x3d_ref"))
    .def("is_active", &DiscreteTimeSwingFootRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
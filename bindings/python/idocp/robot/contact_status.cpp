#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/robot/contact_status.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_status, m) {
  py::class_<ContactStatus>(m, "ContactStatus")
    .def(py::init<const int>())
    .def("is_contact_active", 
          static_cast<bool (ContactStatus::*)(const int) const>(&ContactStatus::isContactActive))
    .def("is_contact_active", 
          static_cast<const std::vector<bool>& (ContactStatus::*)() const>(&ContactStatus::isContactActive))
    .def("activate_contact", &ContactStatus::activateContact)
    .def("deactivate_contact", &ContactStatus::deactivateContact)
    .def("activate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<int>& contact_indices)>(&ContactStatus::activateContacts))
    .def("deactivate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<int>& contact_indices)>(&ContactStatus::deactivateContacts))
    .def("activate_contacts", 
          static_cast<void (ContactStatus::*)()>(&ContactStatus::activateContacts))
    .def("deactivate_contacts", 
          static_cast<void (ContactStatus::*)()>(&ContactStatus::deactivateContacts))
    .def("contact_point", &ContactStatus::contactPoint)
    .def("contact_points", &ContactStatus::contactPoints);
}

} // namespace python
} // namespace idocp
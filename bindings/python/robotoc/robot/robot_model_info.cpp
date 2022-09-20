#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/robot_model_info.hpp"
#include "robotoc/utils/pybind11_macros.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(robot_model_info, m) {
  py::enum_<BaseJointType>(m, "BaseJointType", py::arithmetic())
    .value("FixedBase", BaseJointType::FixedBase)
    .value("FloatingBase", BaseJointType::FloatingBase)
    .export_values();

  py::class_<RobotModelInfo>(m, "RobotModelInfo")
    .def(py::init<const std::string&, const BaseJointType, 
                  const std::vector<ContactModelInfo>&, 
                  const std::vector<ContactModelInfo>&, const double>(),
          py::arg("urdf_path"), py::arg("base_joint_type"), 
          py::arg("point_contacts"), py::arg("surface_contacts"),
          py::arg("contact_inv_damping"))
    .def(py::init<>())
    .def_readwrite("urdf_path", &RobotModelInfo::urdf_path)
    .def_readwrite("base_joint_type", &RobotModelInfo::base_joint_type)
    .def_readwrite("point_contacts", &RobotModelInfo::point_contacts)
    .def_readwrite("surface_contacts", &RobotModelInfo::surface_contacts)
    .def_readwrite("contact_inv_damping", &RobotModelInfo::contact_inv_damping)
    .def_static("Manipulator", &RobotModelInfo::Manipulator,
                 py::arg("urdf_path"))
    .def_static("Quadruped", &RobotModelInfo::Quadruped,
                 py::arg("urdf_path"), py::arg("point_contacts"))
    .def_static("Humanoid", &RobotModelInfo::Humanoid,
                 py::arg("urdf_path"), py::arg("surface_contacts"))
     DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(RobotModelInfo);
}

} // namespace python
} // namespace robotoc
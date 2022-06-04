#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/robot.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(robot, m) {
  py::enum_<BaseJointType>(m, "BaseJointType", py::arithmetic())
    .value("FixedBase", BaseJointType::FixedBase)
    .value("FloatingBase", BaseJointType::FloatingBase)
    .export_values();

  py::class_<Robot>(m, "Robot")
    .def(py::init<const std::string&, const BaseJointType&>(), 
          py::arg("path_to_urdf"),
          py::arg("base_joint_type")=BaseJointType::FixedBase)
    .def(py::init<const std::string&, const BaseJointType&, 
                  const std::vector<int>&, const std::vector<ContactType>&, 
                  const std::pair<double, double>&, const double>(),
          py::arg("path_to_urdf"), py::arg("base_joint_type"),
          py::arg("contact_frames"), py::arg("contact_types"), 
          py::arg("baumgarte_weights"), py::arg("contact_inv_damping")=0.)
    .def(py::init<const std::string&, const BaseJointType&, 
                  const std::vector<std::string>&, const std::vector<ContactType>&, 
                  const std::pair<double, double>&, const double>(),
          py::arg("path_to_urdf"), py::arg("base_joint_type"),
          py::arg("contact_frame_names"), py::arg("contact_types"), 
          py::arg("baumgarte_weights"), py::arg("contact_inv_damping")=0.)
    .def(py::init<const std::string&, const BaseJointType&, 
                  const std::vector<int>&, const std::vector<ContactType>&, 
                  const double, const double>(),
          py::arg("path_to_urdf"), py::arg("base_joint_type"),
          py::arg("contact_frames"), py::arg("contact_types"), 
          py::arg("baumgarte_time_step"), py::arg("contact_inv_damping")=0.)
    .def(py::init<const std::string&, const BaseJointType&, 
                  const std::vector<std::string>&, const std::vector<ContactType>&, 
                  const double, const double>(),
          py::arg("path_to_urdf"), py::arg("base_joint_type"),
          py::arg("contact_frame_names"), py::arg("contact_types"), 
          py::arg("baumgarte_time_step"), py::arg("contact_inv_damping")=0.)
    .def("integrate_coeff_wise_jacobian", [](const Robot& self, const Eigen::VectorXd& q, Eigen::MatrixXd& J) {
        self.integrateCoeffWiseJacobian(q, J);
      },  py::arg("q"), py::arg("J"))
    .def("forward_kinematics", [](Robot& self, const Eigen::VectorXd& q) {
        self.updateFrameKinematics(q);
      },  py::arg("q"))
    .def("frame_position", &Robot::framePosition,
          py::arg("frame_id"))
    .def("frame_rotation", &Robot::frameRotation,
          py::arg("frame_id"))
    .def("frame_placement", &Robot::framePlacement,
          py::arg("frame_id"))
    .def("com", &Robot::CoM)
    .def("transform_from_local_to_world", [](const Robot& self, 
                                             const int frame_id, 
                                             const Eigen::Vector3d& vec_local,
                                             Eigen::Vector3d& vec_world) {
        self.transformFromLocalToWorld(frame_id, vec_local, vec_world);
      },  py::arg("frame_id"), py::arg("vec_local"), py::arg("vec_world"))
    .def("generate_feasible_configuration", &Robot::generateFeasibleConfiguration)
    .def("normalize_configuration", [](const Robot& self, Eigen::VectorXd& q) {
        self.normalizeConfiguration(q);
      },  py::arg("q"))
    .def("create_contact_status", &Robot::createContactStatus)
    .def("create_impulse_status", &Robot::createImpulseStatus)
    .def("frame_id", &Robot::frameId,
          py::arg("frame_name"))
    .def("total_weight", &Robot::totalWeight)
    .def("dimq", &Robot::dimq)
    .def("dimv", &Robot::dimv)
    .def("dimu", &Robot::dimu)
    .def("max_dimf", &Robot::max_dimf)
    .def("dim_passive", &Robot::dim_passive)
    .def("max_num_contacts", &Robot::maxNumContacts)
    .def("contact_type", &Robot::contactType,
          py::arg("contact_index"))
    .def("contact_types", &Robot::contactTypes)
    .def("contact_frames", &Robot::contactFrames)
    .def("point_contact_frames", &Robot::pointContactFrames)
    .def("surface_contact_frames", &Robot::surfaceContactFrames)
    .def("set_joint_effort_limit", &Robot::setJointEffortLimit,
          py::arg("joint_effort_limit"))
    .def("set_joint_velocity_limit", &Robot::setJointVelocityLimit,
          py::arg("joint_velocity_limit"))
    .def("set_lower_joint_position_limit", &Robot::setLowerJointPositionLimit,
          py::arg("lower_joint_position_limit"))
    .def("set_upper_joint_position_limit", &Robot::setUpperJointPositionLimit,
          py::arg("upper_joint_position_limit"))
    .def("__str__", [](const Robot& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc
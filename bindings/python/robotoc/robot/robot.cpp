#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/robot.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

std::vector<int> empty_contact_frames = {};
std::pair<double, double> zero_baumgarte_weights = {0., 0.};

PYBIND11_MODULE(robot, m) {
  py::enum_<BaseJointType>(m, "BaseJointType", py::arithmetic())
    .value("FixedBase",  BaseJointType::FixedBase)
    .value("FloatingBase", BaseJointType::FloatingBase)
    .export_values();

  py::class_<Robot>(m, "Robot")
    .def(py::init<const std::string&, const BaseJointType&, 
                  const std::vector<int>&, 
                  const std::pair<double, double>&>(),
         py::arg("path_to_urdf"),
         py::arg("base_joint_type")=BaseJointType::FixedBase,
         py::arg("contact_frames")=empty_contact_frames,
         py::arg("baumgarte_weights")=zero_baumgarte_weights)
    .def(py::init<const std::string&, const BaseJointType&, 
                  const std::vector<int>&, const double>())
    .def("forward_kinematics", [](Robot& self, const Eigen::VectorXd& q) {
        self.updateFrameKinematics(q);
      })
    .def("frame_position", &Robot::framePosition)
    .def("frame_rotation", &Robot::frameRotation)
    .def("frame_placement", &Robot::framePlacement)
    .def("com", &Robot::CoM)
    .def("transform_from_local_to_world", [](const Robot& self, 
                                             const int frame_id, 
                                             const Eigen::Vector3d& vec_local,
                                             Eigen::Vector3d& vec_world) {
        self.transformFromLocalToWorld(frame_id, vec_local, vec_world);
      })
    .def("generate_feasible_configuration", &Robot::generateFeasibleConfiguration)
    .def("normalize_configuration", [](const Robot& self, Eigen::VectorXd& q) {
        self.normalizeConfiguration(q);
      })
    .def("create_contact_status", &Robot::createContactStatus)
    .def("create_impulse_status", &Robot::createImpulseStatus)
    .def("total_weight", &Robot::totalWeight)
    .def("dimq", &Robot::dimq)
    .def("dimv", &Robot::dimv)
    .def("dimu", &Robot::dimu)
    .def("max_dimf", &Robot::max_dimf)
    .def("dim_passive", &Robot::dim_passive)
    .def("max_point_contacts", &Robot::maxPointContacts)
    .def("contact_frames", &Robot::contactFrames)
    .def("set_joint_effort_limit", &Robot::setJointEffortLimit)
    .def("set_joint_velocity_limit", &Robot::setJointVelocityLimit)
    .def("set_lower_joint_position_limit", &Robot::setLowerJointPositionLimit)
    .def("set_upper_joint_position_limit", &Robot::setUpperJointPositionLimit)
    .def("__str__", [](const Robot& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc
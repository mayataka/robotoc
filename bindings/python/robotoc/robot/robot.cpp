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
    .def("clone", &Robot::clone)
    .def("integrate_configuration", [](const Robot& self, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, const double dt) {
        Eigen::VectorXd q_ret = Eigen::VectorXd::Zero(self.dimq());
        self.integrateConfiguration(q, v, dt, q_ret);
        return q_ret;
     }, py::arg("q"), py::arg("v"), py::arg("dt"))
    .def("subtract_configuration", [](const Robot& self, const Eigen::VectorXd& qf, 
                                      const Eigen::VectorXd& q0) {
        Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(self.dimv());
        self.subtractConfiguration(qf, q0, qdiff);
        return qdiff;
     }, py::arg("qf"), py::arg("q0"))
    .def("integrate_coeff_wise_jacobian", [](const Robot& self, const Eigen::VectorXd& q) {
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(self.dimq(), self.dimv());
        self.integrateCoeffWiseJacobian(q, J);
        return J;
      },  py::arg("q"))
    .def("forward_kinematics", [](Robot& self, const Eigen::VectorXd& q) {
        self.updateFrameKinematics(q);
      },  py::arg("q"))
    .def("forward_kinematics", [](Robot& self, const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& v) {
        self.updateFrameKinematics(q, v);
      },  py::arg("q"), py::arg("v"))
    .def("forward_kinematics", [](Robot& self, const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& v, 
                                  const Eigen::VectorXd& a) {
        self.updateFrameKinematics(q, v, a);
      },  py::arg("q"), py::arg("v"), py::arg("a"))
    .def("frame_position", 
          static_cast<const Eigen::Vector3d& (Robot::*)(const int) const>(&Robot::framePosition),
          py::arg("frame_id"))
    .def("frame_position", 
          static_cast<const Eigen::Vector3d& (Robot::*)(const std::string&) const>(&Robot::framePosition),
          py::arg("frame_name"))
    .def("frame_rotation", 
          static_cast<const Eigen::Matrix3d& (Robot::*)(const int) const>(&Robot::frameRotation),
          py::arg("frame_id"))
    .def("frame_rotation", 
          static_cast<const Eigen::Matrix3d& (Robot::*)(const std::string&) const>(&Robot::frameRotation),
          py::arg("frame_name"))
    .def("frame_placement", 
          static_cast<const SE3& (Robot::*)(const int) const>(&Robot::framePlacement),
          py::arg("frame_id"))
    .def("frame_placement", 
          static_cast<const SE3& (Robot::*)(const std::string&) const>(&Robot::framePlacement),
          py::arg("frame_name"))
    .def("com", &Robot::CoM)
    .def("frame_linear_velocity", [](const Robot& self, const int frame_id) {
        return self.frameLinearVelocity(frame_id);
      }, py::arg("frame_id"))
    .def("frame_linear_velocity", [](const Robot& self, const std::string& frame_name) {
        return self.frameLinearVelocity(frame_name);
      }, py::arg("frame_name"))
    .def("frame_angular_velocity", [](const Robot& self, const int frame_id) {
        return self.frameAngularVelocity(frame_id);
      }, py::arg("frame_id"))
    .def("frame_angular_velocity", [](const Robot& self, const std::string& frame_name) {
        return self.frameAngularVelocity(frame_name);
      }, py::arg("frame_name"))
    .def("frame_spatial_velocity", [](const Robot& self, const int frame_id) {
        return self.frameSpatialVelocity(frame_id);
      }, py::arg("frame_id"))
    .def("frame_spatial_velocity", [](const Robot& self, const std::string& frame_name) {
        return self.frameSpatialVelocity(frame_name);
      }, py::arg("frame_name"))
    .def("com_velocity", &Robot::CoMVelocity)
    .def("transform_from_local_to_world", [](const Robot& self, 
                                             const int frame_id, 
                                             const Eigen::Vector3d& vec_local) {
        Eigen::Vector3d vec_world = Eigen::Vector3d::Zero();
        self.transformFromLocalToWorld(frame_id, vec_local, vec_world);
        return vec_world;
      },  py::arg("frame_id"), py::arg("vec_local"))
    .def("generate_feasible_configuration", &Robot::generateFeasibleConfiguration)
    .def("normalize_configuration", [](const Robot& self, Eigen::VectorXd& q) {
        self.normalizeConfiguration(q);
      },  py::arg("q"))
    .def("create_contact_status", &Robot::createContactStatus)
    .def("create_impulse_status", &Robot::createImpulseStatus)
    .def("frame_id", &Robot::frameId,
          py::arg("frame_name"))
    .def("frame_name", &Robot::frameName,
          py::arg("frame_id"))
    .def("total_weight", &Robot::totalWeight)
    .def("dimq", &Robot::dimq)
    .def("dimv", &Robot::dimv)
    .def("dimu", &Robot::dimu)
    .def("max_dimf", &Robot::max_dimf)
    .def("dim_passive", &Robot::dim_passive)
    .def("has_floating_base", &Robot::hasFloatingBase)
    .def("max_num_contacts", &Robot::maxNumContacts)
    .def("max_num_point_contacts", &Robot::maxNumPointContacts)
    .def("max_num_surface_contacts", &Robot::maxNumSurfaceContacts)
    .def("contact_type", &Robot::contactType,
          py::arg("contact_index"))
    .def("contact_types", &Robot::contactTypes)
    .def("contact_frames", &Robot::contactFrames)
    .def("contact_frame_names", &Robot::contactFrameNames)
    .def("point_contact_frames", &Robot::pointContactFrames)
    .def("point_contact_frame_names", &Robot::pointContactFrameNames)
    .def("surface_contact_frames", &Robot::surfaceContactFrames)
    .def("surface_contact_frame_names", &Robot::surfaceContactFrameNames)
    .def("generalized_momentum_bias", &Robot::generalizedMomentumBias)
    .def("set_generalized_momentum_bias", &Robot::setGeneralizedMomentumBias,
          py::arg("generalized_momentum_bias"))
    .def("joint_effort_limit", &Robot::jointEffortLimit)
    .def("joint_velocity_limit", &Robot::jointVelocityLimit)
    .def("lower_joint_position_limit", &Robot::lowerJointPositionLimit)
    .def("upper_joint_position_limit", &Robot::upperJointPositionLimit)
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
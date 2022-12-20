#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/constraints/pdipm.hpp"
#include "robotoc/constraints/friction_cone.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class FrictionConeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier_param = 1.0e-03;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    mu = 0.7;
    fraction_to_boundary_rule = 0.995;
    cone.resize(5, 3);
    cone.setZero();
    cone <<  0,  0, -1, 
             1,  0, -(mu/std::sqrt(2)),
            -1,  0, -(mu/std::sqrt(2)),
             0,  1, -(mu/std::sqrt(2)),
             0, -1, -(mu/std::sqrt(2));
    contact_surface = Eigen::Quaterniond::UnitRandom().toRotationMatrix();
    cone_surface_local = cone * contact_surface.transpose();
    grid_info = GridInfo::Random();
  }

  virtual void TearDown() {
  }

  void test_kinematics(Robot& robot, const ContactStatus& contact_status) const;
  void test_isFeasible(Robot& robot, const ContactStatus& contact_status) const;
  void test_setSlack(Robot& robot, const ContactStatus& contact_status) const;
  void test_evalConstraint(Robot& robot, const ContactStatus& contact_status) const;
  void test_evalDerivatives(Robot& robot, const ContactStatus& contact_status) const;
  void test_condenseSlackAndDual(Robot& robot, 
                                const ContactStatus& contact_status) const;
  void test_expandSlackAndDual(Robot& robot, const ContactStatus& contact_status) const;

  double barrier_param, dt, mu, fraction_to_boundary_rule;
  Eigen::MatrixXd cone, cone_surface_local;
  Eigen::Matrix3d contact_surface;
  GridInfo grid_info;
};


void FrictionConeTest::test_kinematics(Robot& robot, 
                                      const ContactStatus& contact_status) const {
  FrictionCone constr(robot); 
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void FrictionConeTest::test_isFeasible(Robot& robot, 
                                      const ContactStatus& contact_status) const {
  FrictionCone constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  EXPECT_EQ(constr.dimc(), 5*contact_status.maxNumContacts());
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateFrameKinematics(s.q);
  if (contact_status.hasActiveContacts()) {
    bool feasible = true;
    for (int i=0; i<contact_status.maxNumContacts(); ++i) {
      if (contact_status.isContactActive(i)) {
        Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
        robot.transformFromLocalToWorld(robot.contactFrames()[i], s.f[i].template head<3>(), f_world);
        Eigen::VectorXd res =  Eigen::VectorXd::Zero(5);
        FrictionCone::frictionConeResidual(mu, f_world, contact_surface, res);
        if (res.maxCoeff() > 0) {
          feasible = false;
        }
      }
    }
    EXPECT_EQ(constr.isFeasible(robot, contact_status, grid_info, s, data), feasible);
  }
}


void FrictionConeTest::test_setSlack(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam()), data_ref(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  constr.allocateExtraData(data_ref);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateFrameKinematics(s.q);
  constr.setSlack(robot, contact_status, grid_info, s, data);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
    robot.transformFromLocalToWorld(robot.contactFrames()[i], s.f[i].template head<3>(), f_world);
    FrictionCone::frictionConeResidual(mu, f_world, contact_surface, data_ref.residual.segment(5*i, 5));
    data_ref.slack.segment(5*i, 5) = - data_ref.residual.segment(5*i, 5);
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void FrictionConeTest::test_evalConstraint(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone constr(robot); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  constr.evalConstraint(robot, contact_status, grid_info, s, data);
  data_ref.residual.setZero();
  data_ref.cmpl.setZero();
  data_ref.log_barrier = 0;
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      robot.transformFromLocalToWorld(robot.contactFrames()[i], s.f[i].template head<3>(), f_world);
      FrictionCone::frictionConeResidual(mu, f_world, contact_surface, data_ref.residual.segment(5*i, 5));
      data_ref.residual.template segment<5>(5*i) += data_ref.slack.segment(5*i, 5);
      for (int j=0; j<5; ++j) {
        data_ref.cmpl.coeffRef(5*i+j) 
            = data_ref.slack.coeff(5*i+j) * data_ref.dual.coeff(5*i+j) - barrier_param;
      }
      data_ref.log_barrier += pdipm::logBarrier(barrier_param, data_ref.slack.segment(5*i, 5));
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void FrictionConeTest::test_evalDerivatives(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, contact_status, grid_info, s, data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  constr.evalConstraint(robot, contact_status, grid_info, s, data);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      robot.transformFromLocalToWorld(robot.contactFrames()[i], s.f[i].template head<3>(), f_world);
      Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, robot.dimv());
      robot.getFrameJacobian(robot.contactFrames()[i], J);
      Eigen::MatrixXd dfW_dq = Eigen::MatrixXd::Zero(3, robot.dimv());
      for (int j=0; j<robot.dimv(); ++j) {
        dfW_dq.col(j) = J.template bottomRows<3>().col(j).cross(f_world);
      }
      const Eigen::MatrixXd dg_dq = cone_surface_local * dfW_dq;
      const Eigen::MatrixXd dg_df = cone_surface_local * robot.frameRotation(robot.contactFrames()[i]);
      kkt_res_ref.lq().noalias() 
          += dg_dq.transpose() * data_ref.dual.segment(5*i, 5);
      kkt_res_ref.lf().segment(dimf_stack, 3)
          += dg_df.transpose() * data_ref.dual.segment(5*i, 5);
      switch (robot.contactType(i)) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void FrictionConeTest::test_condenseSlackAndDual(Robot& robot, 
                                                const ContactStatus& contact_status) const {
  FrictionCone constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, contact_status, grid_info, s, data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto kkt_mat = SplitKKTMatrix::Random(robot, contact_status);
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  constr.evalConstraint(robot, contact_status, grid_info, s, data);
  constr.evalDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  auto data_ref = data;
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(contact_status, grid_info, data, kkt_mat, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      robot.transformFromLocalToWorld(robot.contactFrames()[i], s.f[i].template head<3>(), f_world);
      Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, robot.dimv());
      robot.getFrameJacobian(robot.contactFrames()[i], J);
      Eigen::MatrixXd dfW_dq = Eigen::MatrixXd::Zero(3, robot.dimv());
      for (int j=0; j<robot.dimv(); ++j) {
        dfW_dq.col(j) = J.template bottomRows<3>().col(j).cross(f_world);
      }
      const Eigen::MatrixXd dg_dq = cone_surface_local * dfW_dq;
      const Eigen::MatrixXd dg_df = cone_surface_local * robot.frameRotation(robot.contactFrames()[i]);
      Eigen::VectorXd r(5);
      r.array() = (data_ref.dual.segment(5*i, 5).array()*data_ref.residual.segment(5*i, 5).array()-data_ref.cmpl.segment(5*i, 5).array()) 
                  / data_ref.slack.segment(5*i, 5).array();
      kkt_res_ref.lq() += dg_dq.transpose() * r;
      kkt_res_ref.lf().template segment<3>(dimf_stack) += dg_df.transpose() * r;
      r.array() = data_ref.dual.segment(5*i, 5).array() 
                   / data_ref.slack.segment(5*i, 5).array();
      kkt_mat_ref.Qqq()
          += dg_dq.transpose() * r.asDiagonal() * dg_dq;
      kkt_mat_ref.Qqf().middleCols(dimf_stack, 3)
          += dg_dq.transpose() * r.asDiagonal() * dg_df;
      kkt_mat_ref.Qff().block(dimf_stack, dimf_stack, 3, 3)
          += dg_df.transpose() * r.asDiagonal() * dg_df;
      switch (robot.contactType(i)) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void FrictionConeTest::test_expandSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  constr.setSlack(robot, contact_status, grid_info, s, data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  data.dslack.setRandom();
  data.ddual.setRandom();
  auto kkt_mat = SplitKKTMatrix::Random(robot, contact_status);
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  constr.evalConstraint(robot, contact_status, grid_info, s, data);
  constr.evalDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  constr.condenseSlackAndDual(contact_status, grid_info, data, kkt_mat, kkt_res);
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot, contact_status);
  constr.expandSlackAndDual(contact_status, grid_info, d, data);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      robot.transformFromLocalToWorld(robot.contactFrames()[i], s.f[i].template head<3>(), f_world);
      Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, robot.dimv());
      robot.getFrameJacobian(robot.contactFrames()[i], J);
      Eigen::MatrixXd dfW_dq = Eigen::MatrixXd::Zero(3, robot.dimv());
      for (int j=0; j<robot.dimv(); ++j) {
        dfW_dq.col(j) = J.template bottomRows<3>().col(j).cross(f_world);
      }
      const Eigen::MatrixXd dg_dq = cone_surface_local * dfW_dq;
      const Eigen::MatrixXd dg_df = cone_surface_local * robot.frameRotation(robot.contactFrames()[i]);
      data_ref.dslack.segment(5*i, 5)
          = - dg_dq * d.dq() - dg_df * d.df().segment(dimf_stack, 3) 
            - data_ref.residual.segment(5*i, 5);
      for (int j=0; j<5; ++j) {
        data_ref.ddual(5*i+j) 
          = - (data_ref.dual(5*i+j)*data_ref.dslack(5*i+j)+data_ref.cmpl(5*i+j))
              / data_ref.slack(5*i+j);
      }
      switch (robot.contactType(i)) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(FrictionConeTest, frictionConeResidual) {
  const Eigen::Vector3d f = Eigen::Vector3d::Random();
  const Eigen::VectorXd res_ref = cone_surface_local * f;
  Eigen::VectorXd res = Eigen::VectorXd::Zero(5);
  FrictionCone::frictionConeResidual(mu, f, contact_surface, res);
  EXPECT_TRUE(res.isApprox(res_ref));
}


TEST_F(FrictionConeTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random(), contact_surface);
  }
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
  contact_status.activateContact(0);
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
}


TEST_F(FrictionConeTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random(), contact_surface);
  }
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
  contact_status.setRandom();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random(), contact_surface);
  }
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
}


TEST_F(FrictionConeTest, humanoidRobot) {
  auto robot = testhelper::CreateHumanoidRobot(dt);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random(), contact_surface);
  }
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
  contact_status.setRandom();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random(), contact_surface);
  }
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
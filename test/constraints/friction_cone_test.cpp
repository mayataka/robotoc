#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/pdipm.hpp"
#include "idocp/constraints/friction_cone.hpp"

#include "robot_factory.hpp"

namespace idocp {

class FrictionConeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier = 1.0e-04;
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
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot, const ContactStatus& contact_status) const;
  void testfLocal2World(Robot& robot, const ContactStatus& contact_status) const;
  void testIsFeasible(Robot& robot, const ContactStatus& contact_status) const;
  void testSetSlack(Robot& robot, const ContactStatus& contact_status) const;
  void testComputePrimalAndDualResidual(Robot& robot, 
                                        const ContactStatus& contact_status) const;
  void testComputePrimalResidualDerivatives(Robot& robot, 
                                            const ContactStatus& contact_status) const;
  void testCondenseSlackAndDual(Robot& robot, 
                                const ContactStatus& contact_status) const;
  void testExpandSlackAndDual(Robot& robot, const ContactStatus& contact_status) const;

  double barrier, dt, mu, fraction_to_boundary_rule;
  Eigen::MatrixXd cone;
};


void FrictionConeTest::testKinematics(Robot& robot, 
                                      const ContactStatus& contact_status) const {
  FrictionCone constr(robot, mu); 
  EXPECT_TRUE(constr.useKinematics());
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void FrictionConeTest::testfLocal2World(Robot& robot, 
                                        const ContactStatus& contact_status) const {
  const Eigen::Vector3d f_local = Eigen::Vector3d::Random();
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  robot.updateFrameKinematics(q);
  for (const auto frame : robot.contactFrames()) {
    const Eigen::Vector3d f_world_ref = robot.frameRotation(frame) * f_local;
    Eigen::Vector3d f_world = Eigen::Vector3d::Random();
    FrictionCone::fLocal2World(robot, frame, f_local, f_world);
    EXPECT_TRUE(f_world.isApprox(f_world_ref));
  }
}


void FrictionConeTest::testIsFeasible(Robot& robot, 
                                      const ContactStatus& contact_status) const {
  FrictionCone constr(robot, mu); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  constr.allocateExtraData(data);
  EXPECT_EQ(constr.dimc(), 5*contact_status.maxPointContacts());
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateFrameKinematics(s.q);
  if (contact_status.hasActiveContacts()) {
    bool feasible = true;
    for (int i=0; i<contact_status.maxPointContacts(); ++i) {
      if (contact_status.isContactActive(i)) {
        Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
        FrictionCone::fLocal2World(robot, robot.contactFrames()[i], s.f[i], f_world);
        Eigen::VectorXd res =  Eigen::VectorXd::Zero(5);
        FrictionCone::frictionConeResidual(mu, f_world, res);
        if (res.maxCoeff() > 0) {
          feasible = false;
        }
      }
    }
    EXPECT_EQ(constr.isFeasible(robot, data, s), feasible);
  }
}


void FrictionConeTest::testSetSlack(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone constr(robot, mu); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter()), data_ref(constr.dimc(), constr.barrierParameter());
  constr.allocateExtraData(data);
  constr.allocateExtraData(data_ref);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateFrameKinematics(s.q);
  constr.setSlack(robot, data, s);
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
    FrictionCone::fLocal2World(robot, robot.contactFrames()[i], s.f[i], f_world);
    FrictionCone::frictionConeResidual(mu, f_world, data_ref.residual.segment(5*i, 5));
    data_ref.slack.segment(5*i, 5) = - data_ref.residual.segment(5*i, 5);
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void FrictionConeTest::testComputePrimalAndDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone constr(robot, mu); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  constr.allocateExtraData(data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  constr.computePrimalAndDualResidual(robot, data, s);
  data_ref.residual.setZero();
  data_ref.cmpl.setZero();
  data_ref.log_barrier = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      FrictionCone::fLocal2World(robot, robot.contactFrames()[i], s.f[i], f_world);
      FrictionCone::frictionConeResidual(mu, f_world, data_ref.residual.segment(5*i, 5));
      data_ref.residual.template segment<5>(5*i) += data_ref.slack.segment(5*i, 5);
      for (int j=0; j<5; ++j) {
        data_ref.cmpl.coeffRef(5*i+j) 
            = data_ref.slack.coeff(5*i+j) * data_ref.dual.coeff(5*i+j) - barrier;
      }
      data_ref.log_barrier += pdipm::LogBarrier(barrier, data_ref.slack.segment(5*i, 5));
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void FrictionConeTest::testComputePrimalResidualDerivatives(Robot& robot, 
                                                            const ContactStatus& contact_status) const {
  FrictionCone constr(robot, mu); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  constr.computePrimalAndDualResidual(robot, data, s);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  auto kkt_res_ref = kkt_res;
  constr.computePrimalResidualDerivatives(robot, data, dt, s, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      FrictionCone::fLocal2World(robot, robot.contactFrames()[i], s.f[i], f_world);
      Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, robot.dimv());
      robot.getFrameJacobian(robot.contactFrames()[i], J);
      Eigen::MatrixXd dfW_dq = Eigen::MatrixXd::Zero(3, robot.dimv());
      for (int j=0; j<robot.dimv(); ++j) {
        dfW_dq.col(j) = J.template bottomRows<3>().col(j).cross(f_world);
      }
      const Eigen::MatrixXd dg_dq = cone * dfW_dq;
      const Eigen::MatrixXd dg_df = cone * robot.frameRotation(robot.contactFrames()[i]);
      kkt_res_ref.lq().noalias() 
          += dt * dg_dq.transpose() * data_ref.dual.segment(5*i, 5);
      kkt_res_ref.lf().segment(dimf_stack, 3)
          += dt * dg_df.transpose() * data_ref.dual.segment(5*i, 5);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void FrictionConeTest::testCondenseSlackAndDual(Robot& robot, 
                                                const ContactStatus& contact_status) const {
  FrictionCone constr(robot, mu); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto kkt_mat = SplitKKTMatrix::Random(robot, contact_status);
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  constr.computePrimalAndDualResidual(robot, data, s);
  constr.computePrimalResidualDerivatives(robot, data, dt, s, kkt_res);
  auto data_ref = data;
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      FrictionCone::fLocal2World(robot, robot.contactFrames()[i], s.f[i], f_world);
      Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, robot.dimv());
      robot.getFrameJacobian(robot.contactFrames()[i], J);
      Eigen::MatrixXd dfW_dq = Eigen::MatrixXd::Zero(3, robot.dimv());
      for (int j=0; j<robot.dimv(); ++j) {
        dfW_dq.col(j) = J.template bottomRows<3>().col(j).cross(f_world);
      }
      const Eigen::MatrixXd dg_dq = cone * dfW_dq;
      const Eigen::MatrixXd dg_df = cone * robot.frameRotation(robot.contactFrames()[i]);
      Eigen::VectorXd r(5);
      r.array() = (data_ref.dual.segment(5*i, 5).array()*data_ref.residual.segment(5*i, 5).array()-data_ref.cmpl.segment(5*i, 5).array()) 
                  / data_ref.slack.segment(5*i, 5).array();
      kkt_res_ref.lq() += dt * dg_dq.transpose() * r;
      kkt_res_ref.lf().template segment<3>(dimf_stack) += dt * dg_df.transpose() * r;
      r.array() = data_ref.dual.segment(5*i, 5).array() 
                   / data_ref.slack.segment(5*i, 5).array();
      kkt_mat_ref.Qqq()
          += dt * dg_dq.transpose() * r.asDiagonal() * dg_dq;
      kkt_mat_ref.Qqf().middleCols(dimf_stack, 3)
          += dt * dg_dq.transpose() * r.asDiagonal() * dg_df;
      kkt_mat_ref.Qff().block(dimf_stack, dimf_stack, 3, 3)
          += dt * dg_df.transpose() * r.asDiagonal() * dg_df;
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void FrictionConeTest::testExpandSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone constr(robot, mu); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  constr.setSlack(robot, data, s);
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
  constr.computePrimalAndDualResidual(robot, data, s);
  constr.computePrimalResidualDerivatives(robot, data, dt, s, kkt_res);
  constr.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot, contact_status);
  constr.expandSlackAndDual(data, s, d);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::Vector3d f_world = Eigen::Vector3d::Zero();
      FrictionCone::fLocal2World(robot, robot.contactFrames()[i], s.f[i], f_world);
      Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, robot.dimv());
      robot.getFrameJacobian(robot.contactFrames()[i], J);
      Eigen::MatrixXd dfW_dq = Eigen::MatrixXd::Zero(3, robot.dimv());
      for (int j=0; j<robot.dimv(); ++j) {
        dfW_dq.col(j) = J.template bottomRows<3>().col(j).cross(f_world);
      }
      const Eigen::MatrixXd dg_dq = cone * dfW_dq;
      const Eigen::MatrixXd dg_df = cone * robot.frameRotation(robot.contactFrames()[i]);
      data_ref.dslack.segment(5*i, 5)
          = - dg_dq * d.dq() - dg_df * d.df().segment(dimf_stack, 3) 
            - data_ref.residual.segment(5*i, 5);
      for (int j=0; j<5; ++j) {
        data_ref.ddual(5*i+j) 
          = - (data_ref.dual(5*i+j)*data_ref.dslack(5*i+j)+data_ref.cmpl(5*i+j))
              / data_ref.slack(5*i+j);
      }
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(FrictionConeTest, frictionConeResidual) {
  const Eigen::Vector3d f = Eigen::Vector3d::Random();
  Eigen::VectorXd res_ref = Eigen::VectorXd::Zero(5);
  res_ref(0) = - f(2);
  res_ref(1) =   f(0) - mu * f(2) / std::sqrt(2);
  res_ref(2) = - f(0) - mu * f(2) / std::sqrt(2);
  res_ref(3) =   f(1) - mu * f(2) / std::sqrt(2);
  res_ref(4) = - f(1) - mu * f(2) / std::sqrt(2);
  Eigen::VectorXd res = Eigen::VectorXd::Zero(5);
  FrictionCone::frictionConeResidual(mu, f, res);
  EXPECT_TRUE(res.isApprox(res_ref));
}


TEST_F(FrictionConeTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  testKinematics(robot, contact_status);
  testfLocal2World(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlack(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testComputePrimalResidualDerivatives(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testExpandSlackAndDual(robot, contact_status);
  contact_status.activateContact(0);
  testKinematics(robot, contact_status);
  testfLocal2World(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlack(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testComputePrimalResidualDerivatives(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testExpandSlackAndDual(robot, contact_status);
}


TEST_F(FrictionConeTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  testKinematics(robot, contact_status);
  testfLocal2World(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlack(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testComputePrimalResidualDerivatives(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testExpandSlackAndDual(robot, contact_status);
  contact_status.setRandom();
  testKinematics(robot, contact_status);
  testfLocal2World(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlack(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testComputePrimalResidualDerivatives(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testExpandSlackAndDual(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

#include "joint_space_constraints.hpp"

namespace idocp {

class FixedBaseConstraintsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    std::vector<int> contact_frames = {18};
    const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
    const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
    robot = Robot(urdf, contact_frames, baum_a, baum_b);
    std::random_device rnd;
    contact_status.push_back(rnd()%2==0);
    robot.setContactStatus(contact_status);
    s = SplitSolution(robot);
    s.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s.q);
    s.v = Eigen::VectorXd::Random(robot.dimv());
    s.a = Eigen::VectorXd::Random(robot.dimv());
    s.u = Eigen::VectorXd::Random(robot.dimv());
    s.f = Eigen::VectorXd::Random(robot.max_dimf());
    s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s.lmd = Eigen::VectorXd::Random(robot.dimv());
    s.gmm = Eigen::VectorXd::Random(robot.dimv());
    d = SplitDirection(robot);
    d.dq() = Eigen::VectorXd::Random(robot.dimv());
    d.dv() = Eigen::VectorXd::Random(robot.dimv());
    d.da() = Eigen::VectorXd::Random(robot.dimv());
    d.df() = Eigen::VectorXd::Random(robot.dimf());
    d.du = Eigen::VectorXd::Random(robot.dimv());
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    constraints = std::make_shared<Constraints>();
    auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
    auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
    auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
    auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
    auto torques_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
    auto torques_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
    constraints->push_back(joint_upper_limit); 
    constraints->push_back(joint_lower_limit);
    constraints->push_back(velocity_lower_limit); 
    constraints->push_back(velocity_upper_limit);
    constraints->push_back(torques_lower_limit); 
    constraints->push_back(torques_upper_limit);
    data = constraints->createConstraintsData(robot);
    constraints_ref = pdipmtest::JointSpaceConstraints(robot);
    kkt_matrix = KKTMatrix(robot);
    kkt_residual = KKTResidual(robot);
    lq_ref = Eigen::VectorXd::Zero(robot.dimv());
    lv_ref = Eigen::VectorXd::Zero(robot.dimv());
    la_ref = Eigen::VectorXd::Zero(robot.dimv());
    lu_ref = Eigen::VectorXd::Zero(robot.dimv());
    Qqq_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
    Qvv_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
    Qaa_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
    Quu_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  std::vector<bool> contact_status;
  std::shared_ptr<Constraints> constraints;
  ConstraintsData data;
  SplitSolution s;
  SplitDirection d;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  pdipmtest::JointSpaceConstraints constraints_ref;
  Eigen::VectorXd lq_ref, lv_ref, la_ref, lu_ref;
  Eigen::MatrixXd Qqq_ref, Qvv_ref, Qaa_ref, Quu_ref;
};


TEST_F(FixedBaseConstraintsTest, useKinematics) {
  EXPECT_FALSE(constraints->useKinematics());
}


TEST_F(FixedBaseConstraintsTest, isFeasible) {
  ASSERT_EQ(data.size(), 6);
  for (int i=0; i<6; ++i) {
    ASSERT_EQ(data[i].slack.size(), robot.dimv());
    ASSERT_EQ(data[i].dual.size(), robot.dimv());
    ASSERT_EQ(data[i].residual.size(), robot.dimv());
    ASSERT_EQ(data[i].duality.size(), robot.dimv());
    ASSERT_EQ(data[i].dslack.size(), robot.dimv());
    ASSERT_EQ(data[i].ddual.size(), robot.dimv());
    ASSERT_TRUE(data[i].slack.isZero());
    ASSERT_TRUE(data[i].dual.isZero());
    ASSERT_TRUE(data[i].residual.isZero());
    ASSERT_TRUE(data[i].duality.isZero());
    ASSERT_TRUE(data[i].dslack.isZero());
    ASSERT_TRUE(data[i].ddual.isZero());
  }
  EXPECT_EQ(constraints->isFeasible(robot, data, s),
            constraints_ref.isFeasible(s.q, s.v, s.a, s.u));
}


TEST_F(FixedBaseConstraintsTest, augmentDualResidual) {
  constraints->setSlackAndDual(robot, data, dtau, s);
  constraints_ref.setSlackAndDual(dtau, s.q, s.v, s.a, s.u);
  constraints->augmentDualResidual(robot, data, dtau, kkt_residual.lu);
  constraints_ref.augmentDualResidual(dtau, lu_ref);
  EXPECT_TRUE(kkt_residual.lu.isApprox(lu_ref));
  constraints->augmentDualResidual(robot, data, dtau, kkt_residual);
  constraints_ref.augmentDualResidual(dtau, lq_ref, lv_ref, la_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isApprox(lv_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lu.isApprox(lu_ref));
  EXPECT_DOUBLE_EQ(constraints->residualL1Nrom(robot, data, dtau, s),
                   constraints_ref.residualL1Nrom(dtau, s.q, s.v, s.a, s.u));
  EXPECT_DOUBLE_EQ(constraints->squaredKKTErrorNorm(robot, data, dtau, s),
                   constraints_ref.residualSquaredNrom(dtau, s.q, s.v, s.a, s.u));
}


TEST_F(FixedBaseConstraintsTest, condenseSlackAndDual) {
  constraints->setSlackAndDual(robot, data, dtau, s);
  constraints_ref.setSlackAndDual(dtau, s.q, s.v, s.a, s.u);
  constraints->condenseSlackAndDual(robot, data, dtau, s.u, kkt_matrix.Quu, 
                                    kkt_residual.lu);
  constraints_ref.condenseSlackAndDual(dtau, s.u, Quu_ref, lu_ref);
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(Quu_ref));
  EXPECT_TRUE(kkt_residual.lu.isApprox(lu_ref));
  constraints->condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, 
                                    kkt_residual);
  constraints_ref.condenseSlackAndDual(dtau, s.q, s.v, s.a, Qqq_ref, Qvv_ref, 
                                       Qaa_ref, lq_ref, lv_ref, la_ref);
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(Qqq_ref));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(Qvv_ref));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(Qaa_ref));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(Quu_ref));
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isApprox(lv_ref));
  EXPECT_TRUE(kkt_residual.la().isApprox(la_ref));
  EXPECT_TRUE(kkt_residual.lu.isApprox(lu_ref));
}


TEST_F(FixedBaseConstraintsTest, updateSlackAndDualDirection) {
  constraints->setSlackAndDual(robot, data, dtau, s);
  constraints_ref.setSlackAndDual(dtau, s.q, s.v, s.a, s.u);
  constraints->computeSlackAndDualDirection(robot, data, dtau, d);
  constraints_ref.computeSlackAndDualDirection(dtau, d.dq(), d.dv(), d.da(), d.du);
  EXPECT_DOUBLE_EQ(constraints->maxSlackStepSize(data), 
                   constraints_ref.maxSlackStepSize());
  EXPECT_DOUBLE_EQ(constraints->maxDualStepSize(data), 
                   constraints_ref.maxDualStepSize());
  const double slack_step_size = constraints->maxSlackStepSize(data);
  const double dual_step_size = constraints->maxDualStepSize(data);
  EXPECT_DOUBLE_EQ(constraints->costSlackBarrier(data), 
                   constraints_ref.costSlackBarrier());
  EXPECT_DOUBLE_EQ(constraints->costSlackBarrier(data, slack_step_size), 
                   constraints_ref.costSlackBarrier(slack_step_size));
  constraints->updateSlack(data, slack_step_size);
  constraints_ref.updateSlack(slack_step_size);
  constraints->updateDual(data, dual_step_size);
  constraints_ref.updateDual(dual_step_size);
  constraints->augmentDualResidual(robot, data, dtau, kkt_residual.lu);
  constraints_ref.augmentDualResidual(dtau, lu_ref);
  EXPECT_TRUE(kkt_residual.lu.isApprox(lu_ref));
  constraints->augmentDualResidual(robot, data, dtau, kkt_residual);
  constraints_ref.augmentDualResidual(dtau, lq_ref, lv_ref, la_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isApprox(lv_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lu.isApprox(lu_ref));
  EXPECT_DOUBLE_EQ(constraints->residualL1Nrom(robot, data, dtau, s),
                   constraints_ref.residualL1Nrom(dtau, s.q, s.v, s.a, s.u));
  EXPECT_DOUBLE_EQ(constraints->squaredKKTErrorNorm(robot, data, dtau, s),
                   constraints_ref.residualSquaredNrom(dtau, s.q, s.v, s.a, s.u));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
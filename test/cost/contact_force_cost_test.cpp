#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"

#include "idocp/utils/derivative_checker.hpp"

#include "robot_factory.hpp"

namespace idocp {

class ContactForceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot) const;
  void testTerminalCost(Robot& robot) const;
  void testImpulseCost(Robot& robot) const;

  double dt, t;
};


void ContactForceCostTest::testStageCost(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref, fi_weight, fi_ref;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
    fi_weight.push_back(Eigen::Vector3d::Random());
    fi_ref.push_back(Eigen::Vector3d::Random());
  }
  auto cost = std::make_shared<ContactForceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_f_weight(f_weight);
  cost->set_f_ref(f_ref);
  // cost->set_f_ref(robot);
  cost->set_fi_weight(fi_weight);
  cost->set_fi_ref(fi_ref);
  ContactStatus contact_status = robot.createContactStatus();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    contact_status.deactivateContact(i);
  }
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_mat.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  EXPECT_DOUBLE_EQ(cost->computeStageCost(robot, data, t, dt, s), 0);
  cost->computeStageCostDerivatives(robot, data, t, dt, s, kkt_res);
  cost->computeStageCostHessian(robot, data, t, dt, s, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  contact_status.setRandom();
  s.setRandom(robot, contact_status);
  double l_ref = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      l_ref += (f_weight[i].array() * (s.f[i].array()-f_ref[i].array()) 
                                    * (s.f[i].array()-f_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost->computeStageCost(robot, data, t, dt, s), 0.5*dt*l_ref);
  kkt_res.setContactStatus(contact_status);
  kkt_mat.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  cost->computeStageCostDerivatives(robot, data, t, dt, s, kkt_res);
  cost->computeStageCostHessian(robot, data, t, dt, s, kkt_mat);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_res_ref.lf().segment<3>(dimf_stack).array()
          += dt * f_weight[i].array() * (s.f[i].array()-f_ref[i].array());
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_mat_ref.Qff().diagonal().segment<3>(dimf_stack) += dt * f_weight[i];
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost, contact_status));
  EXPECT_TRUE(derivative_checker.checkSecondOrderStageCostDerivatives(cost, contact_status));
}


void ContactForceCostTest::testTerminalCost(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref, fi_weight, fi_ref;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
    fi_weight.push_back(Eigen::Vector3d::Random());
    fi_ref.push_back(Eigen::Vector3d::Random());
  }
  auto cost = std::make_shared<ContactForceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_f_weight(f_weight);
  cost->set_f_ref(f_ref);
  cost->set_fi_weight(fi_weight);
  cost->set_fi_ref(fi_ref);
  ContactStatus contact_status = robot.createContactStatus();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    contact_status.deactivateContact(i);
  }
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_mat.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  EXPECT_DOUBLE_EQ(cost->computeTerminalCost(robot, data, t, s), 0);
  cost->computeTerminalCostDerivatives(robot, data, t, s, kkt_res);
  cost->computeTerminalCostHessian(robot, data, t, s, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  contact_status.setRandom();
  s.setRandom(robot, contact_status);
  EXPECT_DOUBLE_EQ(cost->computeTerminalCost(robot, data, t, s), 0);
  kkt_res.setContactStatus(contact_status);
  kkt_mat.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  cost->computeTerminalCostDerivatives(robot, data, t, s, kkt_res);
  cost->computeTerminalCostHessian(robot, data, t, s, kkt_mat);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
  EXPECT_TRUE(derivative_checker.checkSecondOrderTerminalCostDerivatives(cost));
}


void ContactForceCostTest::testImpulseCost(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref, fi_weight, fi_ref;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
    fi_weight.push_back(Eigen::Vector3d::Random());
    fi_ref.push_back(Eigen::Vector3d::Random());
  }
  auto cost = std::make_shared<ContactForceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_f_weight(f_weight);
  cost->set_f_ref(f_ref);
  cost->set_fi_weight(fi_weight);
  cost->set_fi_ref(fi_ref);
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    impulse_status.deactivateImpulse(i);
  }
  ImpulseSplitKKTMatrix kkt_mat(robot);
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_mat.setImpulseStatus(impulse_status);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  EXPECT_DOUBLE_EQ(cost->computeImpulseCost(robot, data, t, s), 0);
  cost->computeImpulseCostDerivatives(robot, data, t, s, kkt_res);
  cost->computeImpulseCostHessian(robot, data, t, s, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  impulse_status.setRandom();
  s.setRandom(robot, impulse_status);
  double l_ref = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      l_ref += (fi_weight[i].array() * (s.f[i].array()-fi_ref[i].array()) 
                                     * (s.f[i].array()-fi_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost->computeImpulseCost(robot, data, t, s), 0.5*l_ref);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_mat.setImpulseStatus(impulse_status);
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  ImpulseSplitKKTResidual kkt_res_ref = kkt_res;
  ImpulseSplitKKTMatrix kkt_mat_ref = kkt_mat;
  cost->computeImpulseCostDerivatives(robot, data, t, s, kkt_res);
  cost->computeImpulseCostHessian(robot, data, t, s, kkt_mat);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_res_ref.lf().segment<3>(dimf_stack).array()
          += fi_weight[i].array() * (s.f[i].array()-fi_ref[i].array());
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_mat_ref.Qff().diagonal().segment<3>(dimf_stack) += fi_weight[i];
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost, impulse_status));
  EXPECT_TRUE(derivative_checker.checkSecondOrderImpulseCostDerivatives(cost, impulse_status));
}


TEST_F(ContactForceCostTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}


TEST_F(ContactForceCostTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
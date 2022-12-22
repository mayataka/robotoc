#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/cost/local_contact_force_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class LocalContactForceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    t = grid_info.t;
    dt = grid_info.dt;
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot) const;
  void testTerminalCost(Robot& robot) const;
  void testImpactCost(Robot& robot) const;

  GridInfo grid_info;
  double dt, t;
};


void LocalContactForceCostTest::testStageCost(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref, fi_weight, fi_ref;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
    fi_weight.push_back(Eigen::Vector3d::Random());
    fi_ref.push_back(Eigen::Vector3d::Random());
  }
  auto cost = std::make_shared<LocalContactForceCost>(robot);
  CostFunctionData data(robot);
  cost->set_f_weight(f_weight);
  cost->set_f_ref(f_ref);
  cost->set_fi_weight(fi_weight);
  cost->set_fi_ref(fi_ref);
  ContactStatus contact_status = robot.createContactStatus();
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    contact_status.deactivateContact(i);
  }
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactDimension(contact_status.dimf());
  kkt_mat.setContactDimension(contact_status.dimf());
  auto s = SplitSolution::Random(robot, contact_status);
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, grid_info, s, data), 0);
  cost->evalStageCostDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, grid_info, s, data, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  contact_status.setRandom();
  s.setRandom(robot, contact_status);
  double l_ref = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      l_ref += (f_weight[i].array() * (fl.array()-f_ref[i].array()) 
                                    * (fl.array()-f_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, grid_info, s, data), 0.5*l_ref);
  kkt_res.setContactDimension(contact_status.dimf());
  kkt_mat.setContactDimension(contact_status.dimf());
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  auto kkt_res_ref = kkt_res;
  auto kkt_mat_ref = kkt_mat;
  cost->evalStageCostDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, grid_info, s, data, kkt_mat);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      kkt_res_ref.lf().segment<3>(dimf_stack).array()
          += f_weight[i].array() * (fl.array()-f_ref[i].array());
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
  dimf_stack = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_mat_ref.Qff().diagonal().segment<3>(dimf_stack) += f_weight[i];
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
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost, contact_status));
  EXPECT_TRUE(derivative_checker.checkSecondOrderStageCostDerivatives(cost, contact_status));
}


void LocalContactForceCostTest::testTerminalCost(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref, fi_weight, fi_ref;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
    fi_weight.push_back(Eigen::Vector3d::Random());
    fi_ref.push_back(Eigen::Vector3d::Random());
  }
  auto cost = std::make_shared<LocalContactForceCost>(robot);
  CostFunctionData data(robot);
  cost->set_f_weight(f_weight);
  cost->set_f_ref(f_ref);
  cost->set_fi_weight(fi_weight);
  cost->set_fi_ref(fi_ref);
  ContactStatus contact_status = robot.createContactStatus();
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    contact_status.deactivateContact(i);
  }
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactDimension(contact_status.dimf());
  kkt_mat.setContactDimension(contact_status.dimf());
  auto s = SplitSolution::Random(robot, contact_status);
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, grid_info, s, data), 0);
  cost->evalTerminalCostDerivatives(robot, grid_info, s, data, kkt_res);
  cost->evalTerminalCostHessian(robot, grid_info, s, data, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  contact_status.setRandom();
  s.setRandom(robot, contact_status);
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, grid_info, s, data), 0);
  kkt_res.setContactDimension(contact_status.dimf());
  kkt_mat.setContactDimension(contact_status.dimf());
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  auto kkt_res_ref = kkt_res;
  auto kkt_mat_ref = kkt_mat;
  cost->evalTerminalCostDerivatives(robot, grid_info, s, data, kkt_res);
  cost->evalTerminalCostHessian(robot, grid_info, s, data, kkt_mat);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
  EXPECT_TRUE(derivative_checker.checkSecondOrderTerminalCostDerivatives(cost));
}


void LocalContactForceCostTest::testImpactCost(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref, fi_weight, fi_ref;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
    fi_weight.push_back(Eigen::Vector3d::Random());
    fi_ref.push_back(Eigen::Vector3d::Random());
  }
  auto cost = std::make_shared<LocalContactForceCost>(robot);
  CostFunctionData data(robot);
  cost->set_f_weight(f_weight);
  cost->set_f_ref(f_ref);
  cost->set_fi_weight(fi_weight);
  cost->set_fi_ref(fi_ref);
  ImpactStatus impact_status = robot.createImpactStatus();
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    impact_status.deactivateImpact(i);
  }
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactDimension(impact_status.dimf());
  kkt_mat.setContactDimension(impact_status.dimf());
  SplitSolution s = SplitSolution::Random(robot, impact_status);
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, grid_info, s, data), 0);
  cost->evalImpactCostDerivatives(robot, impact_status, grid_info, s, data, kkt_res);
  cost->evalImpactCostHessian(robot, impact_status, grid_info, s, data, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  impact_status.setRandom();
  s.setRandom(robot, impact_status);
  double l_ref = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (impact_status.isImpactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      l_ref += (fi_weight[i].array() * (fl.array()-fi_ref[i].array()) 
                                     * (fl.array()-fi_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, grid_info, s, data), 0.5*l_ref);
  kkt_res.setContactDimension(impact_status.dimf());
  kkt_mat.setContactDimension(impact_status.dimf());
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  auto kkt_res_ref = kkt_res;
  auto kkt_mat_ref = kkt_mat;
  cost->evalImpactCostDerivatives(robot, impact_status, grid_info, s, data, kkt_res);
  cost->evalImpactCostHessian(robot, impact_status, grid_info, s, data, kkt_mat);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (impact_status.isImpactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      kkt_res_ref.lf().segment<3>(dimf_stack).array()
          += fi_weight[i].array() * (fl.array()-fi_ref[i].array());
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
  dimf_stack = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (impact_status.isImpactActive(i)) {
      kkt_mat_ref.Qff().diagonal().segment<3>(dimf_stack) += fi_weight[i];
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
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpactCostDerivatives(cost, impact_status));
  EXPECT_TRUE(derivative_checker.checkSecondOrderImpactCostDerivatives(cost, impact_status));
}


TEST_F(LocalContactForceCostTest, defaultConstructor) {
  EXPECT_NO_THROW(
    auto cost = std::make_shared<LocalContactForceCost>();
  );
}


TEST_F(LocalContactForceCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpactCost(robot);
}


TEST_F(LocalContactForceCostTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpactCost(robot);
}


TEST_F(LocalContactForceCostTest, humanoidRobot) {
  auto robot = testhelper::CreateHumanoidRobot(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpactCost(robot);
}
} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
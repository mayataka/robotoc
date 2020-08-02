#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "robot/robot.hpp"
#include "ocp/split_terminal_ocp.hpp"
#include "manipulator/cost_function.hpp"
#include "manipulator/constraints.hpp"
#include "quadruped/cost_function.hpp"
#include "quadruped/constraints.hpp"


namespace idocp {

class SplitTerminalOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    manipulator_urdf_ = "../../urdf/iiwa14/iiwa14.urdf";
    quadruped_urdf_ = "../../urdf/anymal/anymal.urdf";
    manipulator_ = Robot(manipulator_urdf_);
    quadruped_ = Robot(quadruped_urdf_);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double t_, dtau_;
  std::string manipulator_urdf_, quadruped_urdf_;
  Robot manipulator_, quadruped_;
};


TEST_F(SplitTerminalOCPTest, isFeasible) {
  manipulator::CostFunction cost(manipulator_);
  manipulator::Constraints constraints(manipulator_);
  SplitTerminalOCP ocp(manipulator_, &cost, &constraints);
  const Eigen::VectorXd q = Eigen::VectorXd::Random(manipulator_.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(manipulator_.dimv());
  EXPECT_TRUE(ocp.isFeasible(q, v));
  std::random_device rnd;
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
}


TEST_F(SplitTerminalOCPTest, linearizeOCPFixedBase) {
  manipulator::CostFunction cost(manipulator_);
  manipulator::Constraints constraints(manipulator_);
  SplitTerminalOCP ocp(manipulator_, &cost, &constraints);
  const int dimq = manipulator_.dimq();
  const int dimv = manipulator_.dimv();
  const Eigen::VectorXd q = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  Eigen::MatrixXd Qqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv = Eigen::VectorXd::Zero(dimv);
  std::random_device rnd;
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  ocp.linearizeOCP(manipulator_, t_, lmd, gmm, q, v, Qqq, Qqv, Qvq, Qvv, Qq, Qv);
  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv_ref = Eigen::VectorXd::Zero(dimv);
  cost.phiq(t_, q, v, Qq_ref);
  cost.phiv(t_, q, v, Qv_ref);
  Qq_ref = - Qq_ref + lmd;
  Qv_ref = - Qv_ref + gmm;
  cost.phiqq(t_, q, v, Qqq_ref);
  cost.phivv(t_, q, v, Qvv_ref);
  EXPECT_TRUE(Qqq.isApprox(Qqq_ref));
  EXPECT_TRUE(Qqv.isApprox(Qqv_ref));
  EXPECT_TRUE(Qvq.isApprox(Qvq_ref));
  EXPECT_TRUE(Qvv.isApprox(Qvv_ref));
  EXPECT_TRUE(Qq.isApprox(Qq_ref));
  EXPECT_TRUE(Qv.isApprox(Qv_ref));
}


TEST_F(SplitTerminalOCPTest, linearizeOCPFloatingBase) {
  quadruped::CostFunction cost(quadruped_);
  quadruped::Constraints constraints(quadruped_);
  SplitTerminalOCP ocp(quadruped_, &cost, &constraints);
  const int dimq = quadruped_.dimq();
  const int dimv = quadruped_.dimv();
  pinocchio::Model model;
  pinocchio::urdf::buildModel(quadruped_urdf_, model);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(dimq);
  quadruped_.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  Eigen::MatrixXd Qqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv = Eigen::VectorXd::Zero(dimv);
  std::random_device rnd;
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  ocp.linearizeOCP(quadruped_, t_, lmd, gmm, q, v, Qqq, Qqv, Qvq, Qvv, Qq, Qv);
  cost.setConfigurationJacobian(quadruped_, q);
  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv_ref = Eigen::VectorXd::Zero(dimv);
  cost.phiq(t_, q, v, Qq_ref);
  cost.phiv(t_, q, v, Qv_ref);
  Qq_ref = - Qq_ref + lmd;
  Qv_ref = - Qv_ref + gmm;
  cost.phiqq(t_, q, v, Qqq_ref);
  cost.phivv(t_, q, v, Qvv_ref);
  EXPECT_TRUE(Qqq.isApprox(Qqq_ref));
  EXPECT_TRUE(Qqv.isApprox(Qqv_ref));
  EXPECT_TRUE(Qvq.isApprox(Qvq_ref));
  EXPECT_TRUE(Qvv.isApprox(Qvv_ref));
  EXPECT_TRUE(Qq.isApprox(Qq_ref));
  EXPECT_TRUE(Qv.isApprox(Qv_ref));
}


TEST_F(SplitTerminalOCPTest, terminalCostFixedBase) {
  manipulator::CostFunction cost(manipulator_);
  manipulator::Constraints constraints(manipulator_);
  SplitTerminalOCP ocp(manipulator_, &cost, &constraints);
  const int dimq = manipulator_.dimq();
  const int dimv = manipulator_.dimv();
  std::random_device rnd;
  const Eigen::VectorXd q = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv);
  double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  while (step_size > 1) {
    step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  }
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  const double terminal_cost
      = ocp.terminalCost(t_, q, v);
  const double terminal_cost_with_step
      = ocp.terminalCost(manipulator_, step_size, t_, q, v, dq, dv);
  const double terminal_cost_ref = cost.phi(t_, q, v);
  const double terminal_cost_with_step_ref 
      = cost.phi(t_, q+step_size*dq, v+step_size*dv);
  EXPECT_DOUBLE_EQ(terminal_cost, terminal_cost_ref);
  EXPECT_DOUBLE_EQ(terminal_cost_with_step, terminal_cost_with_step_ref);
}


TEST_F(SplitTerminalOCPTest, terminalCostFloatingBase) {
  manipulator::CostFunction cost(quadruped_);
  manipulator::Constraints constraints(quadruped_);
  SplitTerminalOCP ocp(quadruped_, &cost, &constraints);
  const int dimq = quadruped_.dimq();
  const int dimv = quadruped_.dimv();
  std::random_device rnd;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(dimq);
  quadruped_.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv);
  double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  while (step_size > 1) {
    step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  }
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  const double terminal_cost
      = ocp.terminalCost(t_, q, v);
  const double terminal_cost_with_step
      = ocp.terminalCost(quadruped_, step_size, t_, q, v, dq, dv);
  const double terminal_cost_ref = cost.phi(t_, q, v);
  Eigen::VectorXd q_tmp = q;
  quadruped_.integrateConfiguration(dq, step_size, q_tmp);
  const double terminal_cost_with_step_ref 
      = cost.phi(t_, q_tmp, v+step_size*dv);
  EXPECT_DOUBLE_EQ(terminal_cost, terminal_cost_ref);
  EXPECT_DOUBLE_EQ(terminal_cost_with_step, terminal_cost_with_step_ref);
}


TEST_F(SplitTerminalOCPTest, updateFixedBase) {
  manipulator::CostFunction cost(manipulator_);
  manipulator::Constraints constraints(manipulator_);
  SplitTerminalOCP ocp(manipulator_, &cost, &constraints);
  const int dimq = manipulator_.dimq();
  const int dimv = manipulator_.dimv();
  Eigen::VectorXd q = Eigen::VectorXd::Random(dimq);
  Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd q_ref = q;
  Eigen::VectorXd v_ref = v;
  Eigen::VectorXd lmd_ref = lmd;
  Eigen::VectorXd gmm_ref = gmm;
  std::random_device rnd;
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  const Eigen::MatrixXd Pqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::VectorXd sq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd sv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv);
  double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  while (step_size > 1) {
    step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  }
  ocp.updateDual(step_size);
  ocp.updatePrimal(manipulator_, step_size, Pqq, Pqv, Pvq, Pvv, sq, sv, dq, dv, 
                   lmd, gmm, q, v);
  lmd_ref += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm_ref += step_size * (Pvq * dq + Pvv * dv - sv);
  q_ref += step_size * dq;
  v_ref += step_size * dv;
  EXPECT_TRUE(lmd.isApprox(lmd_ref));
  EXPECT_TRUE(gmm.isApprox(gmm_ref));
  EXPECT_TRUE(q.isApprox(q_ref));
  EXPECT_TRUE(v.isApprox(v_ref));
}


TEST_F(SplitTerminalOCPTest, updateFloatingdBase) {
  manipulator::CostFunction cost(quadruped_);
  manipulator::Constraints constraints(quadruped_);
  SplitTerminalOCP ocp(quadruped_, &cost, &constraints);
  const int dimq = quadruped_.dimq();
  const int dimv = quadruped_.dimv();
  Eigen::VectorXd q = Eigen::VectorXd::Zero(dimq);
  quadruped_.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd q_ref = q;
  Eigen::VectorXd v_ref = v;
  Eigen::VectorXd lmd_ref = lmd;
  Eigen::VectorXd gmm_ref = gmm;
  std::random_device rnd;
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  const Eigen::MatrixXd Pqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::VectorXd sq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd sv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv);
  double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  while (step_size > 1) {
    step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  }
  ocp.updateDual(step_size);
  ocp.updatePrimal(quadruped_, step_size, Pqq, Pqv, Pvq, Pvv, sq, sv, dq, dv, 
                   lmd, gmm, q, v);
  lmd_ref += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm_ref += step_size * (Pvq * dq + Pvv * dv - sv);
  quadruped_.integrateConfiguration(dq, step_size, q_ref);
  v_ref += step_size * dv;
  EXPECT_TRUE(lmd.isApprox(lmd_ref));
  EXPECT_TRUE(gmm.isApprox(gmm_ref));
  EXPECT_TRUE(q.isApprox(q_ref));
  EXPECT_TRUE(v.isApprox(v_ref));
}


TEST_F(SplitTerminalOCPTest, squaredKKTErrorNormFixedBase) {
  manipulator::CostFunction cost(manipulator_);
  manipulator::Constraints constraints(manipulator_);
  SplitTerminalOCP ocp(manipulator_, &cost, &constraints);
  const int dimq = manipulator_.dimq();
  const int dimv = manipulator_.dimv();
  const Eigen::VectorXd q = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  std::random_device rnd;
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  const double result 
      = ocp.squaredKKTErrorNorm(manipulator_, t_, lmd, gmm, q, v);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv);
  cost.phiq(t_, q, v, lq_ref);
  cost.phiv(t_, q, v, lv_ref);
  lq_ref -= lmd;
  lv_ref -= gmm;
  const double result_ref = lq_ref.squaredNorm() + lv_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(result, result_ref);
}


TEST_F(SplitTerminalOCPTest, squaredKKTErrorNormFloatingBase) {
  quadruped::CostFunction cost(quadruped_);
  quadruped::Constraints constraints(quadruped_);
  SplitTerminalOCP ocp(quadruped_, &cost, &constraints);
  const int dimq = quadruped_.dimq();
  const int dimv = quadruped_.dimv();
  Eigen::VectorXd q = Eigen::VectorXd::Zero(dimq);
  quadruped_.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  std::random_device rnd;
  ocp.initConstraints(rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(q, v));
  const double result 
      = ocp.squaredKKTErrorNorm(quadruped_, t_, lmd, gmm, q, v);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv);
  cost.setConfigurationJacobian(quadruped_, q);
  cost.phiq(t_, q, v, lq_ref);
  cost.phiv(t_, q, v, lv_ref);
  lq_ref -= lmd;
  lv_ref -= gmm;
  const double result_ref = lq_ref.squaredNorm() + lv_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(result, result_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
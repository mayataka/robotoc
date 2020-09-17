#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"


namespace idocp {

class FixedBaseTerminalOCPTest : public ::testing::Test {
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
    s.f_stack() = Eigen::VectorXd::Random(robot.dimf());
    s.set_f();
    s.mu_stack() = Eigen::VectorXd::Random(robot.dim_passive()+robot.dimf());
    s.lmd = Eigen::VectorXd::Random(robot.dimv());
    s.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_tmp= SplitSolution(robot);
    s_tmp.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s_tmp.q);
    s_tmp.v = Eigen::VectorXd::Random(robot.dimv());
    s_tmp.a = Eigen::VectorXd::Random(robot.dimv());
    s_tmp.f_stack() = Eigen::VectorXd::Random(robot.dimf());
    s_tmp.mu_stack() = Eigen::VectorXd::Random(robot.dim_passive()+robot.dimf());
    s_tmp.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_tmp.gmm = Eigen::VectorXd::Random(robot.dimv());
    d = SplitDirection(robot);
    d.dq() = Eigen::VectorXd::Random(robot.dimv());
    d.dv() = Eigen::VectorXd::Random(robot.dimv());
    d.da() = Eigen::VectorXd::Random(robot.dimv());
    d.df() = Eigen::VectorXd::Random(robot.dimf());
    d.du = Eigen::VectorXd::Random(robot.dimv());
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    auto joint_cost = std::make_shared<JointSpaceCost>(robot);
    auto contact_cost = std::make_shared<ContactCost>(robot);
    const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_ref);
    const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    std::vector<Eigen::Vector3d> f_weight, f_ref;
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      f_weight.push_back(Eigen::Vector3d::Random());
      f_ref.push_back(Eigen::Vector3d::Random());
    }
    joint_cost->set_q_weight(q_weight);
    joint_cost->set_q_ref(q_ref);
    joint_cost->set_v_weight(v_weight);
    joint_cost->set_v_ref(v_ref);
    joint_cost->set_a_weight(a_weight);
    joint_cost->set_a_ref(a_ref);
    joint_cost->set_u_weight(u_weight);
    joint_cost->set_u_ref(u_ref);
    joint_cost->set_qf_weight(qf_weight);
    joint_cost->set_vf_weight(vf_weight);
    contact_cost->set_f_weight(f_weight);
    contact_cost->set_f_ref(f_ref);
    cost = std::make_shared<CostFunction>();
    cost->push_back(joint_cost);
    cost->push_back(contact_cost);
    cost_data = CostFunctionData(robot);
    constraints = std::make_shared<Constraints>();
    auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
    auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
    auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
    auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
    constraints->push_back(joint_upper_limit); 
    constraints->push_back(joint_lower_limit);
    constraints->push_back(velocity_lower_limit); 
    constraints->push_back(velocity_upper_limit);
    constraints_data = constraints->createConstraintsData(robot);
    kkt_matrix = KKTMatrix(robot);
    kkt_residual = KKTResidual(robot);
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  std::vector<bool> contact_status;
  std::shared_ptr<CostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<Constraints> constraints;
  ConstraintsData constraints_data;
  SplitSolution s, s_tmp;
  SplitDirection d;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  RiccatiFactorization riccati;
};


TEST_F(FixedBaseTerminalOCPTest, linearizeOCP) {
  s.setContactStatus(robot);
  TerminalOCP ocp(robot, cost, constraints);
  ocp.linearizeOCP(robot, t, s, riccati);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  EXPECT_TRUE(riccati.sq.isApprox(-1*kkt_residual.lq()+s.lmd));
  EXPECT_TRUE(riccati.sv.isApprox(-1*kkt_residual.lv()+s.gmm));
  cost->computeTerminalCostHessian(robot, cost_data, t, s, kkt_matrix);
  EXPECT_TRUE(riccati.Pqq.isApprox(kkt_matrix.Qqq()));
  EXPECT_TRUE(riccati.Pvv.isApprox(kkt_matrix.Qvv()));

  double KKT_ref = 0;
  KKT_ref += (-1*kkt_residual.lq()+s.lmd).squaredNorm();
  KKT_ref += (-1*kkt_residual.lv()+s.gmm).squaredNorm();
  EXPECT_DOUBLE_EQ(KKT_ref, ocp.squaredKKTErrorNorm(robot, t, s));
}


TEST_F(FixedBaseTerminalOCPTest, terminalCost) {
  TerminalOCP ocp(robot, cost, constraints);
  EXPECT_DOUBLE_EQ(ocp.terminalCost(robot, t, s), 
                   cost->phi(robot, cost_data, t, s));
}


TEST_F(FixedBaseTerminalOCPTest, terminalCostWithStepSize) {
  TerminalOCP ocp(robot, cost, constraints);
  const double step_size = 0.3;
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
                   cost->phi(robot, cost_data, t, s_tmp);
  EXPECT_DOUBLE_EQ(ocp.terminalCost(robot, step_size, t, s, d), 
                   cost->phi(robot, cost_data, t, s_tmp));
}


TEST_F(FixedBaseTerminalOCPTest, updatePrimal) {
  TerminalOCP ocp(robot, cost, constraints);
  const double step_size = 0.3;
  riccati.Pqq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati.Pqv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati.Pvq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati.Pvv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati.sq = Eigen::VectorXd::Random(robot.dimv());
  riccati.sv = Eigen::VectorXd::Random(robot.dimv());
  d.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
  s_tmp.lmd = s.lmd + step_size * d.dlmd();
  s_tmp.gmm = s.gmm + step_size * d.dgmm();
  s_tmp.q = s.q + step_size * d.dq();
  s_tmp.v = s.v + step_size * d.dv();
  ocp.updatePrimal(robot, step_size, riccati, d, s);
  EXPECT_TRUE(s.lmd.isApprox(s_tmp.lmd));
  EXPECT_TRUE(s.gmm.isApprox(s_tmp.gmm));
  EXPECT_TRUE(s.q.isApprox(s_tmp.q));
  EXPECT_TRUE(s.v.isApprox(s_tmp.v));
}


TEST_F(FixedBaseTerminalOCPTest, computeSquaredKKTErrorNorm) {
  s.setContactStatus(robot);
  TerminalOCP ocp(robot, cost, constraints);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  kkt_residual.lq() -= s.lmd;
  kkt_residual.lv() -= s.gmm;
  const double kkt_error_ref = kkt_residual.lq().squaredNorm() 
                                + kkt_residual.lv().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error_ref, ocp.computeSquaredKKTErrorNorm(robot, t, s));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
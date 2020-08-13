#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/ocp/kkt_composition.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"


namespace idocp {

class ParNMPCLinearizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(ParNMPCLinearizerTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.set(robot);
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimq());
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf());
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(robot.max_dimf());
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  CostFunctionData data(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  parnmpclinearizer::linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual, lu);
  EXPECT_TRUE(kkt_residual.lq().isApprox(-dtau_*q_weight.asDiagonal()*q_ref));
  EXPECT_TRUE(kkt_residual.lv().isApprox(-dtau_*v_weight.asDiagonal()*v_ref));
  EXPECT_TRUE(kkt_residual.la().isApprox(-dtau_*a_weight.asDiagonal()*a_ref));
  EXPECT_TRUE(kkt_residual.lf().isApprox((-dtau_*f_weight.asDiagonal()*f_ref).head(robot.dimf())));
  EXPECT_TRUE(lu.isApprox(-dtau_*u_weight.asDiagonal()*u_ref));
}


TEST_F(ParNMPCLinearizerTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.set(robot);
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimq());
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf());
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(robot.max_dimf());
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  CostFunctionData data(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  parnmpclinearizer::linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual, lu);
  EXPECT_TRUE(kkt_residual.lq().isApprox(-dtau_*q_weight.asDiagonal()*q_ref));
  EXPECT_TRUE(kkt_residual.lv().isApprox(-dtau_*v_weight.asDiagonal()*v_ref));
  EXPECT_TRUE(kkt_residual.la().isApprox(-dtau_*a_weight.asDiagonal()*a_ref));
  EXPECT_TRUE(kkt_residual.lf().isApprox((-dtau_*f_weight.asDiagonal()*f_ref).head(robot.dimf())));
  EXPECT_TRUE(lu.isApprox(-dtau_*u_weight.asDiagonal()*u_ref));
  std::cout << kkt_residual.lf().transpose() << std::endl;
  std::cout << (-dtau_*f_weight.asDiagonal()*f_ref).head(robot.dimf()).transpose() << std::endl;
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
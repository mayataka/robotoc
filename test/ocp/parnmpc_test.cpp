#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/ocp/kkt_composition.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/ocp/parnmpc.hpp"


namespace idocp {

class ParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    T_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    N_ = 10;
    num_proc_ = 1;
    dtau_ = T_ / N_;
  }

  virtual void TearDown() {
  }

  double dtau_, t_, T_;
  int N_, num_proc_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(ParNMPCTest, updateSolutionFixedBase) {
  std::vector<int> contact_frames = {18};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  // Robot robot(fixed_base_urdf_, contact_frames, baum_a, baum_b);
  Robot robot(fixed_base_urdf_);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  // robot.setContactStatus(contact_status);
  // std::vector<std::vector<bool>> contact_sequence = std::vector<std::vector<bool>>(N_, contact_status);
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.generateFeasibleConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.1);
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Zero(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Zero(robot.dimv());
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Constant(robot.max_dimf(), 0.01);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(robot.max_dimf());
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_qf_weight(qf_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_vf_weight(vf_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  std::shared_ptr<Constraints> constraints = std::make_shared<Constraints>();
  ParNMPC parnmpc(robot, cost, constraints, T_, N_, num_proc_);
  // parnmpc.setContactSequence(contact_sequence);
  Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::cout << "initial KKT error = " << parnmpc.KKTError(t_, q, v) << std::endl;
  const int num_itr = 10;
  for (int i=0; i<num_itr; ++i) {
    parnmpc.updateSolution(t_, q, v, false);
    std::cout << "KKT error after " << (i+1) << "iteration = " << parnmpc.KKTError(t_, q, v) << std::endl;
  }
}


TEST_F(ParNMPCTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  Robot robot(floating_base_urdf_, contact_frames, baum_a, baum_b);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  std::vector<std::vector<bool>> contact_sequence = std::vector<std::vector<bool>>(N_, contact_status);
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.generateFeasibleConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.1);
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Zero(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Constant(robot.dimv(), 0.1);
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Constant(robot.max_dimf(), 0.01);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(robot.max_dimf());
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_qf_weight(qf_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_vf_weight(vf_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  std::shared_ptr<Constraints> constraints = std::make_shared<Constraints>();
  ParNMPC parnmpc(robot, cost, constraints, T_, N_, num_proc_);
  parnmpc.setContactSequence(contact_sequence);
  Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::cout << "initial KKT error = " << parnmpc.KKTError(t_, q, v) << std::endl;
  const int num_itr = 10;
  for (int i=0; i<num_itr; ++i) {
    parnmpc.updateSolution(t_, q, v, false);
    std::cout << "KKT error after " << (i+1) << "iteration = " << parnmpc.KKTError(t_, q, v) << std::endl;
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
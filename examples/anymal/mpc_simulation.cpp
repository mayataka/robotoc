#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/parnmpc.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

#include "idocp/utils/joint_constraints_factory.hpp"
#include "idocp/utils/quadruped_simulator.hpp"

namespace mpcsimulation {
namespace anymal {

void SimulateWithContactsByOCP() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames);
  auto cost = std::make_shared<idocp::CostFunction>();
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  Eigen::VectorXd q_ref(robot.dimq());
  // q_ref << 0, 0, 0.48, 0, 0, 0, 1, 
  // q_ref << 0, 0, 0.48, 0, -50, 0, 1, 
  // q_ref << 0, 0, 0.48, 0, 50, 0, 1, 
  q_ref << 0, 0, 0.48, -1, 0, 1, 0.5, 
  // q_ref << 0, 0, 0.48, 1, 0, -1, 0.5, 
           0.0315, 0.4, -0.8, 
           0.0315, -0.4, 0.8, 
           -0.0315, 0.4, -0.8,
           -0.0315, -0.4, 0.8;
  robot.normalizeConfiguration(q_ref);
  std::cout << q_ref.transpose() << std::endl;
  joint_cost->set_q_ref(q_ref);
  Eigen::VectorXd q_weight(robot.dimv());
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.0));
  cost->push_back(joint_cost);
  // std::vector<Eigen::Vector3d> f_weight;
  // for (int i=0; i<contact_frames.size(); ++i) {
  //   f_weight.push_back(Eigen::Vector3d::Constant(0.0));
  // }
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();
  const double T = 1;
  const int N = 20;
  const int num_proc = 4;
  const double t = 0;
  Eigen::VectorXd q(robot.dimq());
  q << 0, 0, 0.48, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactStatus({true, true, true, true});
  robot.updateKinematics(q, v, Eigen::VectorXd::Zero(robot.dimv()));
  robot.setContactPointsByCurrentKinematics();
  idocp::MPC<idocp::OCP> mpc(robot, cost, constraints, T, N, num_proc);
  mpc.activateContacts({0, 1, 2, 3}, 0, N);
  mpc.setContactPointByKinematics(q);
  mpc.initializeSolution(t, q, v, 100);
  const std::string urdf_for_raisim_file_name = "../anymal/anymal_for_raisim.urdf";
  idocp::QuadrupedSimulator simulator(urdf_for_raisim_file_name, "../sim_result", "forward");
  simulator.viz(mpc, 5, 0.0025, 0, q, v, false);
}

} // namespace anymal 
} // namespace ocpbenchmark


int main() {
  mpcsimulation::anymal::SimulateWithContactsByOCP();
  return 0;
}

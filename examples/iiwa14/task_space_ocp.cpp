#include <string>
#include <memory>

#include "Eigen/Core"

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/solver/unconstr_ocp_solver.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/task_space_6d_cost.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "robotoc/utils/trajectory_viewer.hpp"
#endif 


class TaskSpace6DRef final : public robotoc::TaskSpace6DRefBase {
public:
  TaskSpace6DRef() 
    : TaskSpace6DRefBase() {
    rotm_  <<  0, 0, 1, 
               0, 1, 0,
              -1, 0, 0;
    pos0_ << 0.546, 0, 0.76;
    radius_ = 0.05;
  }

  ~TaskSpace6DRef() {}

  void updateRef(const robotoc::GridInfo& grid_info, robotoc::SE3& ref) const override {
    Eigen::Vector3d pos(pos0_);
    pos.coeffRef(1) += radius_ * sin(M_PI*grid_info.t);
    pos.coeffRef(2) += radius_ * cos(M_PI*grid_info.t);
    ref = robotoc::SE3(rotm_, pos);
  }

  bool isActive(const robotoc::GridInfo& grid_info) const override {
    return true;
  }

private:
  double radius_;
  Eigen::Matrix3d rotm_;
  Eigen::Vector3d pos0_;
};


int main(int argc, char *argv[]) {
  // Create a robot.
  robotoc::RobotModelInfo model_info;
  model_info.urdf_path = "../iiwa_description/urdf/iiwa14.urdf";
  robotoc::Robot robot(model_info);

  const std::string ee_frame = "iiwa_link_ee_kuka";

  // Change the limits from the default parameters.
  robot.setJointEffortLimit(Eigen::VectorXd::Constant(robot.dimu(), 50));
  robot.setJointVelocityLimit(Eigen::VectorXd::Constant(robot.dimv(), M_PI_2));

  // Create a cost function.
  auto cost = std::make_shared<robotoc::CostFunction>();
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  config_cost->set_q_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.0001));
  config_cost->set_v_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 0.0001));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.0001));
  cost->push_back(config_cost);
  auto x6d_ref = std::make_shared<TaskSpace6DRef>();
  auto task_cost = std::make_shared<robotoc::TaskSpace6DCost>(robot, ee_frame, x6d_ref);
  task_cost->set_weight(Eigen::Vector3d::Constant(1000), Eigen::Vector3d::Constant(1000));
  task_cost->set_weight_terminal(Eigen::Vector3d::Constant(1000), Eigen::Vector3d::Constant(1000));
  cost->push_back(task_cost);

  // Create joint constraints.
  const double barrier_param = 1.0e-03;
  const double fraction_to_boundary_rule = 0.995;
  auto constraints = std::make_shared<robotoc::Constraints>(barrier_param, fraction_to_boundary_rule);
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);

  // Create the OCP solver for unconstrained rigid-body systems.
  const double T = 6;
  const int N = 120;
  robotoc::OCP ocp(robot, cost, constraints, T, N);
  auto solver_options = robotoc::SolverOptions();
  solver_options.nthreads = 4;
  robotoc::UnconstrOCPSolver ocp_solver(ocp, solver_options);

  // Initial time and initial state
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, M_PI_2, 0, M_PI_2, 0, M_PI_2, 0;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());

  // Solves the OCP.
  ocp_solver.discretize(t);
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  ocp_solver.initConstraints();
  std::cout << "Initial KKT error: " << ocp_solver.KKTError(t, q, v) << std::endl;
  ocp_solver.solve(t, q, v);
  std::cout << "KKT error after convergence: " << ocp_solver.KKTError(t, q, v) << std::endl;
  std::cout << ocp_solver.getSolverStatistics() << std::endl;

#ifdef ENABLE_VIEWER
  robotoc::TrajectoryViewer viewer(path_to_urdf);
  const double dt = T/N;
  viewer.display(ocp_solver.getSolution("q"), dt);
#endif 

  return 0;
}
#include "ocp/MPC.hpp"


namespace idocp {


MPC::MPC(const Robot& robot, const CostFunctionInterface* cost,
         const ConstraintsInterface* constraints, const double T, 
         const unsigned int N, const unsigned int num_proc) 
  : ocp_(robot, cost, constraints, T, N, num_proc) {
}


void MPC::initializeSolution(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, 
                             const unsigned int max_itr) {
  ocp_.setStateTrajectory(q, v);
  for (int i=0; i<max_iter; ++i) {
    ocp_.solveSQP(t, q, v, true);
  }
}


void MPC::updateSolution(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v) {
  ocp_.solveSQP(t, q, v, false);
}


void MPC::getControlInput(Eigen::VectorXd& u) {
  ocp_.getInitialControlInput(u);
}


double MPC::optimalityError(const double t, const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v) {
  ocp_.optimalityError(t, q, v);
}

} // namespace idocp
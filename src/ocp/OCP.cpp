#include "ocp/OCP.hpp"

#include <assert.h>


namespace invdynocp {

OCP::OCP(const Robot& robot, const CostFunctionInterface* cost,
         const ConstraintsInterface* constraints, const double T, 
         const unsigned int N, const unsigned int num_proc)
  : split_ocps_(N, SplitOCP(robot, cost, constraints)),
    robots_(num_proc_, robot),
    cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc),
    x_(N+1, Eigen::VectorXd::Zero(robot.dimq()+robot.dimv())),
    a_(N, Eigen::VectorXd::Zero(robot.dimv())),
    lmd_(N+1, Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(N+1, Eigen::VectorXd::Zero(2*robot.dimv())),
    da_(N, Eigen::VectorXd::Zero(robot.dimv())),
    dlmd_(N+1, Eigen::VectorXd::Zero(2*robot.dimv())),
    s_(N+1, Eigen::VectorXd::Zero(2*robot.dimv())),
    P_(N+1, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
  if (num_proc_ <= 0) {
    num_proc_ = 1;
  }
}


void OCP::solveSQP(const double t, const Eigen::VectorXd& x) {
  int time_step;
  #pragma omp parallel num_threads(num_proc_)
  {
    #pragma omp for
    for (time_step=0; time_step<N_; ++time_step) {
      split_ocps_[time_step].linearizeOCP(robots_[omp_get_num_threads()], 
                                          t+time_step*dtau_, dtau_, 
                                          lmd_[time_step], x_[time_step], 
                                          a_[time_step], lmd_[time_step+1], 
                                          x_[time_step+1]);
    }
  }
  split_ocps_[N_].linearizeTerminalCost(robots_[num_proc_-1], t+T_, x_[N_], 
                                        P_[N_], s_[N_]);
  for (time_step=N_-1; time_step>=0; --time_step) {
    split_ocps_[time_step].backwardRecursion(dtau_, P_[time_step+1], 
                                             s_[time_step+1], P_[time_step], 
                                             s_[time_step]);
  }
  dx_[0].setZero();
  for (time_step=0; time_step<N_; ++time_step) {
    split_ocps_[time_step].forwardRecursion(dtau_, dx_[time_step], 
                                            da_[time_step],
                                            dx_[time_step+1]);
  }
  lmd_[0].noalias() -= s_[0];
  #pragma omp parallel num_threads(num_proc_)
  {
    #pragma omp for
    for (time_step=0; time_step<N_; ++time_step) {
      split_ocps_[time_step].updateOCP(robots_[omp_get_num_threads()], 
                                       dx_[time_step], da_[time_step], 
                                       P_[time_step+1], s_[time_step+1], 
                                       x_[time_step], a_[time_step], 
                                       lmd_[time_step+1]);
    }
  }
  x_[N_].noalias() += dx_[N_];
}

} // namespace invdynocp
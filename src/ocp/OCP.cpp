#include "ocp/OCP.hpp"

#include <assert.h>


namespace invdynocp {

OCP::OCP(const Robot& robot, const CostFunctionInterface* cost,
         const ConstraintsInterface* constraints, const double T, 
         const unsigned int N, const unsigned int num_proc)
  : split_OCPs_(N, SplitOCP(robot, cost, constraints)),
    robots_(num_proc, robot),
    cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc),
    q_(N+1, Eigen::VectorXd::Zero(robot.dimq())),
    v_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    a_(N, Eigen::VectorXd::Zero(robot.dimv())),
    lmd_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    gmm_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    dq_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    dv_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    da_(N, Eigen::VectorXd::Zero(robot.dimv())),
    dlmd_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    dgmm_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    sq_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    sv_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    Pqq_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pqv_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pvq_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pvv_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
  if (num_proc_ == 0) {
    num_proc_ = 1;
  }
}


void OCP::solveSQP(const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v) {
  int time_step;
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for  
    for (time_step=0; time_step<N_; ++time_step) {
      split_OCPs_[time_step].linearizeOCP(robots_[omp_get_thread_num()], 
                                          t+time_step*dtau_, dtau_, 
                                          lmd_[time_step], gmm_[time_step],
                                          q_[time_step], v_[time_step], 
                                          a_[time_step], 
                                          lmd_[time_step+1], gmm_[time_step+1], 
                                          q_[time_step+1], v_[time_step+1]);
    }
  }
  split_OCPs_[N_-1].linearizeTerminalCost(robots_[num_proc_-1], t+T_, q_[N_], 
                                          v_[N_], Pqq_[N_], Pqv_[N_], Pvq_[N_], 
                                          Pvv_[N_], sq_[N_], sv_[N_]);
  for (time_step=N_-1; time_step>=0; --time_step) {
    split_OCPs_[time_step].backwardRecursion(dtau_, Pqq_[time_step+1], 
                                             Pqv_[time_step+1], 
                                             Pvq_[time_step+1], 
                                             Pvv_[time_step+1], sq_[time_step+1], 
                                             sv_[time_step+1], Pqq_[time_step], 
                                             Pqv_[time_step], Pvq_[time_step], 
                                             Pvv_[time_step], sq_[time_step], 
                                             sv_[time_step]);
  }
  dq_[0] = q - q_[0];
  dv_[0] = v - v_[0];
  for (time_step=0; time_step<N_; ++time_step) {
    split_OCPs_[time_step].forwardRecursion(dtau_, dq_[time_step], 
                                            dv_[time_step], da_[time_step],
                                            dq_[time_step+1], dv_[time_step+1]);
  }
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for 
    for (time_step=0; time_step<N_; ++time_step) {
      if (time_step < N_) {
        split_OCPs_[time_step].updateOCP(robots_[omp_get_thread_num()], 
                                        dq_[time_step], dv_[time_step], 
                                        da_[time_step], Pqq_[time_step], 
                                        Pqv_[time_step], Pvq_[time_step], 
                                        Pvv_[time_step], sq_[time_step], 
                                        sv_[time_step], q_[time_step], 
                                        v_[time_step], a_[time_step], 
                                        lmd_[time_step], gmm_[time_step]);
      }
      else {
        split_OCPs_[time_step].updateOCP(robots_[omp_get_thread_num()], 
                                        dq_[time_step], dv_[time_step], 
                                        Pqq_[time_step], Pqv_[time_step], 
                                        Pvq_[time_step], Pvv_[time_step], 
                                        sq_[time_step], sv_[time_step], 
                                        q_[time_step], v_[time_step], 
                                        lmd_[time_step], gmm_[time_step]);
      }
    }
  }
}


void OCP::getInitialControlInput(Eigen::VectorXd& u) {
  robots_[0].RNEA(q_[0], v_[0], a_[0], u);
}


double OCP::optimalityError(const double t, const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v) {
  double error = 0;
  error += (q-q_[0]).norm();
  error += (v-v_[0]).norm();
  for (int i=0; i<N_; ++i) {
    error += split_OCPs_[i].optimalityError(robots_[0], t+i*dtau_, dtau_, 
                                            lmd_[i], gmm_[i], q_[i], v_[i], 
                                            a_[i], lmd_[i+1], gmm_[i+1], 
                                            q_[i+1], v_[i+1]);
  }
  error += split_OCPs_[N_-1].terminalError(robots_[0], t+T_, lmd_[N_], gmm_[N_], 
                                           q_[N_], v_[N_]);
  return error;
}


void OCP::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "time step: " << i << std::endl;
    std::cout << "q: " << q_[i].transpose() << std::endl;
    std::cout << "v: " << v_[i].transpose() << std::endl;
    std::cout << "a: " << a_[i].transpose() << std::endl;
  }
  std::cout << "time step: " << N_ << std::endl;
  std::cout << "q: " << q_[N_].transpose() << std::endl;
  std::cout << "v: " << v_[N_].transpose() << std::endl;
}

} // namespace invdynOCP
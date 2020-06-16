#include "ocp/OCP.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {

OCP::OCP(const Robot& robot, const CostFunctionInterface* cost,
         const ConstraintsInterface* constraints, const double T, 
         const unsigned int N, const unsigned int num_proc)
  : split_OCPs_(N+1, SplitOCP(robot, cost, constraints)),
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
    u_(N, Eigen::VectorXd::Zero(robot.dimv())),
    beta_(N, Eigen::VectorXd::Zero(robot.dimv())),
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
    Pvv_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    step_size_(Eigen::VectorXd::Ones(N)) {
  if (num_proc_ == 0) {
    num_proc_ = 1;
  }
  int time_step;
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for  
    for (time_step=0; time_step<N_; ++time_step) {
      split_OCPs_[time_step].initConstraints(robots_[omp_get_thread_num()], 
                                             dtau_, q_[time_step], 
                                             v_[time_step], a_[time_step], 
                                             u_[time_step]);
    }
  }
}


void OCP::solveLQR(const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v) {
  int time_step;
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for  
    for (time_step=0; time_step<N_; ++time_step) {
      if (time_step < N_) {
        split_OCPs_[time_step].linearizeOCP(robots_[omp_get_thread_num()], 
                                            t+time_step*dtau_, dtau_, 
                                            lmd_[time_step], gmm_[time_step],
                                            q_[time_step], v_[time_step], 
                                            a_[time_step], u_[time_step], 
                                            beta_[time_step],
                                            lmd_[time_step+1], 
                                            gmm_[time_step+1], q_[time_step+1], 
                                            v_[time_step+1]);
      }
      else {
        split_OCPs_[N_].linearizeOCP(robots_[omp_get_thread_num()], t+T_, 
                                     lmd_[N_], gmm_[N_], q_[N_], v_[N_], 
                                     Pqq_[N_], Pqv_[N_], Pvq_[N_], Pvv_[N_], 
                                     sq_[N_], sv_[N_]);
      }
    }
  }
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
      step_size_[time_step] 
          = split_OCPs_[time_step].computeMaxStepSize(
              robots_[omp_get_thread_num()], dtau_, dq_[time_step], 
              dv_[time_step], da_[time_step]);
    }
  }
  const double step_size = step_size_.minCoeff();
  // std::cout << "step size = " << step_size_.transpose() << std::endl;
  // std::cout << "step size = " << step_size << std::endl;
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for 
    for (time_step=0; time_step<=N_; ++time_step) {
      if (time_step < N_) {
        split_OCPs_[time_step].updateOCP(robots_[omp_get_thread_num()], 
                                         step_size, dtau_, dq_[time_step], 
                                         dv_[time_step], da_[time_step], 
                                         Pqq_[time_step], Pqv_[time_step], 
                                         Pvq_[time_step], Pvv_[time_step], 
                                         sq_[time_step], sv_[time_step], 
                                         q_[time_step], v_[time_step], 
                                         a_[time_step], u_[time_step], 
                                         beta_[time_step], lmd_[time_step], 
                                         gmm_[time_step]);
      }
      else {
        split_OCPs_[N_].updateOCP(robots_[omp_get_thread_num()], step_size, 
                                  dq_[N_], dv_[N_], Pqq_[N_], Pqv_[N_], 
                                  Pvq_[N_], Pvv_[N_], sq_[N_], sv_[N_], 
                                  q_[N_], v_[N_], lmd_[N_], gmm_[N_]);
      }
    }
  }
  // #pragma omp parallel num_threads(num_proc_) 
  // {
  //   #pragma omp for 
  //   for (time_step=0; time_step<=N_; ++time_step) {
  //     if (time_step < N_) {
  //       split_OCPs_[time_step].updateOCP(robots_[omp_get_thread_num()], 
  //                                        dtau_, dq_[time_step], dv_[time_step], 
  //                                        da_[time_step], Pqq_[time_step], 
  //                                        Pqv_[time_step], Pvq_[time_step], 
  //                                        Pvv_[time_step], sq_[time_step], 
  //                                        sv_[time_step], q_[time_step], 
  //                                        v_[time_step], a_[time_step], 
  //                                        u_[time_step], beta_[time_step],
  //                                        lmd_[time_step], gmm_[time_step]);
  //     }
  //     else {
  //       split_OCPs_[N_].updateOCP(robots_[omp_get_thread_num()], dq_[N_], 
  //                                 dv_[N_], Pqq_[N_], Pqv_[N_], Pvq_[N_], 
  //                                 Pvv_[N_], sq_[N_], sv_[N_], q_[N_], v_[N_], 
  //                                 lmd_[N_], gmm_[N_]);
  //     }
  //   }
  // }
}


void OCP::getInitialControlInput(Eigen::VectorXd& u) {
  u = u_[0];
}

void OCP::setStateTrajectory(const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  for (int i=0; i<=N_; ++i) {
    v_[i] = v;
  }
  for (int i=0; i<=N_; ++i) {
    q_[i] = q;
  }
}


void OCP::setStateTrajectory(const Eigen::VectorXd& q0, 
                             const Eigen::VectorXd& v0, 
                             const Eigen::VectorXd& qN, 
                             const Eigen::VectorXd& vN) {
  Eigen::VectorXd dv = (vN-v0) / N_;
  for (int i=0; i<=N_; ++i) {
    v_[i] = v0 + i * dv;
  }
  Eigen::VectorXd dq = (qN-q0) / N_;
  for (int i=0; i<=N_; ++i) {
    q_[i] = q0 + i * dq;
  }
}


double OCP::optimalityError(const double t, const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v) {
  double error = 0;
  error += (q-q_[0]).squaredNorm();
  error += (v-v_[0]).squaredNorm();
  for (int i=0; i<N_; ++i) {
    error += split_OCPs_[i].squaredOCPErrorNorm(robots_[0], t+i*dtau_, dtau_, 
                                                lmd_[i], gmm_[i], q_[i], v_[i], 
                                                a_[i], u_[i], beta_[i], 
                                                lmd_[i+1], gmm_[i+1], 
                                                q_[i+1], v_[i+1]);
  }
  error += split_OCPs_[N_].squaredOCPErrorNorm(robots_[0], t+T_, lmd_[N_], 
                                               gmm_[N_], q_[N_], v_[N_]);
  return std::sqrt(error);
}


void OCP::printSolution() {
  for (int i=0; i<N_; ++i) {
    std::cout << "time step: " << i << std::endl;
    std::cout << "q: " << q_[i].transpose() << std::endl;
    std::cout << "v: " << v_[i].transpose() << std::endl;
    // std::cout << "a: " << a_[i].transpose() << std::endl;
    std::cout << "u: " << u_[i].transpose() << std::endl;
  }
  std::cout << "time step: " << N_ << std::endl;
  std::cout << "q: " << q_[N_].transpose() << std::endl;
  std::cout << "v: " << v_[N_].transpose() << std::endl;
}

} // namespace idocp
#include "idocp/ocp/ocp.hpp"

#include <cmath>
#include <omp.h>
#include <cassert>


namespace idocp {

OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
         const std::shared_ptr<Constraints>& constraints, const double T, 
         const int N, const int max_num_impulse, const int num_proc)
  : split_ocps_(N, SplitOCP(robot, cost, constraints), max_num_impulse, 
                SplitImpulseOCP(robot, cost->getImpulseCostFunction(), 
                                constraints->getImpulseConstraints())),
    terminal_ocp_(robot, cost, constraints),
    robots_(num_proc, robot),
    filter_(),
    T_(T),
    dtau_(T/N),
    step_size_reduction_rate_(0.75),
    min_step_size_(0.05),
    N_(N),
    num_proc_(num_proc),
    s_(N+1, SplitSolution(robot), max_num_impulse, ImpulseSplitSolution(robot)),
    d_(N+1, SplitDirection(robot), max_num_impulse, ImpulseSplitDirection(robot)),
    riccati_(N+1, RiccatiFactorization(robot), max_num_impulse, RiccatiFactorization(robot)),
    constraint_factorization_(max_num_impulse, 
                              StateConstraintRiccatiFactorization(robot, N, max_num_impulse)),
    constraint_factorizer_(robot, max_num_impulse, num_proc),
    riccati_recursion_(robot, T, N, max_num_impulse, num_proc),
    primal_step_sizes_(Eigen::VectorXd::Zero(N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N)),
    contact_sequence_(robot, T, N) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N; ++i) {
    robot.normalizeConfiguration(s_[i].q);
  }
  initConstraints();
}


OCP::OCP() 
  : split_ocps_(),
    terminal_ocp_(),
    robots_(),
    contact_sequence_(),
    filter_(),
    T_(),
    dtau_(),
    step_size_reduction_rate_(),
    min_step_size_(),
    N_(),
    num_proc_(),
    s_(),
    d_(),
    riccati_(),
    primal_step_sizes_(),
    dual_step_sizes_() {
}


OCP::~OCP() {
}


void OCP::updateSolution(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v, const bool use_line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());

  const int N_impulse = contact_sequence_.totalNumImpulseStages();
  const int N_lift = contact_sequence_.totalNumLiftStages();
  const int N_all = N_ + 1 + N_impulse + N_lift;

  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i == 0) {
      const int robot_id = omp_get_thread_num();
      if (contact_sequence_.existImpulse(i)) {
        split_ocps_[0].linearizeOCP(robots_[robot_id], 
                                    contact_sequence_.contactStatus(0), 
                                    t+i*dtau_, dtau_, q, s_[0], s_[1],
                                    kkt_matrix_[0], kkt_residual_[0]);
      }
      else if (contact_sequence_.existLift(i)) {
        split_ocps_[0].linearizeOCP(robots_[robot_id], 
                                    contact_sequence_.contactStatus(0), 
                                    t+i*dtau_, dtau_, q, s_[0], s_.impulse[0],
                                    kkt_matrix_[0], kkt_residual_[0]);
      }
      else {
        split_ocps_[0].linearizeOCP(robots_[robot_id], 
                                    contact_sequence_.contactStatus(0), 
                                    t+i*dtau_, dtau_, q, s_[0], s_.lift[0],
                                    kkt_matrix_[0], kkt_residual_[0]);
      }
    }
    else if (i < N_) {
      const int robot_id = omp_get_thread_num();
      if (contact_sequence_.existImpulse(i)) {
        split_ocps_[i].linearizeOCP(robots_[robot_id], 
                                    contact_sequence_.contactStatus(i), 
                                    t+i*dtau_, dtau_, q_prev(i), s_[i], 
                                    s_.impulse[contact_sequence_.impulseIndex(i)],
                                    kkt_matrix_[0], kkt_residual_[0]);
      }
      else if (contact_sequence_.existLift(i)) {
        split_ocps_[i].linearizeOCP(robots_[robot_id], 
                                    contact_sequence_.contactStatus(i), 
                                    t+i*dtau_, dtau_, q_prev(i), s_[i], 
                                    s_.impulse[contact_sequence_.liftIndex(i)],
                                    kkt_matrix_[i], kkt_residual_[i]);
      }
      else {
        split_ocps_[i].linearizeOCP(robots_[robot_id], 
                                    contact_sequence_.contactStatus(i), 
                                    t+i*dtau_, dtau_, q_prev(i), s_[i], s_.lift[i],
                                    kkt_matrix_[i], kkt_residual_[i]);
      }
    }
    else if (i == N_) {
      const int robot_id = omp_get_thread_num();
      terminal_ocp_.linearizeOCP(robots_[robot_id], t+T_, s_[N_]);
      riccati_recursion_.backwardRiccatiRecursionTerminal(
          kkt_matrix_, kkt_residual_, riccati_factorization_);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int robot_id = omp_get_thread_num();
      const int impulse_idx = i - (N_+1);
      const int time_stage_before_impulse = contact_sequence_.timeStageBeforeImpulse(impulse_idx);
      split_ocps_.impulse[i].linearizeOCP(robots_[robot_id], 
                                          contact_sequence_.impulseStatus(i), 
                                          t+contact_sequence_.impulseTime(impulse_idx),  
                                          s_[time_stage_before_impulse].q, 
                                          s_.impulse[i], s_.aux[i]);
      split_ocps_.impulse[i].getStateConstraintFactorization(
          constraint_factorization_[impulse_idx].Eq(), 
          constraint_factorization_[impulse_idx].e());
      constraint_factorization_[impulse_idx].T_impulse(impulse_idx).topRows(robots_[robot_id].dimv())
          = constraint_factorization_[impulse_idx].Eq();
      constraint_factorization_[impulse_idx].T_impulse(impulse_idx).bottomRows(robots_[robot_id].dimv()).setZero();
    }
    else {
      const int robot_id = omp_get_thread_num();
      const int lift_idx = i - (N_+1+N_impulse);
      const int time_stage_before_lift = contact_sequence_.timeStageBeforeLift(lift_idx);
      split_ocps_.lift[i].linearizeOCP(robots_[robot_id], 
                                       contact_sequence_.contactStatus(time_stage_before_lift+1), 
                                       t+contact_sequence_.liftTime(lift_idx),  
                                       dtau_,
                                       s_[time_stage_before_lift].q, 
                                       s_.lift[i],  
                                       s_[time_stage_before_lift+1]);

    }
  }

  riccati_recursion_.backwardRiccatiRecursion(contact_sequence_, kkt_matrix_, 
                                              kkt_residual_, 
                                              riccati_factorization_);
  riccati_recursion_.forwardRiccatiRecursionParallel(contact_sequence_, 
                                                     kkt_matrix_, 
                                                     kkt_residual_);
  riccati_recursion_.forwardRiccatiRecursionSerial(contact_sequence_, 
                                                   kkt_matrix_, kkt_residual_, 
                                                   riccati_factorization_);
  riccati_recursion_.backwardStateConstraintFactorization(kkt_matrix_, 
                                                          constraint_factorization_);
  robots_[0].subtractConfiguration(q, s_[0].q, d_[0].dq());
  d_[0].dv() = v - s_[0].v;
  constraint_factorization_.computeLagrangeMultiplierDirection(
      contact_sequence_, riccati_factorization_.impulse, 
      constraint_factorization_, d_[0].dx(), d_);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      const int robot_id = omp_get_thread_num();
      SplitRiccatiFactorizer::computeStateDirection();
      SplitRiccatiFactorizer::computeCostateDirection();
      SplitRiccatiFactorizer::computeControlInputDirection();
      split_ocps_[i].computeCondensedPrimalDirection(robots_[robot_id], dtau_, d_[i+1], d_[i]);
      
    }
    else {
      const int robot_id = omp_get_thread_num();
      terminal_ocp_.updatePrimal(robots_[robot_id], primal_step_size, d_[N_], s_[N_]);
      terminal_ocp_.updateDual(dual_step_size);
    }
  }

  const double primal_step_size = riccati_recursion_.primalStepSize(contact_sequence_);
  const double dual_step_size = riccati_recursion_.dualStepSize(contact_sequence_);

  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      const int robot_id = omp_get_thread_num();
      split_ocps_[i].computeCondensedDualDirection(robots_[robot_id], dtau_, d_[i+1], d_[i]);
      split_ocps_[i].updatePrimal(robots_[robot_id], primal_step_size, dtau_, d_[i], s_[i]);
      split_ocps_[i].updateDual(dual_step_size);
    }
    else {
      const int robot_id = omp_get_thread_num();
      terminal_ocp_.updatePrimal(robots_[robot_id], primal_step_size, d_[N_], s_[N_]);
      terminal_ocp_.updateDual(dual_step_size);
    }
  }
} 


void OCP::setProblem() {
  
}


void OCP::linearizeSplitOCPs(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  // const int N_impulse = contact_sequence_.totalNumImpulseStages();
  // const int N_lift = contact_sequence_.totalNumLiftStages();
  // const int N_all = N_ + 1 + N_impulse + N_lift;
  // #pragma omp parallel for num_threads(num_proc_)
  // for (int i=0; i<N_all; ++i) {
  //   if (i == 0) {
  //     const int robot_id = omp_get_thread_num();
  //     if (contact_sequence_.existImpulse(0)) {
  //       split_ocps_[0].linearizeOCP(robots_[robot_id], 
  //                                   contact_sequence_.contactStatus(0), 
  //                                   t, dtau_, q, s_[0], s_.aux[0]);
  //     }
  //     else if (contact_sequence_.existLift(0)) {
  //       split_ocps_[0].linearizeOCP(robots_[robot_id], 
  //                                   contact_sequence_.contactStatus(0), 
  //                                   t, dtau_, q, s_[0], s_.lift[0]);
  //     }
  //     else {
  //       split_ocps_[0].linearizeOCP(robots_[robot_id], 
  //                                   contact_sequence_.contactStatus(0), 
  //                                   t, dtau_, q, s_[0], s_[1]);
  //     }
  //   }
  //   else if (i < N_) {
  //     const int robot_id = omp_get_thread_num();
  //     if (contact_sequence_.existImpulse(i)) {
  //       split_ocps_[i].linearizeOCP(robots_[robot_id], 
  //                                   contact_sequence_.contactStatus(i), 
  //                                   t, dtau_, q_prev(i), s_[i], 
  //                                   s_.aux[contact_sequence_.impulseIndex(i)]);
  //     }
  //     else if (contact_sequence_.existLift(i)) {
  //       split_ocps_[i].linearizeOCP(robots_[robot_id], 
  //                                   contact_sequence_.contactStatus(i), 
  //                                   t, dtau_, q_prev(i), s_[i], 
  //                                   s_.lift[contact_sequence_.liftIndex(i)]);
  //     }
  //     else {
  //       split_ocps_[i].linearizeOCP(robots_[robot_id], 
  //                                   contact_sequence_.contactStatus(i), 
  //                                   t, dtau_, q_prev(i), s_[i], s_[i]);
  //     }
  //   }
  //   else if (i == N_) {
  //     const int robot_id = omp_get_thread_num();
  //     terminal_ocp_.linearizeOCP(robots_[robot_id], t+T_, s_[N_]);
  //     terminal_ocp_.backwardRiccatiRecursion(riccati_[N_]);
  //   }
  //   else if (i < N_ + 1 + N_impulse) {
  //     const int robot_id = omp_get_thread_num();
  //     const int impulse_idx = i - (N_+1);
  //     const int time_stage_before_impulse = contact_sequence_.timeStageBeforeImpulse(impulse_idx);
  //     split_ocps_.impulse[i].linearizeOCP(robots_[robot_id], 
  //                                         contact_sequence_.impulseStatus(i), 
  //                                         t+contact_sequence_.impulseTime(impulse_idx),  
  //                                         s_[time_stage_before_impulse].q, 
  //                                         s_.impulse[i], s_.aux[i]);
  //     split_ocps_.impulse[i].getStateConstraintFactorization(
  //         constraint_factorization_[impulse_idx].Eq(), 
  //         constraint_factorization_[impulse_idx].e());
  //     constraint_factorization_[impulse_idx].T_impulse(impulse_idx).topRows(robots_[robot_id].dimv())
  //         = constraint_factorization_[impulse_idx].Eq();
  //     constraint_factorization_[impulse_idx].T_impulse(impulse_idx).bottomRows(robots_[robot_id].dimv()).setZero();
  //   }
  //   else {
  //     const int robot_id = omp_get_thread_num();
  //     const int lift_idx = i - (N_+1+N_impulse);
  //     const int time_stage_before_lift = contact_sequence_.timeStageBeforeLift(lift_idx);
  //     split_ocps_.lift[i].linearizeOCP(robots_[robot_id], 
  //                                      contact_sequence_.contactStatus(time_stage_before_lift+1), 
  //                                      t+contact_sequence_.liftTime(lift_idx),  
  //                                      dtau_,
  //                                      s_[time_stage_before_lift].q, 
  //                                      s_.lift[i],  
  //                                      s_[time_stage_before_lift+1]);
  //   }
  // }

}


const SplitSolution& OCP::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= N_);
  return s_[stage];
}


void OCP::getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                               Eigen::MatrixXd& Kv) const {
  assert(stage >= 0);
  assert(stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  split_ocps_[stage].getStateFeedbackGain(Kq, Kv);
}


bool OCP::setStateTrajectory(const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  Eigen::VectorXd q_normalized = q;
  robots_[0].normalizeConfiguration(q_normalized);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    s_[i].v = v;
    s_[i].q = q_normalized;
  }
  initConstraints();
  const bool feasible = isCurrentSolutionFeasible();
  return feasible;
}


bool OCP::setStateTrajectory(const Eigen::VectorXd& q0, 
                             const Eigen::VectorXd& v0, 
                             const Eigen::VectorXd& qN, 
                             const Eigen::VectorXd& vN) {
  assert(q0.size() == robots_[0].dimq());
  assert(v0.size() == robots_[0].dimv());
  assert(qN.size() == robots_[0].dimq());
  assert(vN.size() == robots_[0].dimv());
  Eigen::VectorXd q0_normalized = q0;
  robots_[0].normalizeConfiguration(q0_normalized);
  Eigen::VectorXd qN_normalized = qN;
  robots_[0].normalizeConfiguration(qN_normalized);
  const Eigen::VectorXd a = (vN-v0) / N_;
  Eigen::VectorXd dqN = Eigen::VectorXd::Zero(robots_[0].dimv());
  robots_[0].subtractConfiguration(qN, q0, dqN);
  const Eigen::VectorXd v = dqN / N_;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    s_[i].a = a;
    s_[i].v = v0 + i * a;
    robots_[0].integrateConfiguration(q0, v, (double)i, s_[i].q);
  }
  initConstraints();
  const bool feasible = isCurrentSolutionFeasible();
  return feasible;
}


void OCP::setContactStatus(const ContactStatus& contact_status) {
  contact_sequence_.setContactStatusUniformly(contact_status);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


ContactSequence& OCP::getContactSequenceHandle() {
  return contact_sequence_;
}


void OCP::setContactPoint(
    const std::vector<Eigen::Vector3d>& contact_points) {
  assert(contact_points.size() == robots_[0].max_point_contacts());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<robots_.size(); ++i) {
    robots_[i].setContactPoints(contact_points);
  }
}


void OCP::setContactPointByKinematics(const Eigen::VectorXd& q) {
  assert(q.size() == robots_[0].dimq());
  const int dimv = robots_[0].dimv();
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<robots_.size(); ++i) {
    robots_[i].updateKinematics(q, Eigen::VectorXd::Zero(dimv), 
                                Eigen::VectorXd::Zero(dimv));
    robots_[i].setContactPointsByCurrentKinematics();
  }
}


void OCP::clearLineSearchFilter() {
  filter_.clear();
}


double OCP::KKTError() {
  Eigen::VectorXd error = Eigen::VectorXd::Zero(N_+1);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i < N_) {
      error(i) = split_ocps_[i].squaredNormKKTResidual(dtau_);
    }
    else {
      error(N_) = terminal_ocp_.squaredNormKKTResidual();
    }
  }
  return std::sqrt(error.sum());
}


void OCP::computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i == 0) {
      split_ocps_[i].computeKKTResidual(robots_[robot_id], 
                                        contact_sequence_.contactStatus(i), 
                                        t+(i+1)*dtau_, dtau_, q, s_[i], s_[i+1]);
    }
    else if (i < N_) {
      split_ocps_[i].computeKKTResidual(robots_[robot_id], 
                                        contact_sequence_.contactStatus(i), 
                                        t+(i+1)*dtau_, dtau_, s_[i-1].q, s_[i], 
                                        s_[i+1]);
    }
    else {
      terminal_ocp_.computeKKTResidual(robots_[robot_id], t+T_, s_[N_]);
    }
  }
}


void OCP::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
    std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
    std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
    std::cout << "f[" << i << "] = ";
    for (int j=0; j<s_[i].f.size(); ++j) {
      std::cout << s_[i].f[j].transpose() << "; ";
    }
    std::cout << std::endl;
    std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    std::cout << "mu[" << i << "] = " << s_[i].mu_stack().transpose() << std::endl;
  }
  std::cout << "q[" << N_ << "] = " << s_[N_].q.transpose() << std::endl;
  std::cout << "v[" << N_ << "] = " << s_[N_].v.transpose() << std::endl;
}


bool OCP::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    const bool feasible = split_ocps_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time step " << i << std::endl;
      return false;
    }
  }
  return true;
}


void OCP::initConstraints() {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].initConstraints(robots_[robot_id], i, dtau_, s_[i]);
  }
}



// void OCP::linearizeOCP(const double t, const Eigen::VectorXd& q, 
//                        const Eigen::VectorXd& v) {
//   assert(q.size() == robots_[0].dimq());
//   assert(v.size() == robots_[0].dimv());
//   const int N_impulse = contact_sequence_.numImpulse();
//   const int N_lift = contact_sequence_.numLifts();
//   const int N = N_ + N_impulse + N_lift;
//   #pragma omp parallel for num_threads(num_proc_)
//   for (int i=0; i<=N; ++i) {
//     if (i == 0) {
//       const int robot_id = omp_get_thread_num();
//       split_ocps_[i].linearizeOCP(robots_[robot_id], 
//                                   contact_sequence_.contactStatus(i), 
//                                   t+i*dtau_, dtau_, q, s_[i], s_[i+1]);
//     }
//     else if (i < N_) {
//       const int robot_id = omp_get_thread_num();
//       split_ocps_[i].linearizeOCP(robots_[robot_id], 
//                                   contact_sequence_.contactStatus(i), 
//                                   t+i*dtau_, dtau_, s_[i-1].q, s_[i], s_[i+1]);
//     }
//     else if (i == N_) {
//       const int robot_id = omp_get_thread_num();
//       terminal_ocp_.linearizeOCP(robots_[robot_id], t+T_, s_[N_]);
//       terminal_ocp_.backwardRiccatiRecursion(riccati_[N_]);
//     }
//     else if (i < N_+N_impulse) {
//       const int robot_id = omp_get_thread_num();
//       split_impulse_ocps_[i].linearizeOCP(robots_[robot_id], 
//                                           contact_sequence_.impulseStatus(i-N_), 
//                                           t+i*dtau_, dtau_, s_[i-1].q, s_[i], s_[i+1]);
//     }
//     else {
//       const int robot_id = omp_get_thread_num();
//       split_lift_ocps_[i].linearizeOCP(robots_[robot_id], contact_sequence_.impulseStatus(i+1), 
//                                        t+i*dtau_, dtau_, s_[i-1].q, s_[i], s_[i+1]);
        
//     }
//   }
// }

} // namespace idocp
#include "idocp/ocp/riccati_direction_calculator.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

RiccatiDirectionCalculator::RiccatiDirectionCalculator(
    const double T, const int N, const int max_num_impulse, const int num_proc) 
  : T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc),
    max_primal_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)), 
    max_dual_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (num_proc <= 0) {
      throw std::out_of_range("invalid value: num_proc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


RiccatiDirectionCalculator::RiccatiDirectionCalculator()
  : T_(0),
    dtau_(0),
    N_(0),
    num_proc_(0),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


RiccatiDirectionCalculator::~RiccatiDirectionCalculator() {
}


void RiccatiDirectionCalculator::computeInitialStateDirection(
    const std::vector<Robot>& robots, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, Direction& d) {
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  robots[0].subtractConfiguration(q, s[0].q, d[0].dq());
  d[0].dv() = v - s[0].v;
}


void RiccatiDirectionCalculator::computeNewtonDirectionFromRiccatiFactorization(
    OCP& ocp, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, 
    const RiccatiFactorizer& factorizer, 
    const RiccatiFactorization& factorization, const Solution& s, Direction& d) {
  assert(robots.size() == num_proc_);
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  N_all_ = N_ + 1 + 2 * N_impulse + N_lift;
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all_; ++i) {
    if (i < N_) {
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d[i], 
                                                      exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          next_riccati_factorization(factorization, contact_sequence, i), d[i], 
          exist_state_constraint);
      ocp[i].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
                                             dtau(contact_sequence, i), s[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) = ocp[i].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = ocp[i].maxDualStepSize();
    }
    else if (i == N_) {
      SplitRiccatiFactorizer::computeCostateDirection(factorization[N_], d[N_], 
                                                      false);
      max_primal_step_sizes_.coeffRef(N_) = ocp.terminal.maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(N_) = ocp.terminal.maxDualStepSize();
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      ImpulseSplitRiccatiFactorizer::computeCostateDirection(
          factorization.impulse[impulse_index], d.impulse[impulse_index], 
          exist_state_constraint);
      ocp.impulse[impulse_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], s.impulse[impulse_index], 
          d.impulse[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) = ocp.impulse[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = ocp.impulse[impulse_index].maxDualStepSize();
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_after_impulse 
          = contact_sequence.timeStageAfterImpulse(impulse_index);
      const double dtau_aux 
          = time_stage_after_impulse * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.aux[impulse_index], d.aux[impulse_index], 
          exist_state_constraint);
      factorizer.aux[impulse_index].computeControlInputDirection(
          factorization[time_stage_after_impulse], d.aux[impulse_index], 
          exist_state_constraint);
      ocp.aux[impulse_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], dtau_aux, s.aux[impulse_index], 
          d.aux[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) = ocp.aux[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = ocp.aux[impulse_index].maxDualStepSize();
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_after_lift 
          = contact_sequence.timeStageAfterLift(lift_index);
      const double dtau_aux
          = time_stage_after_lift * dtau_ 
              - contact_sequence.liftTime(lift_index);
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.lift[lift_index], d.lift[lift_index], 
          exist_state_constraint);
      factorizer.lift[lift_index].computeControlInputDirection(
          factorization[time_stage_after_lift], d.lift[lift_index], 
          exist_state_constraint);
      ocp.lift[lift_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], dtau_aux, s.lift[lift_index], 
          d.lift[lift_index]);
      max_primal_step_sizes_.coeffRef(i) = ocp.lift[lift_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = ocp.lift[lift_index].maxDualStepSize();
    }
  }
}


double RiccatiDirectionCalculator::maxPrimalStepSize() const {
  return max_primal_step_sizes_.head(N_all_).minCoeff();
}


double RiccatiDirectionCalculator::maxDualStepSize() const {
  return max_dual_step_sizes_.head(N_all_).minCoeff();
}

} // namespace idocp 
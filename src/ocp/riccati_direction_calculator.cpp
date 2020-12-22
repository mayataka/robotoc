#include "idocp/ocp/riccati_direction_calculator.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

RiccatiDirectionCalculator::RiccatiDirectionCalculator(
    const int N, const int max_num_impulse, const int num_proc) 
  : N_(N),
    num_proc_(num_proc),
    max_primal_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)), 
    max_dual_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
  try {
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
  : N_(0),
    num_proc_(0),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


RiccatiDirectionCalculator::~RiccatiDirectionCalculator() {
}


void RiccatiDirectionCalculator::computeNewtonDirectionFromRiccatiFactorization(
    OCP& ocp, std::vector<Robot>& robots, const RiccatiFactorizer& factorizer, 
    const RiccatiFactorization& factorization, const Solution& s, Direction& d) {
  assert(robots.size() == num_proc_);
  const int N_impulse = ocp.discrete().numImpulseStages();
  const int N_lift = ocp.discrete().numLiftStages();
  N_all_ = N_ + 1 + 2 * N_impulse + N_lift;
  const bool exist_state_constraint = ocp.discrete().existStateConstraint();
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all_; ++i) {
    if (i < N_) {
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d[i], 
                                                      exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          next_riccati_factorization(ocp.discrete(), factorization, i), d[i], 
          exist_state_constraint);
      ocp[i].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
                                             ocp.discrete().dtau(i), s[i], d[i]);
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
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.impulse[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.impulse[impulse_index].maxDualStepSize();
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.aux[impulse_index], d.aux[impulse_index], 
          exist_state_constraint);
      factorizer.aux[impulse_index].computeControlInputDirection(
          factorization[ocp.discrete().timeStageAfterImpulse(impulse_index)], 
          d.aux[impulse_index], exist_state_constraint);
      ocp.aux[impulse_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], ocp.discrete().dtau_aux(impulse_index), 
          s.aux[impulse_index], d.aux[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.aux[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.aux[impulse_index].maxDualStepSize();
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.lift[lift_index], d.lift[lift_index], 
          exist_state_constraint);
      factorizer.lift[lift_index].computeControlInputDirection(
          factorization[ocp.discrete().timeStageAfterLift(lift_index)], 
          d.lift[lift_index], exist_state_constraint);
      ocp.lift[lift_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], ocp.discrete().dtau_lift(lift_index), 
          s.lift[lift_index], d.lift[lift_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.lift[lift_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.lift[lift_index].maxDualStepSize();
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
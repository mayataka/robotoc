#include "idocp/ocp/riccati_direction_calculator.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

RiccatiDirectionCalculator::RiccatiDirectionCalculator(
    const int N, const int max_num_impulse, const int nthreads) 
  : N_(N),
    nthreads_(nthreads),
    N_all_(N+1),
    max_primal_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)), 
    max_dual_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


RiccatiDirectionCalculator::RiccatiDirectionCalculator()
  : N_(0),
    nthreads_(0),
    N_all_(0),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


RiccatiDirectionCalculator::~RiccatiDirectionCalculator() {
}


void RiccatiDirectionCalculator::computeInitialStateDirection(
    const std::vector<Robot>& robots, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const KKTMatrix& kkt_matrix, const Solution& s, 
    Direction& d) {
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  if (robots[0].hasFloatingBase()) {
    robots[0].subtractConfiguration(q, s[0].q, d[0].dq());
    d[0].dq().template head<6>() 
        = - kkt_matrix[0].Fqq_prev_inv * d[0].dq().template head<6>();
    d[0].dv() = v - s[0].v;
  }
  else {
    d[0].dq() = q - s[0].q;
    d[0].dv() = v - s[0].v;
  }
}


void RiccatiDirectionCalculator::computeNewtonDirectionFromRiccatiFactorization(
    OCP& ocp, std::vector<Robot>& robots, 
    const RiccatiFactorization& factorization, 
    const Solution& s, Direction& d) {
  assert(robots.size() == nthreads_);
  const int N_impulse = ocp.discrete().numImpulseStages();
  const int N_lift = ocp.discrete().numLiftStages();
  N_all_ = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all_; ++i) {
    if (i < N_) {
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d[i]);
      ocp[i].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
                                             ocp.discrete().dtau(i), s[i], d[i]);
      if (ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i+1);
        SplitRiccatiFactorizer::computeLagrangeMultiplierDirection(
            factorization.constraint[impulse_index], d[i]);
      }
      max_primal_step_sizes_.coeffRef(i) = ocp[i].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = ocp[i].maxDualStepSize();
    }
    else if (i == N_) {
      SplitRiccatiFactorizer::computeCostateDirection(factorization[N_], d[N_]);
      max_primal_step_sizes_.coeffRef(N_) = ocp.terminal.maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(N_) = ocp.terminal.maxDualStepSize();
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      ImpulseSplitRiccatiFactorizer::computeCostateDirection(
          factorization.impulse[impulse_index], d.impulse[impulse_index]);
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
          factorization.aux[impulse_index], d.aux[impulse_index]);
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
          factorization.lift[lift_index], d.lift[lift_index]);
      ocp.lift[lift_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], ocp.discrete().dtau_lift(lift_index), 
          s.lift[lift_index], d.lift[lift_index]);
      const int time_stage_after_lift 
          = ocp.discrete().timeStageAfterLift(lift_index);
      if (ocp.discrete().isTimeStageBeforeImpulse(time_stage_after_lift)) {
        const int impulse_index
            = ocp.discrete().impulseIndexAfterTimeStage(time_stage_after_lift);
        SplitRiccatiFactorizer::computeLagrangeMultiplierDirection(
            factorization.constraint[impulse_index], d.lift[lift_index]);
      }
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
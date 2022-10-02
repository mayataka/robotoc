#include "robotoc/riccati/riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>

namespace robotoc {

RiccatiRecursion::RiccatiRecursion(const OCP& ocp, const int nthreads, 
                                   const double max_dts0)
  : nthreads_(nthreads),
    N_all_(ocp.N()+1),
    factorizer_(ocp.robot(), max_dts0),
    lqr_policy_(ocp.N()+3*ocp.reservedNumDiscreteEvents()+1, LQRPolicy(ocp.robot())),
    sto_policy_(ocp.N()+3*ocp.reservedNumDiscreteEvents()+1, STOPolicy(ocp.robot())),
    factorization_m_(ocp.robot()),
    max_primal_step_sizes_(
        Eigen::VectorXd::Zero(ocp.N()+1+3*ocp.reservedNumDiscreteEvents())), 
    max_dual_step_sizes_(
        Eigen::VectorXd::Zero(ocp.N()+1+3*ocp.reservedNumDiscreteEvents())) {
  if (nthreads <= 0) {
    throw std::out_of_range("[RiccatiRecursion] invalid argument: 'nthreads' must be positive!");
  }
}


RiccatiRecursion::RiccatiRecursion()
  : nthreads_(0),
    N_all_(0),
    factorizer_(),
    lqr_policy_(),
    sto_policy_(),
    factorization_m_(),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


void RiccatiRecursion::setRegularization(const double max_dts0) {
  assert(max_dts0 > 0);
  factorizer_.setRegularization(max_dts0);
}


void RiccatiRecursion::reserve(const OCP& ocp) {
  // const int reserved_num_discrete_events = ocp.reservedNumDiscreteEvents();
  // lqr_policy_.reserve(ocp.robot(), reserved_num_discrete_events);
  // while (sto_policy_.size() < 2*reserved_num_discrete_events+1) {
  //   sto_policy_.emplace_back(ocp.robot());
  // }
  // const int N = ocp.timeDiscretization().N();
  // const int max_N_all = N + 1 + 3*reserved_num_discrete_events;
  // if (max_N_all > max_primal_step_sizes_.size()) {
  //   max_primal_step_sizes_.resize(max_N_all);
  //   max_dual_step_sizes_.resize(max_N_all);
  // }
}


void RiccatiRecursion::backwardRiccatiRecursion(
    const TimeDiscretization& time_discretization, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual, RiccatiFactorization& factorization) {
  const int N = time_discretization.N_grids();
  factorization[N].P = kkt_matrix[N].Qxx;
  factorization[N].s = - kkt_residual[N].lx;
  for (int i=N-1; i>=0; --i) {
    const auto& grid = time_discretization.grid(i);
    if (grid.type == GridType::Impulse) {
      if (time_discretization.grid(i-1).sto || grid.sto) {
        factorizer_.backwardRiccatiRecursionPhaseTransition(
            factorization[i+1], factorization_m_, sto_policy_[i], grid.sto_next);
        factorizer_.backwardRiccatiRecursion(factorization_m_, kkt_matrix[i], 
                                             kkt_residual[i], factorization[i],
                                             grid.sto);
      }
      else {
        factorizer_.backwardRiccatiRecursion(factorization[i+1], kkt_matrix[i], 
                                             kkt_residual[i], factorization[i],
                                             grid.sto);
      }
    }
    else if (time_discretization.grid(i+1).type == GridType::Lift) {
      if (grid.sto || grid.sto_next) {
        factorizer_.backwardRiccatiRecursionPhaseTransition(
            factorization[i+1], factorization_m_, sto_policy_[i+1], grid.sto_next);
        factorizer_.backwardRiccatiRecursion(factorization_m_, kkt_matrix[i], 
                                             kkt_residual[i], factorization[i],
                                             lqr_policy_[i], grid.sto, grid.sto_next);
      }
      else {
        factorizer_.backwardRiccatiRecursion(factorization[i+1], kkt_matrix[i], 
                                             kkt_residual[i], factorization[i],
                                             lqr_policy_[i], grid.sto, grid.sto_next);
      }
    }
    else {
      factorizer_.backwardRiccatiRecursion(factorization[i+1], kkt_matrix[i], 
                                           kkt_residual[i], factorization[i], 
                                           lqr_policy_[i], grid.sto, grid.sto_next);
    }
  }
  const auto& grid = time_discretization.grid(0);
  if (grid.sto) {
    factorizer_.backwardRiccatiRecursionPhaseTransition(
        factorization[0], factorization_m_, sto_policy_[0], grid.sto_next);
  }
}


void RiccatiRecursion::forwardRiccatiRecursion(
    const TimeDiscretization& time_discretization, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, const RiccatiFactorization& factorization,
    Direction& d) const {
  const int N = time_discretization.N_grids();
  d[0].dts = 0.0;
  d[0].dts_next = 0.0;
  if (time_discretization.grid(0).sto) {
    constexpr bool sto_prev = false;
    ::robotoc::computeSwitchingTimeDirection(sto_policy_[0], d[0], sto_prev);
  }
  for (int i=0; i<N; ++i) {
    const auto& grid = time_discretization.grid(i);
    if (grid.type == GridType::Impulse) {
      d[i].dts = d[i-1].dts_next;
      d[i].dts_next = 0.0;
      ::robotoc::forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i], d[i], d[i+1]);
      if (grid.sto_next) {
        d[i+1].dts = d[i-1].dts_next;
        d[i+1].dts_next = 0.0;
        ::robotoc::computeSwitchingTimeDirection(sto_policy_[i], d[i+1], grid.sto);
        d[i].dts = d[i+1].dts;
        d[i].dts_next = d[i+1].dts_next;
      }
      ::robotoc::computeCostateDirection(factorization[i], d[i], grid.sto);
    }
    else if (grid.type == GridType::Lift) {
      d[i].dts = d[i-1].dts_next;
      d[i].dts_next = 0.0;
      if (grid.sto_next) {
        ::robotoc::computeSwitchingTimeDirection(sto_policy_[i], d[i], grid.sto);
      }
      ::robotoc::forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i], lqr_policy_[i], 
                                         d[i], d[i+1], grid.sto, grid.sto_next);
      ::robotoc::computeCostateDirection(factorization[i], d[i], grid.sto, grid.sto_next);
    }
    else {
      ::robotoc::forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  lqr_policy_[i], 
                                         d[i], d[i+1], grid.sto, grid.sto_next);
      ::robotoc::computeCostateDirection(factorization[i], d[i], grid.sto, grid.sto_next);
    }
    if (grid.switching_constraint) {
      ::robotoc::computeLagrangeMultiplierDirection(factorization[i], d[i], grid.sto, grid.sto_next);
    }
  }
  constexpr bool sto = false; 
  constexpr bool sto_next = false; 
  ::robotoc::computeCostateDirection(factorization[N], d[N], sto, sto_next);
}


const aligned_vector<LQRPolicy>& RiccatiRecursion::getLQRPolicy() const {
  return lqr_policy_;
}

} // namespace robotoc
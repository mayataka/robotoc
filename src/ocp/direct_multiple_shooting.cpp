#include "robotoc/ocp/direct_multiple_shooting.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc{

DirectMultipleShooting::DirectMultipleShooting(const int nthreads) 
  : nthreads_(nthreads) {
  try {
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


DirectMultipleShooting::DirectMultipleShooting()
  : nthreads_(0) {
}


DirectMultipleShooting::~DirectMultipleShooting() {
}


bool DirectMultipleShooting::isFeasible(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Solution& s) const {
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  std::vector<bool> is_feasible(N_all, true);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      const int contact_phase = ocp.discrete().contactPhase(i);
      is_feasible[i] = ocp[i].isFeasible(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(contact_phase), s[i]);
    }
    else if (i == N) {
      is_feasible[i] = ocp.terminal.isFeasible(robots[omp_get_thread_num()], s[N]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      is_feasible[i] = ocp.impulse[impulse_index].isFeasible(
          robots[omp_get_thread_num()], 
          contact_sequence->impulseStatus(impulse_index), 
          s.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      const int contact_phase = ocp.discrete().contactPhaseAfterImpulse(impulse_index);
      is_feasible[i] = ocp.aux[impulse_index].isFeasible(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(contact_phase), s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int contact_phase = ocp.discrete().contactPhaseAfterLift(lift_index);
      is_feasible[i] = ocp.lift[lift_index].isFeasible(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(lift_index), s.lift[lift_index]);
    }
  }
  for (int i=0; i<N_all; ++i) {
    if (!is_feasible[i]) {
      return false;
    }
  }
  return true;
}


void DirectMultipleShooting::initConstraints(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Solution& s) const {
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      ocp[i].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(ocp.discrete().contactPhase(i)),
          i, s[i]);
    }
    else if (i == N) {
      ocp.terminal.initConstraints(robots[omp_get_thread_num()], N, s[N]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      ocp.impulse[impulse_index].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence->impulseStatus(impulse_index), 
          s.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      ocp.aux[impulse_index].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(
              ocp.discrete().contactPhaseAfterImpulse(impulse_index)), 
          0, s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      ocp.lift[lift_index].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(
              ocp.discrete().contactPhaseAfterLift(lift_index)), 
          0, s.lift[lift_index]);
    }
  }
}


void DirectMultipleShooting::computeKKTResidual(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) const {
  runParallel<internal::ComputeKKTResidual>(ocp, robots, contact_sequence, q, v,  
                                            s, kkt_matrix, kkt_residual);
}


void DirectMultipleShooting::computeKKTSystem(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) const {
  runParallel<internal::ComputeKKTSystem>(ocp, robots, contact_sequence, q, v, 
                                          s, kkt_matrix, kkt_residual);
}


double DirectMultipleShooting::KKTError(const OCP& ocp, 
                                        const KKTResidual& kkt_residual) {
  double kkt_error = 0;
  for (int i=0; i<=ocp.discrete().N(); ++i) {
    kkt_error += kkt_residual[i].kkt_error;
  }
  for (int i=0; i<ocp.discrete().N_impulse(); ++i) {
    kkt_error += kkt_residual.impulse[i].kkt_error;
    kkt_error += kkt_residual.aux[i].kkt_error;
  }
  for (int i=0; i<ocp.discrete().N_lift(); ++i) {
    kkt_error += kkt_residual.lift[i].kkt_error;
  }
  return kkt_error;
}


double DirectMultipleShooting::totalCost(const OCP& ocp) {
  double total_cost = 0;
  for (int i=0; i<ocp.discrete().N(); ++i) {
    total_cost += ocp[i].stageCost();
  }
  total_cost += ocp.terminal.terminalCost();
  for (int i=0; i<ocp.discrete().N_impulse(); ++i) {
    total_cost += ocp.impulse[i].stageCost();
    total_cost += ocp.aux[i].stageCost();
  }
  for (int i=0; i<ocp.discrete().N_lift(); ++i) {
    total_cost += ocp.lift[i].stageCost();
  }
  return total_cost;
}


void DirectMultipleShooting::computeInitialStateDirection(
    const OCP& ocp, const aligned_vector<Robot>& robots, 
    const Eigen::VectorXd& q0, const Eigen::VectorXd& v0, 
    const Solution& s, Direction& d) {
  ocp[0].computeInitialStateDirection(robots[0], q0, v0, s[0], d[0]);
}


void DirectMultipleShooting::integrateSolution(
    OCP& ocp, const aligned_vector<Robot>& robots, 
    const double primal_step_size, const double dual_step_size, 
    const KKTMatrix& kkt_matrix, Direction& d, Solution& s) const {
  assert(robots.size() == nthreads_);
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i);
        ocp[i].expandDual(ocp.discrete().dt(i), d.impulse[impulse_index], d[i], 
                          dts_stage(ocp, d, i));
      }
      else if (ocp.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp.discrete().liftIndexAfterTimeStage(i);
        ocp[i].expandDual(ocp.discrete().dt(i), d.lift[lift_index], d[i], 
                          dts_stage(ocp, d, i));
      }
      else if (ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i+1);
        ocp[i].expandDual(ocp.discrete().dt(i), d[i+1], 
                          kkt_matrix.switching[impulse_index], d[i], 
                          dts_stage(ocp, d, i));
      }
      else {
        ocp[i].expandDual(ocp.discrete().dt(i), d[i+1], d[i], 
                          dts_stage(ocp, d, i));
      }
      ocp[i].updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                          d[i], s[i]);
      ocp[i].updateDual(dual_step_size);
    }
    else if (i == N) {
      ocp.terminal.expandDual(d[N]);
      ocp.terminal.updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                                d[N], s[N]);
      ocp.terminal.updateDual(dual_step_size);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      ocp.impulse[impulse_index].expandDual(d.aux[impulse_index], 
                                            d.impulse[impulse_index]);
      ocp.impulse[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                              primal_step_size, 
                                              d.impulse[impulse_index], 
                                              s.impulse[impulse_index]);
      ocp.impulse[impulse_index].updateDual(dual_step_size);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      ocp.aux[impulse_index].expandDual(
          ocp.discrete().dt_aux(impulse_index), 
          d[ocp.discrete().timeStageAfterImpulse(impulse_index)], 
          d.aux[impulse_index], dts_aux(ocp, d, impulse_index));
      ocp.aux[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                          primal_step_size, 
                                          d.aux[impulse_index], 
                                          s.aux[impulse_index]);
      ocp.aux[impulse_index].updateDual(dual_step_size);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      ocp.lift[lift_index].expandDual(
          ocp.discrete().dt_lift(lift_index), 
          d[ocp.discrete().timeStageAfterLift(lift_index)], 
          d.lift[lift_index], dts_lift(ocp, d, lift_index));
      ocp.lift[lift_index].updatePrimal(robots[omp_get_thread_num()], 
                                        primal_step_size, d.lift[lift_index], 
                                        s.lift[lift_index]);
      ocp.lift[lift_index].updateDual(dual_step_size);
    }
  }
}

} // namespace robotoc
#include "robotoc/ocp/direct_multiple_shooting.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc{

DirectMultipleShooting::DirectMultipleShooting(const int nthreads) 
  : ocp_data_(),
    intermediate_stage_(),
    impact_stage_(),
    terminal_stage_(),
    nthreads_(nthreads) {
  if (nthreads <= 0) {
    throw std::out_of_range("[DirectMultipleShooting] invalid argument: 'nthreads' must be positive!");
  }
}


DirectMultipleShooting::DirectMultipleShooting(
    const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints, 
    const std::shared_ptr<ContactSequence>& contact_sequence, const int nthreads) 
  : ocp_data_(),
    intermediate_stage_(cost, constraints, contact_sequence),
    impact_stage_(cost, constraints, contact_sequence),
    terminal_stage_(cost, constraints, contact_sequence),
    nthreads_(nthreads) {
  if (nthreads <= 0) {
    throw std::out_of_range("[DirectMultipleShooting] invalid argument: 'nthreads' must be positive!");
  }
}


DirectMultipleShooting(const OCPDef& ocp, const int nthreads)
  : DirectMultipleShooting(ocp.cost, ocp.constraints, ocp.contact_sequence, nthreads) {
  ocp_data_.resize(ocp.N+1);
  for (int i=0; i<=ocp.N; ++i) {
    ocp_data_[i] = intermediate_stage_.createData(ocp.robot);
  }
}


DirectMultipleShooting::DirectMultipleShooting()
  : ocp_data_(),
    intermediate_stage_(),
    impact_stage_(),
    terminal_stage_(),
    nthreads_(0) {
}


DirectMultipleShooting::~DirectMultipleShooting() {
}


bool DirectMultipleShooting::isFeasible(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Solution& s) const {
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  std::vector<bool> is_feasible(N_all, true);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      const int contact_phase = ocp.timeDiscretization().contactPhase(i);
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
      const int contact_phase = ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index);
      is_feasible[i] = ocp.aux[impulse_index].isFeasible(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(contact_phase), s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int contact_phase = ocp.timeDiscretization().contactPhaseAfterLift(lift_index);
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
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      ocp[i].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(ocp.timeDiscretization().contactPhase(i)),
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
              ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index)), 
          0, s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      ocp.lift[lift_index].initConstraints(
          robots[omp_get_thread_num()], 
          contact_sequence->contactStatus(
              ocp.timeDiscretization().contactPhaseAfterLift(lift_index)), 
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
  for (int i=0; i<=ocp.timeDiscretization().N(); ++i) {
    kkt_error += kkt_residual[i].kkt_error;
  }
  for (int i=0; i<ocp.timeDiscretization().N_impulse(); ++i) {
    kkt_error += kkt_residual.impulse[i].kkt_error;
    kkt_error += kkt_residual.aux[i].kkt_error;
  }
  for (int i=0; i<ocp.timeDiscretization().N_lift(); ++i) {
    kkt_error += kkt_residual.lift[i].kkt_error;
  }
  return kkt_error;
}


double DirectMultipleShooting::totalCost(const OCP& ocp, 
                                         const bool include_cost_barrier) {
  double total_cost = 0;
  for (int i=0; i<ocp.timeDiscretization().N(); ++i) {
    total_cost += ocp[i].stageCost(include_cost_barrier);
  }
  total_cost += ocp.terminal.terminalCost(include_cost_barrier);
  for (int i=0; i<ocp.timeDiscretization().N_impulse(); ++i) {
    total_cost += ocp.impulse[i].stageCost(include_cost_barrier);
    total_cost += ocp.aux[i].stageCost(include_cost_barrier);
  }
  for (int i=0; i<ocp.timeDiscretization().N_lift(); ++i) {
    total_cost += ocp.lift[i].stageCost(include_cost_barrier);
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
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp.timeDiscretization().impulseIndexAfterTimeStage(i);
        ocp[i].expandDual(ocp.timeDiscretization().gridInfo(i), 
                          d.impulse[impulse_index], d[i], dts_stage(ocp, d, i));
      }
      else if (ocp.timeDiscretization().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp.timeDiscretization().liftIndexAfterTimeStage(i);
        ocp[i].expandDual(ocp.timeDiscretization().gridInfo(i), 
                          d.lift[lift_index], d[i], dts_stage(ocp, d, i));
      }
      else if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp.timeDiscretization().impulseIndexAfterTimeStage(i+1);
        ocp[i].expandDual(ocp.timeDiscretization().gridInfo(i), d[i+1], 
                          kkt_matrix.switching[impulse_index], d[i], 
                          dts_stage(ocp, d, i));
      }
      else {
        ocp[i].expandDual(ocp.timeDiscretization().gridInfo(i), d[i+1], d[i], 
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
          ocp.timeDiscretization().gridInfoAux(impulse_index), 
          d[ocp.timeDiscretization().timeStageAfterImpulse(impulse_index)], 
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
          ocp.timeDiscretization().gridInfoLift(lift_index), 
          d[ocp.timeDiscretization().timeStageAfterLift(lift_index)], 
          d.lift[lift_index], dts_lift(ocp, d, lift_index));
      ocp.lift[lift_index].updatePrimal(robots[omp_get_thread_num()], 
                                        primal_step_size, d.lift[lift_index], 
                                        s.lift[lift_index]);
      ocp.lift[lift_index].updateDual(dual_step_size);
    }
  }
}


bool DirectMultipleShooting::isFeasible(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Solution& s) {
  const int N = time_discretization.N();
  std::vector<bool> is_feasible(N, true);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization.gridInfo(i);
    if (grid.type == GridType::Intermediate) {
      is_feasible[i] = intermediate_stage_.isFeasible(robots[omp_get_thread_num()], 
                                                      grid, s[i], ocp_data_[i]);
    }
    else if (grid.type == GridType::Impulse) {
      is_feasible[i] = impact_stage_.isFeasible(robots[omp_get_thread_num()], 
                                                grid, s[i], ocp_data_[i]);
    }
    else {
      is_feasible[i] = terminal_stage_.isFeasible(robots[omp_get_thread_num()], 
                                                  grid, s[i], ocp_data_[i]);
    }
  }
  for (const auto e : is_feasible) {
    if (!e) return false;
  }
  return true;
}


void DirectMultipleShooting::initConstraints(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Solution& s) {
  const int N = time_discretization.N();
  std::vector<bool> is_feasible(N, true);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization.gridInfo(i);
    if (grid.type == GridType::Intermediate) {
      intermediate_stage_.initConstraints(robots[omp_get_thread_num()], 
                                          grid, s[i], ocp_data_[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.initConstraints(robots[omp_get_thread_num()], 
                                    grid, s[i], ocp_data_[i]);
    }
    else {
      terminal_stage_.initConstraints(robots[omp_get_thread_num()], 
                                      grid, s[i], ocp_data_[i]);
    }
  }
}


void DirectMultipleShooting::evalOCP(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTResidual& kkt_residual) {
  const int N = time_discretization.N();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization.gridInfo(i);
    if (grid.type == GridType::Intermediate) {
      intermediate_stage_.evalOCP(robots[omp_get_thread_num()], grid, s[i], s[i+1], 
                                  ocp_data_[i], kkt_residual[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.evalOCP(robots[omp_get_thread_num()], grid, s[i], s[i+1], 
                            ocp_data_[i], kkt_residual[i]);
    }
    else {
      terminal_stage_.evalOCP(robots[omp_get_thread_num()], grid, s[i],  
                              ocp_data_[i], kkt_residual[i]);
    }
  }
}


void DirectMultipleShooting::evalKKT(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const int N = time_discretization.N();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization.gridInfo(i);
    if (i == 0) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], grid, q, s[i], s[i+1], 
                                  ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else if (grid.type == GridType::Intermediate) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], grid, s[i-1].q, s[i], s[i+1],
                                  ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.evalKKT(robots[omp_get_thread_num()], grid, s[i-1].q, s[i], s[i+1],
                            ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else {
      terminal_stage_.evalKKT(robots[omp_get_thread_num()], grid, s[i-1].q, s[i], 
                              ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
  }
}


void DirectMultipleShooting::computeInitialStateDirection(
    const Robot& robot,  const Eigen::VectorXd& q0, const Eigen::VectorXd& v0, 
    const Solution& s, Direction& d) const {
  ::robotoc::computeInitialStateDirection(robot, q0, v0, s[0], ocp_data_[0], d[0]);
}


PerformanceIndex DirectMultipleShooting::getEval() const {
  PerformanceIndex performance_index;
  for (const auto& e : ocp_data_) {
    performance_index += e.performance_index;
  }
  return performance_index;
}


void DirectMultipleShooting::integrateSolution(
    const aligned_vector<Robot>& robots, 
    const TimeDiscretization& time_discretization, 
    const double primal_step_size, const double dual_step_size, 
    const KKTMatrix& kkt_matrix, Direction& d, Solution& s) {
  const int N = time_discretization.N();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization.gridInfo(i);
    if (grid.type == GridType::Intermediate) {
      intermediate_stage_.expandDual(grid, ocp_data_[i], d[i+1], d[i]);
      intermediate_stage_.updatePrimal(robots[omp_get_thread_num()], 
                                       primal_step_size, d[i], s[i], ocp_data_[i]);
      intermediate_stage_.updateDual(dual_step_size, ocp_data_[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.expandDual(grid, ocp_data_[i], d[i+1], d[i]);
      impact_stage_.updatePrimal(robots[omp_get_thread_num()], 
                                 primal_step_size, d[i], s[i], ocp_data_[i]);
      impact_stage_.updateDual(dual_step_size, ocp_data_[i]);
    }
    else {
      terminal_stage_.expandDual(grid, ocp_data_[i], d[i]);
      terminal_stage_.updatePrimal(robots[omp_get_thread_num()], 
                                   primal_step_size, d[i], s[i], ocp_data_[i]);
      terminal_stage_.updateDual(dual_step_size, ocp_data_[i]);
    }
  }
}

} // namespace robotoc
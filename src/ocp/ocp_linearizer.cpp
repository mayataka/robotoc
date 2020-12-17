#include "idocp/ocp/ocp_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp{

OCPLinearizer::OCPLinearizer(const double T, const int N, 
                             const int max_num_impulse, 
                             const int num_proc) 
  : N_(N),
    num_proc_(num_proc),
    ocp_discretizer_(T, N, max_num_impulse),
    kkt_error_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
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


OCPLinearizer::OCPLinearizer()
  : N_(0),
    num_proc_(0),
    ocp_discretizer_(),
    kkt_error_() {
}


OCPLinearizer::~OCPLinearizer() {
}


void OCPLinearizer::initConstraints(OCP& ocp, std::vector<Robot>& robots, 
                                    const ContactSequence& contact_sequence, 
                                    const double t, const Solution& s) {
  ocp_discretizer_.discretizeOCP(contact_sequence, t);
  const int N_impulse = ocp_discretizer_.numImpulseStages();
  const int N_lift = ocp_discretizer_.numLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      ocp[i].initConstraints(robots[omp_get_thread_num()], i, 
                             ocp_discretizer_.dtau(i), s[i]);
    }
    else if (i == N_) {
      ocp.terminal.initConstraints(robots[omp_get_thread_num()], N_, 
                                   ocp_discretizer_.dtau(N_), s[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      ocp.impulse[impulse_index].initConstraints(robots[omp_get_thread_num()], 
                                                 s.impulse[impulse_index]);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      ocp.aux[impulse_index].initConstraints(
          robots[omp_get_thread_num()], 0, 
          ocp_discretizer_.dtau_aux(impulse_index), s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      ocp.lift[lift_index].initConstraints(
          robots[omp_get_thread_num()], 0, 
          ocp_discretizer_.dtau_lift(lift_index), s.lift[lift_index]);
    }
  }
}


void OCPLinearizer::linearizeOCP(OCP& ocp, std::vector<Robot>& robots, 
                                 const ContactSequence& contact_sequence, 
                                 const double t, const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, const Solution& s, 
                                 KKTMatrix& kkt_matrix, 
                                 KKTResidual& kkt_residual) {
  ocp_discretizer_.discretizeOCP(contact_sequence, t);
  runParallel<internal::LinearizeOCP>(ocp, robots, contact_sequence, t, q, v, s, 
                                      kkt_matrix, kkt_residual);
}


void OCPLinearizer::computeKKTResidual(OCP& ocp, std::vector<Robot>& robots, 
                                       const ContactSequence& contact_sequence, 
                                       const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const Solution& s, KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual) {
  ocp_discretizer_.discretizeOCP(contact_sequence, t);
  runParallel<internal::ComputeKKTResidual>(ocp, robots, contact_sequence, t, 
                                            q, v, s, kkt_matrix, kkt_residual);
}


double OCPLinearizer::KKTError(const OCP& ocp, 
                               const KKTResidual& kkt_residual) {
  const int N_impulse = ocp_discretizer_.numImpulseStages();
  const int N_lift = ocp_discretizer_.numLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      kkt_error_.coeffRef(i) 
          = ocp[i].squaredNormKKTResidual(kkt_residual[i], 
                                          ocp_discretizer_.dtau(i));
    }
    else if (i == N_) {
      kkt_error_.coeffRef(N_) 
          = ocp.terminal.squaredNormKKTResidual(kkt_residual[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = ocp_discretizer_.timeStageBeforeImpulse(impulse_index);
      const bool is_state_constraint_valid = (time_stage_before_impulse > 0);
      kkt_error_.coeffRef(i) 
          = ocp.impulse[impulse_index].squaredNormKKTResidual(
              kkt_residual.impulse[impulse_index], is_state_constraint_valid);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      kkt_error_.coeffRef(i) 
          = ocp.aux[impulse_index].squaredNormKKTResidual(
                kkt_residual.aux[impulse_index], 
                ocp_discretizer_.dtau_aux(impulse_index));
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      kkt_error_.coeffRef(i) 
          = ocp.lift[lift_index].squaredNormKKTResidual(
              kkt_residual.lift[lift_index], 
              ocp_discretizer_.dtau_lift(lift_index));
    }
  }
  return std::sqrt(kkt_error_.head(N_all).sum());
}


void OCPLinearizer::integrateSolution(OCP& ocp, const std::vector<Robot>& robots, 
                                      const KKTMatrix& kkt_matrix, 
                                      const KKTResidual& kkt_residual, 
                                      const double primal_step_size, 
                                      const double dual_step_size, 
                                      Direction& d, Solution& s) const {
  assert(robots.size() == num_proc_);
  const int N_impulse = ocp_discretizer_.numImpulseStages();
  const int N_lift = ocp_discretizer_.numLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      if (ocp_discretizer_.isTimeStageBeforeImpulse(i)) {
        ocp[i].computeCondensedDualDirection(
            robots[omp_get_thread_num()], ocp_discretizer_.dtau(i), 
            kkt_matrix[i], kkt_residual[i], 
            d.impulse[ocp_discretizer_.impulseIndex(i)], d[i]);
      }
      else if (ocp_discretizer_.isTimeStageBeforeLift(i)) {
        ocp[i].computeCondensedDualDirection(
            robots[omp_get_thread_num()], ocp_discretizer_.dtau(i), 
            kkt_matrix[i], kkt_residual[i], 
            d.lift[ocp_discretizer_.liftIndex(i)], d[i]);
      }
      else {
        ocp[i].computeCondensedDualDirection(
            robots[omp_get_thread_num()], ocp_discretizer_.dtau(i), 
            kkt_matrix[i], kkt_residual[i], d[i+1], d[i]);
      }
      ocp[i].updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                          d[i], s[i]);
      ocp[i].updateDual(dual_step_size);
    }
    else if (i == N_) {
      ocp.terminal.updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                                d[N_], s[N_]);
      ocp.terminal.updateDual(dual_step_size);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const bool is_state_constraint_valid 
          = (ocp_discretizer_.timeStageBeforeImpulse(impulse_index) > 0);
      ocp.impulse[impulse_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], kkt_matrix.impulse[impulse_index], 
          kkt_residual.impulse[impulse_index], d.aux[impulse_index], 
          d.impulse[impulse_index]);
      ocp.impulse[impulse_index].updatePrimal(
          robots[omp_get_thread_num()], primal_step_size, 
          d.impulse[impulse_index], s.impulse[impulse_index], 
          is_state_constraint_valid);
      ocp.impulse[impulse_index].updateDual(dual_step_size);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      ocp.aux[impulse_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], ocp_discretizer_.dtau_aux(impulse_index),
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index],
          d[ocp_discretizer_.timeStageAfterImpulse(impulse_index)], 
          d.aux[impulse_index]);
      ocp.aux[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                          primal_step_size, 
                                          d.aux[impulse_index], 
                                          s.aux[impulse_index]);
      ocp.aux[impulse_index].updateDual(dual_step_size);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      ocp.lift[lift_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], ocp_discretizer_.dtau_lift(lift_index), 
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
          d[ocp_discretizer_.timeStageAfterLift(lift_index)], 
          d.lift[lift_index]);
      ocp.lift[lift_index].updatePrimal(robots[omp_get_thread_num()], 
                                        primal_step_size, 
                                        d.lift[lift_index], s.lift[lift_index]);
      ocp.lift[lift_index].updateDual(dual_step_size);
    }
  }
}

} // namespace idocp
#include "idocp/ocp/ocp_linearizer.hpp"


// void OCPLinearizer::computeDirections(
//     std::vector<Robot>& robots, const ContactSequence& contact_sequence, 
//     const HybridRiccatiFactorizer& factorizer, 
//     HybridRiccatiFactorization& factorization, 
//     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization, 
//     const HybridSolution& s, HybridDirection& d) {
//   assert(robots.size() == num_proc_);
//   const int N_impulse = contact_sequence.totalNumImpulseStages();
//   const int N_lift = contact_sequence.totalNumLiftStages();
//   const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
//   const bool exist_state_constraint = contact_sequence.existImpulseStage();
//   #pragma omp parallel for num_threads(num_proc_)
//   for (int i=0; i<N_all; ++i) {
//     if (i == 0) {
//       computePrimalDirectionInitial(factorizer[0], factorization[0], d[0], 
//                                     exist_state_constraint);
//       split_ocps_[0].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
//                                                      dtau(contact_sequence, 0), 
//                                                      s[0], d[0]);
//       max_primal_step_sizes_.coeffRef(0) = split_ocps_[0].maxPrimalStepSize();
//       max_dual_step_sizes_.coeffRef(0) = split_ocps_[0].maxDualStepSize();
//     }
//     else if (i < N_) {
//       aggregateLagrangeMultiplierDirection(contact_sequence, 
//                                            constraint_factorization, d.impulse, 
//                                            i, factorization[i]);
//       computePrimalDirection(factorizer[i], factorization[i], d[0].dx(), d[i], 
//                              exist_state_constraint);
//       split_ocps_[i].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
//                                                      dtau(contact_sequence, i), 
//                                                      s[i], d[i]);
//       max_primal_step_sizes_.coeffRef(i) = split_ocps_[i].maxPrimalStepSize();
//       max_dual_step_sizes_.coeffRef(i) = split_ocps_[i].maxDualStepSize();
//     }
//     else if (i == N_) {
//       computePrimalDirectionTerminal(factorization[N_], d[0].dx(), d[N_]);
//       max_primal_step_sizes_.coeffRef(N_) = 1.0;
//       max_dual_step_sizes_.coeffRef(N_) = 1.0;
//     }
//     else if (i < N_ + 1 + N_impulse) {
//       const int impulse_index  = i - (N_+1);
//       aggregateLagrangeMultiplierDirectionImpulse(contact_sequence, 
//                                                   constraint_factorization, 
//                                                   d.impulse, impulse_index, 
//                                                   factorization.impulse[impulse_index]);
//       computePrimalDirectionImpulse(factorization.impulse[impulse_index], 
//                                     d[0].dx(), d.impulse[impulse_index]);
//       split_ocps_.impulse[impulse_index].computeCondensedPrimalDirection(
//           robots[omp_get_thread_num()], s.impulse[impulse_index], d.impulse[impulse_index]);
//       max_primal_step_sizes_.coeffRef(i) = split_ocps_.impulse[impulse_index].maxPrimalStepSize();
//       max_dual_step_sizes_.coeffRef(i) = split_ocps_.impulse[impulse_index].maxDualStepSize();
//     }
//     else if (i < N_ + 1 + 2*N_impulse) {
//       const int impulse_index  = i - (N_+1+N_impulse);
//       const int time_stage_before_impulse 
//           = contact_sequence.timeStageBeforeImpulse(impulse_index);
//       const double dtau_aux 
//           = (time_stage_before_impulse+1) * dtau_ 
//               - contact_sequence.impulseTime(impulse_index);
//       aggregateLagrangeMultiplierDirectionAux(contact_sequence, 
//                                               constraint_factorization, d.impulse, 
//                                               impulse_index, 
//                                               factorization[impulse_index]);
//       computePrimalDirection(factorizer.aux[impulse_index], 
//                              factorization.aux[impulse_index], d[0].dx(), 
//                              d.aux[impulse_index], exist_state_constraint);
//       split_ocps_.aux[impulse_index].computeCondensedPrimalDirection(
//           robots[omp_get_thread_num()], dtau_aux, s.aux[impulse_index], d.aux[impulse_index]);
//       max_primal_step_sizes_.coeffRef(i) = split_ocps_.aux[impulse_index].maxPrimalStepSize();
//       max_dual_step_sizes_.coeffRef(i) = split_ocps_.aux[impulse_index].maxDualStepSize();
//     }
//     else {
//       const int lift_index = i - (N_+1+2*N_impulse);
//       const int time_stage_before_lift 
//           = contact_sequence.timeStageBeforeLift(lift_index);
//       const double dtau_aux
//           = (time_stage_before_lift+1) * dtau_ 
//               - contact_sequence.liftTime(lift_index);
//       aggregateLagrangeMultiplierDirectionLift(contact_sequence, 
//                                                constraint_factorization, d.impulse, 
//                                                lift_index, 
//                                                factorization[lift_index]);
//       computePrimalDirection(factorizer.lift[lift_index], 
//                              factorization.lift[lift_index], d[0].dx(), 
//                              d.lift[lift_index], exist_state_constraint);
//       split_ocps_.lift[lift_index].computeCondensedPrimalDirection(
//           robots[omp_get_thread_num()], dtau_aux, s.lift[lift_index], d.lift[lift_index]);
//       max_primal_step_sizes_.coeffRef(i) = split_ocps_.lift[lift_index].maxPrimalStepSize();
//       max_dual_step_sizes_.coeffRef(i) = split_ocps_.lift[lift_index].maxDualStepSize();
//     }
//   }
// }


// double OCPLinearizer::maxPrimalStepSize(
//     const ContactSequence& contact_sequence) const {
//   const int N_impulse = contact_sequence.totalNumImpulseStages();
//   const int N_lift = contact_sequence.totalNumLiftStages();
//   const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
//   return max_primal_step_sizes_.head(N_all).minCoeff();
// }


// double OCPLinearizer::maxDualStepSize(
//     const ContactSequence& contact_sequence) const {
//   const int N_impulse = contact_sequence.totalNumImpulseStages();
//   const int N_lift = contact_sequence.totalNumLiftStages();
//   const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
//   return max_dual_step_sizes_.head(N_all).minCoeff();
// }


// void OCPLinearizer::updateSolution(const std::vector<Robot>& robots, 
//                                    const HybridKKTMatrix& kkt_matrix,
//                                    const HybridKKTResidual& kkt_residual,
//                                    const double primal_step_size, 
//                                    const double dual_step_size, 
//                                    HybridDirection& d, 
//                                    HybridSolution& s) {
//   // assert(robots.size() == num_proc_);
//   // assert(q.size() == robots[0].dimq());
//   // assert(v.size() == robots[0].dimv());
//   // const int N_impulse = contact_sequence.totalNumImpulseStages();
//   // const int N_lift = contact_sequence.totalNumLiftStages();
//   // const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
//   // #pragma omp parallel for num_threads(num_proc_)
//   // for (int i=0; i<N_all; ++i) {
//   //   if (i == 0) {
//   //     if (contact_sequence.existImpulseStage(0)) {
//   //       constexpr int impulse_index = 0;
//   //       const double dtau_impulse = contact_sequence.impulseTime(impulse_index);
//   //       assert(dtau_impulse > 0);
//   //       assert(dtau_impulse < dtau_);
//   //       const int robot_id = omp_get_thread_num();
//   //       split_ocps_[0].computeCondensedDualDirection(robots[robot_id], 
//   //                                                    dtau_impulse, 
//   //                                                    kkt_matrix[0],
//   //                                                    kkt_residual[0],
//   //                                                    d.impulse[impulse_index],
//   //                                                    d[0]);
//   //       split_ocps_[0].updatePrimal(robots[robot_id], primal_step_size, 
//   //                                   dtau_impulse, d[0], s[0]);
//   //       split_ocps_[0].updateDual(dual_step_size);
//   //     }
//   //     else if (contact_sequence.existLiftStage(0)) {
//   //       constexpr int lift_index = 0;
//   //       const double dtau_lift = contact_sequence.liftTime(lift_index);
//   //       assert(dtau_lift > 0);
//   //       assert(dtau_lift < dtau_);
//   //       split_ocps_[0].computeKKTResidual(robots[omp_get_thread_num()], 
//   //                                         contact_sequence.contactStatus(0), 
//   //                                         t, dtau_lift, q, s[0], s.lift[lift_index],
//   //                                         kkt_matrix[0], kkt_residual[0]);
//   //     }
//   //     else {
//   //       const int robot_id = omp_get_thread_num();
//   //       split_ocps_[0].computeCondensedDualDirection(robots[robot_id], dtau_,
//   //                                                    kkt_matrix[0], 
//   //                                                    kkt_residual[0],
//   //                                                    d[1], d[0]);
//   //     }
//   //   }
//   //   else if (i < N_) {
//   //     if (contact_sequence.existImpulseStage(i)) {
//   //       const int impulse_index = contact_sequence.impulseIndex(i);
//   //       const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
//   //       assert(dtau_impulse > 0);
//   //       assert(dtau_impulse < dtau_);
//   //       split_ocps_[i].computeKKTResidual(robots[omp_get_thread_num()], 
//   //                                         contact_sequence.contactStatus(i), 
//   //                                         t+i*dtau_, dtau_impulse, 
//   //                                         q_prev(contact_sequence, s, i),  
//   //                                         s[i], s.impulse[impulse_index],
//   //                                         kkt_matrix[i], kkt_residual[i]);
//   //     }
//   //     else if (contact_sequence.existLiftStage(i)) {
//   //       const int lift_index = contact_sequence.liftIndex(i);
//   //       const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
//   //       assert(dtau_lift > 0);
//   //       assert(dtau_lift < dtau_);
//   //       split_ocps_[i].computeKKTResidual(robots[omp_get_thread_num()], 
//   //                                         contact_sequence.contactStatus(i), 
//   //                                         t+i*dtau_, dtau_lift, 
//   //                                         q_prev(contact_sequence, s, i),  
//   //                                         s[i], s.lift[lift_index],
//   //                                         kkt_matrix[i], kkt_residual[i]);
//   //     }
//   //     else {
//   //       split_ocps_[i].computeKKTResidual(robots[omp_get_thread_num()], 
//   //                                         contact_sequence.contactStatus(i), 
//   //                                         t+i*dtau_, dtau_, 
//   //                                         q_prev(contact_sequence, s, i),  
//   //                                         s[i], s[i+1],
//   //                                         kkt_matrix[i], kkt_residual[i]);
//   //     }
//   //   }
//   //   else if (i == N_) {
//   //     terminal_ocp_.computeKKTResidual(robots[omp_get_thread_num()], t+T_, s[N_],
//   //                                      kkt_residual[N_]);
//   //   }
//   //   else if (i < N_ + 1 + N_impulse) {
//   //     const int impulse_index  = i - (N_+1);
//   //     const int time_stage_before_impulse 
//   //         = contact_sequence.timeStageBeforeImpulse(impulse_index);
//   //     split_ocps_.impulse[impulse_index].computeKKTResidual(
//   //         robots[omp_get_thread_num()], 
//   //         contact_sequence.impulseStatus(impulse_index), 
//   //         t+contact_sequence.impulseTime(impulse_index), 
//   //         s[time_stage_before_impulse].q, 
//   //         s.impulse[impulse_index], s.aux[impulse_index], 
//   //         kkt_matrix.impulse[impulse_index], kkt_residual.impulse[impulse_index]);
//   //   }
//   //   else if (i < N_ + 1 + 2*N_impulse) {
//   //     const int impulse_index  = i - (N_+1+N_impulse);
//   //     const int time_stage_before_impulse 
//   //         = contact_sequence.timeStageBeforeImpulse(impulse_index);
//   //     const double dtau_aux 
//   //         = (time_stage_before_impulse+1) * dtau_ 
//   //             - contact_sequence.impulseTime(impulse_index);
//   //     split_ocps_.aux[impulse_index].computeKKTResidual(
//   //         robots[omp_get_thread_num()], 
//   //         contact_sequence.contactStatus(time_stage_before_impulse+1), 
//   //         t+contact_sequence.impulseTime(impulse_index), dtau_aux,
//   //         s.impulse[impulse_index].q, 
//   //         s.aux[impulse_index], s[time_stage_before_impulse+1], 
//   //         kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
//   //   }
//   //   else {
//   //     const int lift_index = i - (N_+1+2*N_impulse);
//   //     const int time_stage_before_lift 
//   //         = contact_sequence.timeStageBeforeLift(lift_index);
//   //     const double dtau_aux
//   //         = (time_stage_before_lift+1) * dtau_ 
//   //             - contact_sequence.liftTime(lift_index);
//   //     split_ocps_.lift[lift_index].computeKKTResidual(
//   //         robots[omp_get_thread_num()], 
//   //         contact_sequence.contactStatus(time_stage_before_lift+1), 
//   //         t+contact_sequence.liftTime(lift_index), dtau_aux, 
//   //         s[time_stage_before_lift].q, 
//   //         s.lift[lift_index], s[time_stage_before_lift+1],
//   //         kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
//   //   }
//   // }
// }
#include "idocp/ocp/ocp_riccati_solver.hpp"


namespace idocp {

OCPRiccatiSolver::OCPRiccatiSolver(const Robot& robot, const double T, 
                                   const int N, 
                                   const int max_num_impulse, 
                                   const int num_proc) 
  : riccati_recursion_(robot, T, N, max_num_impulse, num_proc),
    riccati_factorizer_(N, max_num_impulse, robot),
    riccati_factorization_(N, max_num_impulse, robot),
    constraint_factorizer_(robot, max_num_impulse, num_proc),
    constraint_factorization_(
        max_num_impulse, StateConstraintRiccatiFactorization(robot, N, 
                                                             max_num_impulse)),
    ocp_direction_calculator_(T, N, max_num_impulse, num_proc) {
}


OCPRiccatiSolver::OCPRiccatiSolver() 
  : riccati_recursion_(),
    riccati_factorizer_(),
    riccati_factorization_(),
    constraint_factorizer_(),
    constraint_factorization_(),
    ocp_direction_calculator_() {
}


OCPRiccatiSolver::~OCPRiccatiSolver() {
}


void OCPRiccatiSolver::setConstraintDimensions(
    const ContactSequence& contact_sequence) {
  const int num_impulse_stages = contact_sequence.totalNumImpulseStages();
  for (int i=0; i<num_impulse_stages; ++i) {
    constraint_factorization_[i].setImpulseStatus(
        contact_sequence.impulseStatus(i));
  }
}


void OCPRiccatiSolver::computeDirection(HybridOCP& split_ocps, 
                                        std::vector<Robot>& robots, 
                                        const ContactSequence& contact_sequence, 
                                        const Eigen::VectorXd& q, 
                                        const Eigen::VectorXd& v, 
                                        const HybridSolution& s, 
                                        HybridDirection& d, 
                                        HybridKKTMatrix& kkt_matrix, 
                                        HybridKKTResidual& kkt_residual) {
  riccati_recursion_.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, 
                                                      riccati_factorization_);
  riccati_recursion_.backwardRiccatiRecursion(riccati_factorizer_, 
                                              contact_sequence, kkt_matrix, 
                                              kkt_residual, 
                                              riccati_factorization_);
  riccati_recursion_.forwardRiccatiRecursionParallel(riccati_factorizer_, 
                                                     contact_sequence, 
                                                     kkt_matrix, kkt_residual,
                                                     constraint_factorization_);
  riccati_recursion_.forwardRiccatiRecursionSerial(riccati_factorizer_,
                                                   contact_sequence, 
                                                   kkt_matrix, kkt_residual, 
                                                   riccati_factorization_);

  // riccati_recursion_.backwardStateConstraintFactorization(
  //     contact_sequence, kkt_matrix, constraint_factorization_);
  

  const int dimv = robots[0].dimv();

  const double T = 1.0;
  const int N = 20;
  const double dtau = T / N;
  const double dtau_impulse = contact_sequence.impulseTime(0) - 3 * dtau;
  const double dtau_aux = dtau - dtau_impulse;

  constraint_factorization_[0].T_impulse(0).topRows(dimv) = kkt_matrix.impulse[0].Pq().transpose();
  constraint_factorization_[0].T_impulse(0).bottomRows(dimv).setZero();
  constraint_factorization_[0].Eq() = kkt_matrix.impulse[0].Pq();
  constraint_factorization_[0].e() = kkt_residual.impulse[0].P();
  kkt_matrix[3].Fxx().topRightCorner(dimv, dimv) = dtau_impulse * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix[2].Fxx().topRightCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix[1].Fxx().topRightCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix[0].Fxx().topRightCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  constraint_factorization_[0].T(3) = kkt_matrix[3].Fxx().transpose() * constraint_factorization_[0].T_impulse(0);
  constraint_factorization_[0].T(2) = kkt_matrix[2].Fxx().transpose() * constraint_factorization_[0].T(3);
  constraint_factorization_[0].T(1) = kkt_matrix[1].Fxx().transpose() * constraint_factorization_[0].T(2);
  constraint_factorization_[0].T(0) = kkt_matrix[0].Fxx().transpose() * constraint_factorization_[0].T(1);

  riccati_factorization_[0].N.setZero();
  riccati_factorization_[1].N = kkt_matrix[0].Fxu() * kkt_matrix[0].Quu().inverse() * kkt_matrix[0].Fxu().transpose();
  riccati_factorization_[2].N = kkt_matrix[1].Fxu() * kkt_matrix[1].Quu().inverse() * kkt_matrix[1].Fxu().transpose() 
                                + kkt_matrix[1].Fxx() * riccati_factorization_[1].N * kkt_matrix[1].Fxx().transpose();
  riccati_factorization_[3].N = kkt_matrix[2].Fxu() * kkt_matrix[2].Quu().inverse() * kkt_matrix[2].Fxu().transpose() 
                                + kkt_matrix[2].Fxx() * riccati_factorization_[2].N * kkt_matrix[2].Fxx().transpose();
  riccati_factorization_.impulse[0].N = kkt_matrix[3].Fxu() * kkt_matrix[3].Quu().inverse() * kkt_matrix[3].Fxu().transpose() 
                                        + kkt_matrix[3].Fxx() * riccati_factorization_[3].N * kkt_matrix[3].Fxx().transpose();

  assert(constraint_factorization_[0].T_impulse(0).topRows(dimv).isApprox(kkt_matrix.impulse[0].Pq().transpose()));
  assert(constraint_factorization_[0].T_impulse(0).bottomRows(dimv).isZero());
  assert(constraint_factorization_[0].Eq().isApprox(kkt_matrix.impulse[0].Pq()));
  assert(constraint_factorization_[0].e().isApprox(kkt_residual.impulse[0].P()));

  ocp_direction_calculator_.computeInitialStateDirection(robots, q, v, s, d);
  if (contact_sequence.existImpulseStage()) {
    // constraint_factorizer_.computeLagrangeMultiplierDirection(
    //     contact_sequence, riccati_factorization_.impulse, 
    //     constraint_factorization_, d[0].dx(), d.impulse);

    const Eigen::MatrixXd E = constraint_factorization_[0].T_impulse(0).transpose();
    assert(E.leftCols(dimv).isApprox(kkt_matrix.impulse[0].Pq()));
    assert(E.rightCols(dimv).isZero());
    const Eigen::VectorXd vec_term = E * riccati_factorization_.impulse[0].Pi * d[0].dx() 
                                      + E * riccati_factorization_.impulse[0].pi
                                      + kkt_residual.impulse[0].P();
    const Eigen::MatrixXd ENEt = E * riccati_factorization_.impulse[0].N * E.transpose();
    const Eigen::MatrixXd ENEtinv = ENEt.inverse();
    d.impulse[0].dxi() = ENEtinv * vec_term;

    for (int i=0; i<=3; ++i) {
      riccati_factorization_[i].n.noalias() = constraint_factorization_[0].T(i) * d.impulse[0].dxi();
    }
    riccati_factorization_.impulse[0].n.noalias() 
        = constraint_factorization_[0].T_impulse(0) * d.impulse[0].dxi();
    riccati_factorization_.aux[0].n.setZero();
    for (int i=4; i<=20; ++i) {
      riccati_factorization_[i].n.setZero();
    }

    // ocp_direction_calculator_.aggregateLagrangeMultiplierDirection(
    //     contact_sequence, riccati_factorization_, constraint_factorization_, d);
  }
  ocp_direction_calculator_.computeDirection(split_ocps, robots, 
                                             contact_sequence, 
                                             riccati_factorizer_, 
                                             riccati_factorization_, 
                                             constraint_factorization_, s, d);

  for (int i=0; i<=2; ++i) {
    d[i].du() -= kkt_matrix[i].Quu().inverse() * kkt_matrix[i].Fxu().transpose() * riccati_factorization_[i+1].n;
  }
  d[3].du() -= kkt_matrix[3].Quu().inverse() * kkt_matrix[3].Fxu().transpose() * riccati_factorization_.impulse[0].n;
  for (int i=1; i<=3; ++i) {
    d[i].dx() -= riccati_factorization_[i].N * riccati_factorization_[i].n;
  }
  for (int i=0; i<=3; ++i) {
    d[i].dlmdgmm() += riccati_factorization_[i].n;
  }
  d.impulse[0].dx() -= riccati_factorization_.impulse[0].N * riccati_factorization_.impulse[0].n;
  d.impulse[0].dlmdgmm() += riccati_factorization_.impulse[0].n;

}


double OCPRiccatiSolver::maxPrimalStepSize(
    const ContactSequence& contact_sequence) const {
  return ocp_direction_calculator_.maxPrimalStepSize(contact_sequence);
}


double OCPRiccatiSolver::maxDualStepSize(
    const ContactSequence& contact_sequence) const {
  return ocp_direction_calculator_.maxDualStepSize(contact_sequence);
}

} // namespace idocp
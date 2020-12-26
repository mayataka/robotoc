#include "derivative_checker.hpp"

#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"

#include <cmath>


namespace idocp {

DerivativeChecker::DerivativeChecker(const Robot& robot, 
                                     const double finite_diff, 
                                     const double test_tol)
  : robot_(robot), 
    finite_diff_(finite_diff), 
    test_tol_(test_tol) {
}


DerivativeChecker::~DerivativeChecker() {
}


void DerivativeChecker::setFiniteDifference(const double finite_diff) {
  finite_diff_ = finite_diff;
}


void DerivativeChecker::setTestTolerance(const double test_tol) {
  test_tol_ = test_tol;
}


bool DerivativeChecker::checkFirstOrderStageCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost, 
    const ContactStatus& contact_status) {
  const SplitSolution s = SplitSolution::Random(robot_, contact_status);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  const int dimv = robot_.dimv();
  const int dimu = robot_.dimu();
  const int dimf = contact_status.dimf();
  SplitKKTResidual kkt_residual(robot_);
  kkt_residual.setContactStatus(contact_status);
  CostFunctionData data(robot_);
  cost->computeStageCostDerivatives(robot_, data, t, dtau, s, kkt_residual);
  double cost0 = cost->computeStageCost(robot_, data, t, dtau, s);
  SplitSolution s1 = s;
  Eigen::VectorXd lq_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    lq_ref(i) = (cost->computeStageCost(robot_, data, t, dtau, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lq().isApprox(lq_ref, test_tol_)) {
    std::cout << "lq is not correct! lq - lq_ref = " 
              << (kkt_residual.lq() - lq_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd lv_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.v(i) += finite_diff_;
    lv_ref(i) = (cost->computeStageCost(robot_, data, t, dtau, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lv().isApprox(lv_ref, test_tol_)) {
    std::cout << "lv is not correct! lv - lv_ref = " 
              << (kkt_residual.lv() - lv_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd la_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.a(i) += finite_diff_;
    la_ref(i) = (cost->computeStageCost(robot_, data, t, dtau, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.la.isApprox(la_ref, test_tol_)) {
    std::cout << "la is not correct! la - la_ref = " 
              << (kkt_residual.la - la_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd lu_ref(dimu);
  for (int i=0; i<dimu; ++i) {
    s1 = s;
    s1.u(i) += finite_diff_;
    lu_ref(i) = (cost->computeStageCost(robot_, data, t, dtau, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lu().isApprox(lu_ref, test_tol_)) {
    std::cout << "lu is not correct! lu - lu_ref = " 
              << (kkt_residual.lu() - lu_ref).transpose() << std::endl;
    return false;
  }
  if (dimf > 0) {
    Eigen::VectorXd lf_ref(dimf);
    for (int i=0; i<dimf; ++i) {
      s1 = s;
      s1.f_stack().coeffRef(i) += finite_diff_;
      lf_ref(i) = (cost->computeStageCost(robot_, data, t, dtau, s1) - cost0) / finite_diff_;
    }
    if (!kkt_residual.lf().isApprox(lf_ref, test_tol_)) {
      std::cout << "lf is not correct! lf - lf_ref = " 
                << (kkt_residual.lf() - lf_ref).transpose() << std::endl;
      return false;
    }
  }
  return true;
}


bool DerivativeChecker::checkSecondOrderStageCostDerivatives(
  const std::shared_ptr<CostFunctionComponentBase>& cost,
  const ContactStatus& conatct_status) {
}


bool DerivativeChecker::checkFirstOrderTerminalCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost) {
  const SplitSolution s = SplitSolution::Random(robot_);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const int dimv = robot_.dimv();
  SplitKKTResidual kkt_residual(robot_);
  CostFunctionData data(robot_);
  cost->computeTerminalCostDerivatives(robot_, data, t, s, kkt_residual);
  double cost0 = cost->computeTerminalCost(robot_, data, t, s);
  SplitSolution s1 = s;
  Eigen::VectorXd lq_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    lq_ref(i) = (cost->computeTerminalCost(robot_, data, t, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lq().isApprox(lq_ref, test_tol_)) {
    std::cout << "lq is not correct! lq - lq_ref = " 
              << (kkt_residual.lq() - lq_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd lv_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.v(i) += finite_diff_;
    lv_ref(i) = (cost->computeTerminalCost(robot_, data, t, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lv().isApprox(lv_ref, test_tol_)) {
    std::cout << "lv is not correct! lv - lv_ref = " 
              << (kkt_residual.lv() - lv_ref).transpose() << std::endl;
    return false;
  }
  return true;
}


bool DerivativeChecker::checkSecondOrderTerminalCostDerivatives(
  const std::shared_ptr<CostFunctionComponentBase>& cost) {
}


bool DerivativeChecker::checkFirstOrderImpulseCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost, 
    const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot_, impulse_status);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const int dimv = robot_.dimv();
  const int dimf = impulse_status.dimf();
  ImpulseSplitKKTResidual kkt_residual(robot_);
  kkt_residual.setImpulseStatus(impulse_status);
  CostFunctionData data(robot_);
  cost->computeImpulseCostDerivatives(robot_, data, t, s, kkt_residual);
  double cost0 = cost->computeImpulseCost(robot_, data, t, s);
  ImpulseSplitSolution s1 = s;
  Eigen::VectorXd lq_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    lq_ref(i) = (cost->computeImpulseCost(robot_, data, t, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lq().isApprox(lq_ref, test_tol_)) {
    std::cout << "lq is not correct! lq - lq_ref = " 
              << (kkt_residual.lq() - lq_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd lv_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.v(i) += finite_diff_;
    lv_ref(i) = (cost->computeImpulseCost(robot_, data, t, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lv().isApprox(lv_ref, test_tol_)) {
    std::cout << "lv is not correct! lv - lv_ref = " 
              << (kkt_residual.lv() - lv_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd ldv_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.dv(i) += finite_diff_;
    ldv_ref(i) = (cost->computeImpulseCost(robot_, data, t, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.ldv.isApprox(ldv_ref, test_tol_)) {
    std::cout << "la is not correct! ldv - ldv_ref = " 
              << (kkt_residual.ldv - ldv_ref).transpose() << std::endl;
    return false;
  }
  if (dimf > 0) {
    Eigen::VectorXd lf_ref(dimf);
    for (int i=0; i<dimf; ++i) {
      s1 = s;
      s1.f_stack().coeffRef(i) += finite_diff_;
      lf_ref(i) = (cost->computeImpulseCost(robot_, data, t, s1) - cost0) / finite_diff_;
    }
    if (!kkt_residual.lf().isApprox(lf_ref, test_tol_)) {
      std::cout << "lf is not correct! lf - lf_ref = " 
                << (kkt_residual.lf() - lf_ref).transpose() << std::endl;
      return false;
    }
  }
  return true;
}


bool DerivativeChecker::checkSecondOrderImpulseCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost,
    const ImpulseStatus& impulse_status) {

}

} // namespace idocp 
#include "robotoc/utils/derivative_checker.hpp"

#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/hybrid/grid_info.hpp"

#include <cmath>


namespace robotoc {

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
    const std::shared_ptr<CostFunctionComponentBase>& cost) {
  return checkFirstOrderStageCostDerivatives(cost, robot_.createContactStatus());
}


bool DerivativeChecker::checkSecondOrderStageCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost) {
  return checkSecondOrderStageCostDerivatives(cost, robot_.createContactStatus());
}


bool DerivativeChecker::checkFirstOrderStageCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost, 
    const ContactStatus& contact_status) {
  const auto s = SplitSolution::Random(robot_, contact_status);
  const auto grid_info = GridInfo::Random();
  const int dimv = robot_.dimv();
  const int dimu = robot_.dimu();
  const int dimf = contact_status.dimf();
  SplitKKTResidual kkt_residual(robot_);
  kkt_residual.setContactStatus(contact_status);
  CostFunctionData data(robot_);
  robot_.updateKinematics(s.q, s.v, s.a);
  double cost0 = cost->evalStageCost(robot_, contact_status, data, grid_info, s);
  cost->evalStageCostDerivatives(robot_, contact_status, data, grid_info, s, kkt_residual);
  auto s1 = s;
  Eigen::VectorXd lq_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    robot_.updateKinematics(s1.q, s1.v, s1.a);
    lq_ref(i) = (cost->evalStageCost(robot_, contact_status, data, grid_info, s1) - cost0) / finite_diff_;
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
    robot_.updateKinematics(s1.q, s1.v, s1.a);
    lv_ref(i) = (cost->evalStageCost(robot_, contact_status, data, grid_info, s1) - cost0) / finite_diff_;
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
    robot_.updateKinematics(s1.q, s1.v, s1.a);
    la_ref(i) = (cost->evalStageCost(robot_, contact_status, data, grid_info, s1) - cost0) / finite_diff_;
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
    lu_ref(i) = (cost->evalStageCost(robot_, contact_status, data, grid_info, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.lu.isApprox(lu_ref, test_tol_)) {
    std::cout << "lu is not correct! lu - lu_ref = " 
              << (kkt_residual.lu - lu_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd lf_ref(dimf);
  if (dimf > 0) {
    for (int i=0; i<dimf; ++i) {
      s1 = s;
      s1.f_stack().coeffRef(i) += finite_diff_;
      s1.set_f_vector();
      lf_ref(i) = (cost->evalStageCost(robot_, contact_status, data, grid_info, s1) - cost0) / finite_diff_;
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
  const ContactStatus& contact_status) {
  const auto s = SplitSolution::Random(robot_, contact_status);
  const auto grid_info = GridInfo::Random();
  const int dimv = robot_.dimv();
  const int dimu = robot_.dimu();
  const int dimf = contact_status.dimf();
  SplitKKTResidual kkt_residual0(robot_);
  kkt_residual0.setContactStatus(contact_status);
  SplitKKTMatrix kkt_matrix(robot_);
  kkt_matrix.setContactStatus(contact_status);
  CostFunctionData data(robot_);
  robot_.updateKinematics(s.q, s.v, s.a);
  cost->evalStageCost(robot_, contact_status, data, grid_info, s);
  cost->evalStageCostDerivatives(robot_, contact_status, data, grid_info, s, kkt_residual0);
  cost->evalStageCostHessian(robot_, contact_status, data, grid_info, s, kkt_matrix);
  SplitKKTResidual kkt_residual(robot_);
  kkt_residual.setContactStatus(contact_status);
  auto s1 = s;
  Eigen::MatrixXd Qqq_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    robot_.updateKinematics(s1.q, s1.v, s1.a);
    kkt_residual.lq().setZero();
    cost->evalStageCost(robot_, contact_status, data, grid_info, s1);
    cost->evalStageCostDerivatives(robot_, contact_status, data, grid_info, s1, kkt_residual);
    Qqq_ref.col(i) = (kkt_residual.lq() - kkt_residual0.lq()) / finite_diff_;
  }
  if (!kkt_matrix.Qqq().isApprox(Qqq_ref, test_tol_)) {
    std::cout << "Qqq is not correct! Qqq - Qqq_ref = " 
              << (kkt_matrix.Qqq() - Qqq_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Qvv_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.v(i) += finite_diff_;
    robot_.updateKinematics(s1.q, s1.v, s1.a);
    kkt_residual.lv().setZero();
    cost->evalStageCost(robot_, contact_status, data, grid_info, s1);
    cost->evalStageCostDerivatives(robot_, contact_status, data, grid_info, s1, kkt_residual);
    Qvv_ref.col(i) = (kkt_residual.lv() - kkt_residual0.lv()) / finite_diff_;
  }
  if (!kkt_matrix.Qvv().isApprox(Qvv_ref, test_tol_)) {
    std::cout << "Qvv is not correct! Qvv - Qvv_ref = " 
              << (kkt_matrix.Qvv() - Qvv_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Qaa_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.a(i) += finite_diff_;
    robot_.updateKinematics(s1.q, s1.v, s1.a);
    kkt_residual.la.setZero();
    cost->evalStageCost(robot_, contact_status, data, grid_info, s1);
    cost->evalStageCostDerivatives(robot_, contact_status, data, grid_info, s1, kkt_residual);
    Qaa_ref.col(i) = (kkt_residual.la - kkt_residual0.la) / finite_diff_;
  }
  if (!kkt_matrix.Qaa.isApprox(Qaa_ref, test_tol_)) {
    std::cout << "Qaa is not correct! Qaa - Qaa_ref = " 
              << (kkt_matrix.Qaa - Qaa_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Quu_ref(dimu, dimu);
  for (int i=0; i<dimu; ++i) {
    s1 = s;
    s1.u(i) += finite_diff_;
    kkt_residual.lu.setZero();
    cost->evalStageCost(robot_, contact_status, data, grid_info, s1);
    cost->evalStageCostDerivatives(robot_, contact_status, data, grid_info, s1, kkt_residual);
    Quu_ref.col(i) = (kkt_residual.lu - kkt_residual0.lu) / finite_diff_;
  }
  if (!kkt_matrix.Quu.isApprox(Quu_ref, test_tol_)) {
    std::cout << "Quu is not correct! Quu - Quu_ref = " 
              << (kkt_matrix.Quu - Quu_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Qff_ref(dimf, dimf);
  if (dimf > 0) {
    for (int i=0; i<dimf; ++i) {
      s1 = s;
      s1.f_stack().coeffRef(i) += finite_diff_;
      s1.set_f_vector();
      kkt_residual.lf().setZero();
      cost->evalStageCost(robot_, contact_status, data, grid_info, s1);
      cost->evalStageCostDerivatives(robot_, contact_status, data, grid_info, s1, kkt_residual);
      Qff_ref.col(i) = (kkt_residual.lf() - kkt_residual0.lf()) / finite_diff_;
    }
    if (!kkt_matrix.Qff().isApprox(Qff_ref, test_tol_)) {
      std::cout << "Qff is not correct! Qff - Qff_ref = " 
                << (kkt_matrix.Qff() - Qff_ref).transpose() << std::endl;
      return false;
    }
  }
  return true;
}


bool DerivativeChecker::checkFirstOrderTerminalCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost) {
  const auto s = SplitSolution::Random(robot_);
  const auto grid_info = GridInfo::Random();
  const int dimv = robot_.dimv();
  SplitKKTResidual kkt_residual(robot_);
  CostFunctionData data(robot_);
  robot_.updateKinematics(s.q, s.v);
  double cost0 = cost->evalTerminalCost(robot_, data, grid_info, s);
  cost->evalTerminalCostDerivatives(robot_, data, grid_info, s, kkt_residual);
  auto s1 = s;
  Eigen::VectorXd lq_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    robot_.updateKinematics(s1.q, s1.v);
    lq_ref(i) = (cost->evalTerminalCost(robot_, data, grid_info, s1) - cost0) / finite_diff_;
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
    robot_.updateKinematics(s1.q, s1.v);
    lv_ref(i) = (cost->evalTerminalCost(robot_, data, grid_info, s1) - cost0) / finite_diff_;
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
  const auto s = SplitSolution::Random(robot_);
  const auto grid_info = GridInfo::Random();
  const int dimv = robot_.dimv();
  SplitKKTMatrix kkt_matrix(robot_);
  SplitKKTResidual kkt_residual0(robot_);
  CostFunctionData data(robot_);
  robot_.updateKinematics(s.q, s.v);
  cost->evalTerminalCost(robot_, data, grid_info, s);
  cost->evalTerminalCostDerivatives(robot_, data, grid_info, s, kkt_residual0);
  cost->evalTerminalCostHessian(robot_, data, grid_info, s, kkt_matrix);
  SplitKKTResidual kkt_residual(robot_);
  auto s1 = s;
  Eigen::MatrixXd Qqq_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    robot_.updateKinematics(s1.q, s1.v);
    kkt_residual.lq().setZero();
    cost->evalTerminalCost(robot_, data, grid_info, s1);
    cost->evalTerminalCostDerivatives(robot_, data, grid_info, s1, kkt_residual);
    Qqq_ref.col(i) = (kkt_residual.lq() - kkt_residual0.lq()) / finite_diff_;
  }
  if (!kkt_matrix.Qqq().isApprox(Qqq_ref, test_tol_)) {
    std::cout << "Qqq is not correct! Qqq - Qqq_ref = " 
              << (kkt_matrix.Qqq() - Qqq_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Qvv_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.v(i) += finite_diff_;
    robot_.updateKinematics(s1.q, s1.v);
    kkt_residual.lv().setZero();
    cost->evalTerminalCost(robot_, data, grid_info, s1);
    cost->evalTerminalCostDerivatives(robot_, data, grid_info, s1, kkt_residual);
    Qvv_ref.col(i) = (kkt_residual.lv() - kkt_residual0.lv()) / finite_diff_;
  }
  if (!kkt_matrix.Qvv().isApprox(Qvv_ref, test_tol_)) {
    std::cout << "Qvv is not correct! Qvv - Qvv_ref = " 
              << (kkt_matrix.Qvv() - Qvv_ref).transpose() << std::endl;
    return false;
  }
  return true;
}


bool DerivativeChecker::checkFirstOrderImpulseCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost) {
  return checkFirstOrderImpulseCostDerivatives(cost, robot_.createImpulseStatus());
}


bool DerivativeChecker::checkSecondOrderImpulseCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost) {
  return checkSecondOrderImpulseCostDerivatives(cost, robot_.createImpulseStatus());
}


bool DerivativeChecker::checkFirstOrderImpulseCostDerivatives(
    const std::shared_ptr<CostFunctionComponentBase>& cost, 
    const ImpulseStatus& impulse_status) {
  const auto s = ImpulseSplitSolution::Random(robot_, impulse_status);
  const auto grid_info = GridInfo::Random();
  const int dimv = robot_.dimv();
  const int dimf = impulse_status.dimi();
  ImpulseSplitKKTResidual kkt_residual(robot_);
  kkt_residual.setImpulseStatus(impulse_status);
  CostFunctionData data(robot_);
  robot_.updateKinematics(s.q, s.v);
  double cost0 = cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s);
  cost->evalImpulseCostDerivatives(robot_, impulse_status, data, grid_info, s, kkt_residual);
  auto s1 = s;
  Eigen::VectorXd lq_ref(dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    robot_.updateKinematics(s1.q, s1.v);
    lq_ref(i) = (cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1) - cost0) / finite_diff_;
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
    robot_.updateKinematics(s1.q, s1.v);
    lv_ref(i) = (cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1) - cost0) / finite_diff_;
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
    ldv_ref(i) = (cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1) - cost0) / finite_diff_;
  }
  if (!kkt_residual.ldv.isApprox(ldv_ref, test_tol_)) {
    std::cout << "ldv is not correct! ldv - ldv_ref = " 
              << (kkt_residual.ldv - ldv_ref).transpose() << std::endl;
    return false;
  }
  Eigen::VectorXd lf_ref(dimf);
  if (dimf > 0) {
    for (int i=0; i<dimf; ++i) {
      s1 = s;
      s1.f_stack().coeffRef(i) += finite_diff_;
      s1.set_f_vector();
      lf_ref(i) = (cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1) - cost0) / finite_diff_;
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
  const auto s = ImpulseSplitSolution::Random(robot_, impulse_status);
  const auto grid_info = GridInfo::Random();
  const int dimv = robot_.dimv();
  const int dimf = impulse_status.dimi();
  ImpulseSplitKKTMatrix kkt_matrix(robot_);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual0(robot_);
  kkt_residual0.setImpulseStatus(impulse_status);
  CostFunctionData data(robot_);
  robot_.updateKinematics(s.q, s.v);
  cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s);
  cost->evalImpulseCostDerivatives(robot_, impulse_status, data, grid_info, s, kkt_residual0);
  cost->evalImpulseCostHessian(robot_, impulse_status, data, grid_info, s, kkt_matrix);
  ImpulseSplitKKTResidual kkt_residual(robot_);
  kkt_residual.setImpulseStatus(impulse_status);
  auto s1 = s;
  Eigen::MatrixXd Qqq_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv);
    dq(i) = 1;
    robot_.integrateConfiguration(s.q, dq, finite_diff_, s1.q);
    kkt_residual.lq().setZero();
    robot_.updateKinematics(s1.q, s1.v);
    cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1);
    cost->evalImpulseCostDerivatives(robot_, impulse_status, data, grid_info, s1, kkt_residual);
    Qqq_ref.col(i) = (kkt_residual.lq() - kkt_residual0.lq()) / finite_diff_;
  }
  if (!kkt_matrix.Qqq().isApprox(Qqq_ref, test_tol_)) {
    std::cout << "Qqq is not correct! Qqq - Qqq_ref = " 
              << (kkt_matrix.Qqq() - Qqq_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Qvv_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.v(i) += finite_diff_;
    kkt_residual.lv().setZero();
    robot_.updateKinematics(s1.q, s1.v);
    cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1);
    cost->evalImpulseCostDerivatives(robot_, impulse_status, data, grid_info, s1, kkt_residual);
    Qvv_ref.col(i) = (kkt_residual.lv() - kkt_residual0.lv()) / finite_diff_;
  }
  if (!kkt_matrix.Qvv().isApprox(Qvv_ref, test_tol_)) {
    std::cout << "Qvv is not correct! Qvv - Qvv_ref = " 
              << (kkt_matrix.Qvv() - Qvv_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Qdvdv_ref(dimv, dimv);
  for (int i=0; i<dimv; ++i) {
    s1 = s;
    s1.dv(i) += finite_diff_;
    kkt_residual.ldv.setZero();
    cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1);
    cost->evalImpulseCostDerivatives(robot_, impulse_status, data, grid_info, s1, kkt_residual);
    Qdvdv_ref.col(i) = (kkt_residual.ldv - kkt_residual0.ldv) / finite_diff_;
  }
  if (!kkt_matrix.Qdvdv.isApprox(Qdvdv_ref, test_tol_)) {
    std::cout << "Qdvdv is not correct! Qdvdv - Qdvdv_ref = " 
              << (kkt_matrix.Qdvdv - Qdvdv_ref).transpose() << std::endl;
    return false;
  }
  Eigen::MatrixXd Qff_ref(dimf, dimf);
  if (dimf > 0) {
    for (int i=0; i<dimf; ++i) {
      s1 = s;
      s1.f_stack().coeffRef(i) += finite_diff_;
      s1.set_f_vector();
      kkt_residual.lf().setZero();
      cost->evalImpulseCost(robot_, impulse_status, data, grid_info, s1);
      cost->evalImpulseCostDerivatives(robot_, impulse_status, data, grid_info, s1, kkt_residual);
      Qff_ref.col(i) = (kkt_residual.lf() - kkt_residual0.lf()) / finite_diff_;
    }
    if (!kkt_matrix.Qff().isApprox(Qff_ref, test_tol_)) {
      std::cout << "Qff is not correct! Qff - Qff_ref = " 
                << (kkt_matrix.Qff() - Qff_ref).transpose() << std::endl;
      return false;
    }
  }
  return true;
}

} // namespace robotoc 
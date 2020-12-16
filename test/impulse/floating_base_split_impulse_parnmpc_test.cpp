#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/split_impulse_parnmpc.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/impulse_constraints.hpp"


namespace idocp {

class FloatingBaseSplitParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/anymal/anymal.urdf";
    std::vector<int> contact_frames = {14, 24, 34, 44};
    robot = Robot(urdf, contact_frames);
    std::random_device rnd;
    std::vector<bool> is_contact_active;
    for (const auto frame : contact_frames) {
      is_contact_active.push_back(rnd()%2==0);
    }
    contact_status = ContactStatus(robot.maxPointContacts());
    contact_status.setContactStatus(is_contact_active);
    if (!contact_status.hasActiveContacts()) {
      contact_status.activateContact(0);
    }
    s = ImpulseSplitSolution::Random(robot, contact_status);
    s_prev = SplitSolution::Random(robot, contact_status);
    s_prev_new = SplitSolution::Random(robot, contact_status);
    s_next = SplitSolution::Random(robot, contact_status);
    s_next_new = SplitSolution::Random(robot, contact_status);
    s_tmp = ImpulseSplitSolution::Random(robot, contact_status);
    s_new = ImpulseSplitSolution::Random(robot, contact_status);
    d = ImpulseSplitDirection::Random(robot, contact_status);
    d_ref = d;
    d_prev = SplitDirection::Random(robot, contact_status);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    q_prev = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_prev);
    v_prev = Eigen::VectorXd::Random(robot.dimv());
    dq_prev = Eigen::VectorXd::Random(robot.dimv());
    dv_prev = Eigen::VectorXd::Random(robot.dimv());
    robot.normalizeConfiguration(s_next.q);
    auto joint_cost = std::make_shared<JointSpaceImpulseCost>(robot);
    auto impulse_cost = std::make_shared<ImpulseForceCost>(robot);
    const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_ref);
    const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd dv_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd dv_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    std::vector<Eigen::Vector3d> f_weight, f_ref;
    for (int i=0; i<robot.maxPointContacts(); ++i) {
      f_weight.push_back(Eigen::Vector3d::Random());
      f_ref.push_back(Eigen::Vector3d::Random());
    }
    joint_cost->set_q_weight(q_weight);
    joint_cost->set_q_ref(q_ref);
    joint_cost->set_v_weight(v_weight);
    joint_cost->set_v_ref(v_ref);
    joint_cost->set_dv_weight(dv_weight);
    joint_cost->set_dv_ref(dv_ref);
    impulse_cost->set_f_weight(f_weight);
    impulse_cost->set_f_ref(f_ref);
    cost = std::make_shared<ImpulseCostFunction>();
    cost->push_back(joint_cost);
    cost->push_back(impulse_cost);
    cost_data = CostFunctionData(robot);
    constraints = std::make_shared<ImpulseConstraints>();
    kkt_matrix = ImpulseSplitKKTMatrix(robot);
    kkt_residual = ImpulseSplitKKTResidual(robot);
    kkt_matrix.setContactStatus(contact_status);
    kkt_residual.setContactStatus(contact_status);
    impulse_dynamics = ImpulseDynamicsBackwardEuler(robot);
  }

  virtual void TearDown() {
  }

  double t;
  std::string urdf;
  Robot robot;
  ContactStatus contact_status;
  std::shared_ptr<ImpulseCostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<ImpulseConstraints> constraints;
  ConstraintsData constraints_data;
  ImpulseSplitSolution s, s_tmp, s_new;
  SplitSolution s_prev, s_prev_new, s_next, s_next_new;
  ImpulseSplitDirection d, d_ref;
  SplitDirection d_prev;
  ImpulseSplitKKTMatrix kkt_matrix;
  ImpulseSplitKKTResidual kkt_residual;
  ImpulseDynamicsBackwardEuler impulse_dynamics;
  Eigen::VectorXd q_prev, v_prev, dq_prev, dv_prev;
};


TEST_F(FloatingBaseSplitParNMPCTest, isFeasible) {
  SplitImpulseParNMPC parnmpc(robot, cost, constraints);
  EXPECT_EQ(parnmpc.isFeasible(robot, s), 
            constraints->isFeasible(robot, constraints_data, s));
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormStateEquation) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  SplitImpulseParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  s.beta.setZero();
  s.mu_stack().setZero();
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const double kkt_error = parnmpc.squaredNormKKTResidual();
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  robot.setContactForces(contact_status, s.f);
  robot.updateKinematics(s.q, s.v);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.ldv.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.dv_res.squaredNorm()
                         + kkt_residual.C().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  const auto pair = parnmpc.costAndConstraintViolation(robot, t, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormStateEquationAndRobotDynamics) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  SplitImpulseParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const double kkt_error = parnmpc.squaredNormKKTResidual();
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  Eigen::MatrixXd ddv_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd ddv_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd ddv_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.RNEAImpulseDerivatives(s.q, s.dv, ddv_dq, ddv_ddv);
  robot.updateKinematics(s.q, s.v);
  robot.dRNEAPartialdFext(contact_status, ddv_df);
  kkt_residual.lq() += ddv_dq.transpose() * s.beta;
  kkt_residual.ldv += ddv_ddv.transpose() * s.beta;
  kkt_residual.lf() += ddv_df.transpose() * s.beta;
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix.Cq(), kkt_matrix.Cv());
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_stack();
  kkt_residual.ldv += kkt_matrix.Cv().transpose() * s.mu_stack();
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.ldv.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.dv_res.squaredNorm()
                         + kkt_residual.C().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  const auto pair = parnmpc.costAndConstraintViolation(robot, t, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + kkt_residual.dv_res.lpNorm<1>()
                                + kkt_residual.C().lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  SplitImpulseParNMPC parnmpc(robot, empty_cost, constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const double kkt_error = parnmpc.squaredNormKKTResidual();
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics.linearizeImpulseDynamics(robot, contact_status, s, 
                                            kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.ldv.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += impulse_dynamics.squaredNormImpulseDynamicsResidual(kkt_residual);
  kkt_error_ref += constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  SplitImpulseParNMPC parnmpc(robot, cost, empty_constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const double kkt_error = parnmpc.squaredNormKKTResidual();
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics.linearizeImpulseDynamics(robot, contact_status, s, 
                                            kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.ldv.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += impulse_dynamics.squaredNormImpulseDynamicsResidual(kkt_residual);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNorm) {
  SplitImpulseParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const double kkt_error = parnmpc.squaredNormKKTResidual();
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics.linearizeImpulseDynamics(robot, contact_status, s, 
                                            kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.ldv.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += impulse_dynamics.squaredNormImpulseDynamicsResidual(kkt_residual);
  kkt_error_ref += constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, costAndConstraintViolation) {
  SplitImpulseParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const auto pair = parnmpc.costAndConstraintViolation(robot, t, s); 
  const double cost_ref 
      = cost->l(robot, cost_data, t, s) + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fv() = v_prev - s.v + s.dv;
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.updateKinematics(s.q, s.v);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  const double violation_ref 
      = kkt_residual.Fx().lpNorm<1>() 
          + kkt_residual.dv_res.lpNorm<1>() 
          + kkt_residual.C().lpNorm<1>()
          + constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, costAndConstraintViolationWithStepSizeInitial) {
  SplitImpulseParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const double step_size = 0.3;
  const auto pair = parnmpc.costAndConstraintViolation(robot, contact_status, 
                                                       step_size, t, s_prev.q, 
                                                       s_prev.v, s, d, s_tmp); 
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.dv = s.dv + step_size * d.ddv;
  s_tmp.f_stack() = s.f_stack() + step_size * d.df();
  s_tmp.set_f();
  constraints->computePrimalAndDualResidual(robot, constraints_data, s_tmp);
  const double cost_ref 
      = cost->l(robot, cost_data, t, s_tmp) 
          + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s_prev.q, s_tmp.q, kkt_residual.Fq());
  kkt_residual.Fv() = s_prev.v - s_tmp.v + s_tmp.dv;
  robot.setContactForces(contact_status, s_tmp.f);
  robot.RNEAImpulse(s_tmp.q, s_tmp.dv, kkt_residual.dv_res);
  robot.updateKinematics(s_tmp.q, s_tmp.v);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fx().lpNorm<1>()  
          + kkt_residual.dv_res.lpNorm<1>() 
          + kkt_residual.C().lpNorm<1>()
          + constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, costAndConstraintViolationWithStepSize) {
  SplitImpulseParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  parnmpc.computeKKTResidual(robot, contact_status, t, q_prev, v_prev, s, s_next);
  const double step_size = 0.3;
  const auto pair = parnmpc.costAndConstraintViolation(robot, contact_status, 
                                                       step_size, t, s_prev, 
                                                       d_prev, s, d, s_tmp); 
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.dv = s.dv + step_size * d.ddv;
  s_tmp.f_stack() = s.f_stack() + step_size * d.df();
  s_tmp.set_f();
  constraints->computePrimalAndDualResidual(robot, constraints_data, s_tmp);
  const double cost_ref 
      = cost->l(robot, cost_data, t, s_tmp) 
          + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s_prev.q, s_tmp.q, kkt_residual.Fq());
  kkt_residual.Fq() += step_size * d_prev.dq();
  kkt_residual.Fv() = s_prev.v + step_size * d_prev.dv() - s_tmp.v + s_tmp.dv;
  robot.setContactForces(contact_status, s_tmp.f);
  robot.RNEAImpulse(s_tmp.q, s_tmp.dv, kkt_residual.dv_res);
  robot.updateKinematics(s_tmp.q, s_tmp.v);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fx().lpNorm<1>()  
          + kkt_residual.dv_res.lpNorm<1>() 
          + kkt_residual.C().lpNorm<1>()
          + constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, coarseUpdate) {
  SplitImpulseParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const int dimx = 2 * robot.dimv();
  Eigen::MatrixXd aux_mat_seed = Eigen::MatrixXd::Random(dimx, dimx);
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed.transpose() * aux_mat_seed;
  ImpulseSplitSolution s_new_coarse(robot);
  s_new_coarse.setContactStatus(contact_status);
  // coarse update 
  parnmpc.coarseUpdate(robot, contact_status, t, q_prev, v_prev, s, s_next, 
                       aux_mat_next, d, s_new_coarse);
  Eigen::MatrixXd aux_mat = Eigen::MatrixXd::Zero(dimx, dimx);
  parnmpc.getAuxiliaryMatrix(aux_mat);
  // coarse update ref
  cost->computeStageCostHessian(robot, cost_data, t, s, kkt_matrix);
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics.condenseImpulseDynamics(robot, contact_status, s, 
                                           kkt_matrix, kkt_residual);
  kkt_matrix.Qxx().noalias() += aux_mat_next;
  kkt_matrix.symmetrize(); 
  const int dimKKT = kkt_residual.dimKKT();
  Eigen::MatrixXd kkt_matrix_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  kkt_matrix.invert(kkt_matrix_inv);
  d_ref.split_direction() = kkt_matrix_inv * kkt_residual.KKT_residual();
  ImpulseSplitSolution s_new_coarse_ref(robot);
  s_new_coarse_ref.setContactStatus(contact_status);
  s_new_coarse_ref.lmd = s.lmd - d_ref.dlmd();
  s_new_coarse_ref.gmm = s.gmm - d_ref.dgmm();
  s_new_coarse_ref.mu_stack() = s.mu_stack() - d_ref.dmu();
  s_new_coarse_ref.f_stack() = s.f_stack() - d_ref.df();
  robot.integrateConfiguration(s.q, d_ref.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v = s.v - d_ref.dv();
  Eigen::MatrixXd aux_mat_ref = - kkt_matrix_inv.topLeftCorner(dimx, dimx);
  EXPECT_TRUE(aux_mat.isApprox(aux_mat_ref));
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu_stack().isApprox(s_new_coarse_ref.mu_stack()));
  EXPECT_TRUE(s_new_coarse.f_stack().isApprox(s_new_coarse_ref.f_stack()));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  // backward correction serial 
  parnmpc.backwardCorrectionSerial(robot, s_next, s_next_new, s_new_coarse);
  // backward correction serial ref
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(dimx);
  x_res.head(robot.dimv()) = s_next_new.lmd - s_next.lmd;
  x_res.tail(robot.dimv()) = s_next_new.gmm - s_next.gmm;
  Eigen::VectorXd dx = kkt_matrix_inv.topRightCorner(dimx, dimx) * x_res;
  s_new_coarse_ref.lmd -= dx.head(robot.dimv());
  s_new_coarse_ref.gmm -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu_stack().isApprox(s_new_coarse_ref.mu_stack()));
  EXPECT_TRUE(s_new_coarse.f_stack().isApprox(s_new_coarse_ref.f_stack()));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  // backward correction parallel
  parnmpc.backwardCorrectionParallel(robot, d, s_new_coarse);
  // backward correction parallel ref
  d_ref.split_direction().tail(dimKKT-dimx)
      = kkt_matrix_inv.bottomRightCorner(dimKKT-dimx, dimx) * x_res;
  s_new_coarse_ref.mu_stack().noalias() -= d_ref.dmu();
  s_new_coarse_ref.f_stack().noalias() -= d_ref.df();
  robot.integrateConfiguration(d_ref.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v.noalias() -= d_ref.dv();
  // forward correction serial
  parnmpc.forwardCorrectionSerial(robot, s_prev, s_prev_new, s_new_coarse);
  // forward correction serial ref
  robot.subtractConfiguration(s_prev_new.q, s_prev.q, x_res.head(robot.dimv()));
  x_res.tail(robot.dimv()) = s_prev_new.v - s_prev.v;
  dx = kkt_matrix_inv.bottomLeftCorner(dimx, dimx) * x_res;
  robot.integrateConfiguration(dx.head(robot.dimv()), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu_stack().isApprox(s_new_coarse_ref.mu_stack()));
  EXPECT_TRUE(s_new_coarse.f_stack().isApprox(s_new_coarse_ref.f_stack()));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  // forward correction parallel
  parnmpc.forwardCorrectionParallel(robot, d, s_new_coarse);
  // forward correction parallel ref
  d_ref.split_direction().head(dimKKT-dimx) = kkt_matrix_inv.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_coarse_ref.lmd -= d_ref.dlmd();
  s_new_coarse_ref.gmm -= d_ref.dgmm();
  s_new_coarse_ref.mu_stack() -= d_ref.dmu();
  s_new_coarse_ref.f_stack() -= d_ref.df();
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu_stack().isApprox(s_new_coarse_ref.mu_stack()));
  EXPECT_TRUE(s_new_coarse.f_stack().isApprox(s_new_coarse_ref.f_stack()));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  // compute direction
  parnmpc.computePrimalAndDualDirection(robot, s, s_new_coarse, d);
  // compute direction ref
  d_ref.dlmd() = s_new_coarse_ref.lmd - s.lmd;
  d_ref.dgmm() = s_new_coarse_ref.gmm - s.gmm;
  d_ref.dmu() = s_new_coarse_ref.mu_stack() - s.mu_stack();
  d_ref.df() = s_new_coarse_ref.f_stack() - s.f_stack();
  robot.subtractConfiguration(s_new_coarse_ref.q, s.q, d_ref.dq());
  d_ref.dv() = s_new_coarse_ref.v - s.v;
  impulse_dynamics.computeCondensedDirection(kkt_matrix, kkt_residual, d_ref);
  EXPECT_TRUE(d.dlmd().isApprox(d_ref.dlmd()));
  EXPECT_TRUE(d.dgmm().isApprox(d_ref.dgmm()));
  EXPECT_TRUE(d.dmu().isApprox(d_ref.dmu()));
  EXPECT_TRUE(d.df().isApprox(d_ref.df()));
  EXPECT_TRUE(d.dq().isApprox(d_ref.dq()));
  EXPECT_TRUE(d.dv().isApprox(d_ref.dv()));
  EXPECT_TRUE(d.ddv.isApprox(d_ref.ddv));
  EXPECT_TRUE(d.dbeta.isApprox(d_ref.dbeta));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
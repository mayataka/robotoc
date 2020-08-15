#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/ocp/kkt_composition.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"
#include "idocp/ocp/inverse_dynamics_condenser.hpp"
#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"


namespace idocp {

class FixedBaseSplitParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    std::vector<int> contact_frames = {18};
    const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
    const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
    robot = Robot(urdf, contact_frames, baum_a, baum_b);
    std::random_device rnd;
    std::vector<bool> contact_status = {rnd()%2==0};
    robot.setContactStatus(contact_status);
    s = SplitSolution(robot);
    s.set(robot);
    robot.generateFeasibleConfiguration(s.q);
    s.v = Eigen::VectorXd::Random(robot.dimv());
    s.a = Eigen::VectorXd::Random(robot.dimv());
    s.f = Eigen::VectorXd::Random(robot.max_dimf());
    s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s.lmd = Eigen::VectorXd::Random(robot.dimv());
    s.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_tmp = SplitSolution(robot);
    s_old = SplitSolution(robot);
    s_old.set(robot);
    robot.generateFeasibleConfiguration(s_old.q);
    s_old.v = Eigen::VectorXd::Random(robot.dimv());
    s_old.a = Eigen::VectorXd::Random(robot.dimv());
    s_old.f = Eigen::VectorXd::Random(robot.max_dimf());
    s_old.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s_old.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_old.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_new = SplitSolution(robot);
    s_new.set(robot);
    robot.generateFeasibleConfiguration(s_new.q);
    s_new.v = Eigen::VectorXd::Random(robot.dimv());
    s_new.a = Eigen::VectorXd::Random(robot.dimv());
    s_new.f = Eigen::VectorXd::Random(robot.max_dimf());
    s_new.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s_new.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_new.gmm = Eigen::VectorXd::Random(robot.dimv());
    d = SplitDirection(robot);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    q_prev = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_prev);
    v_prev = Eigen::VectorXd::Random(robot.dimv());
    lmd_next = Eigen::VectorXd::Random(robot.dimv());
    gmm_next = Eigen::VectorXd::Random(robot.dimv());
    q_next = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_next);
    std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
    std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
    const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_ref);
    const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf()).array().abs();
    const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(robot.max_dimf());
    joint_cost->set_q_weight(q_weight);
    joint_cost->set_q_ref(q_ref);
    joint_cost->set_v_weight(v_weight);
    joint_cost->set_v_ref(v_ref);
    joint_cost->set_a_weight(a_weight);
    joint_cost->set_a_ref(a_ref);
    joint_cost->set_u_weight(u_weight);
    joint_cost->set_u_ref(u_ref);
    joint_cost->set_qf_weight(qf_weight);
    joint_cost->set_vf_weight(vf_weight);
    contact_cost->set_f_weight(f_weight);
    contact_cost->set_f_ref(f_ref);
    cost = std::make_shared<CostFunction>();
    cost->push_back(joint_cost);
    cost->push_back(contact_cost);
    cost_data = CostFunctionData(robot);
    constraints = std::make_shared<Constraints>();
    std::shared_ptr<JointPositionLowerLimit> joint_lower_limit 
        = std::make_shared<JointPositionLowerLimit>(robot);
    std::shared_ptr<JointPositionUpperLimit> joint_upper_limit 
        = std::make_shared<JointPositionUpperLimit>(robot);
    constraints->push_back(joint_upper_limit);
    constraints->push_back(joint_lower_limit);
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  std::shared_ptr<CostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<Constraints> constraints;
  SplitSolution s, s_tmp, s_old, s_new;
  SplitDirection d;
  Eigen::VectorXd q_prev, v_prev, lmd_next, gmm_next, q_next;
};


TEST_F(FixedBaseSplitParNMPCTest, initconstraints) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  if (parnmpc.isFeasible(robot, s)) {
    parnmpc.initConstraints(robot, 2, dtau, s);
    EXPECT_TRUE(parnmpc.isFeasible(robot, s));
  }
}

TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNorm) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  pdipm::JointSpaceConstraints constraints(robot);
  constraints.setSlackAndDual(dtau, s.q, s.v, s.a, s.u);
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  ParNMPCLinearizer linearizer(robot);
  InverseDynamicsCondenser condenser(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  linearizer.linearizeStageCost(robot, cost, cost_data, t, dtau, s, kkt_residual);
  linearizer.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, lmd_next,
                                    gmm_next, q_next, kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau, kkt_residual, kkt_matrix);
  constraints.augmentDualResidual(dtau, kkt_residual.lq(), kkt_residual.lv(), 
                                  kkt_residual.la());
  condenser.setContactStatus(robot);
  condenser.linearizeStageCost(robot, cost, cost_data, t, dtau, s);
  condenser.linearizeInequalityConstraints(robot, constraints, t, dtau, s);
  condenser.linearizeInverseDynamics(robot, dtau, s);
  condenser.augmentInverseDynamicsDerivatives(robot, dtau, s, kkt_residual);
  condenser.augmentInverseDynamicsDerivatives(dtau, s);
  condenser.linearizeFloatingBaseConstraint(dtau, s);
  condenser.augmentFloatingBaseConstraint(dtau, s, kkt_residual);
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_active();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_active();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual.lf().noalias() += kkt_matrix.Cf().transpose() * s.mu_active();
  }
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  kkt_error_ref += condenser.squaredKKTErrorNorm();
  kkt_error_ref += constraints.residualSquaredNrom(dtau, s.q, s.v, s.a, s.u);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNormTerminal) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s);
  pdipm::JointSpaceConstraints constraints(robot);
  constraints.setSlackAndDual(dtau, s.q, s.v, s.a, s.u);
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  ParNMPCLinearizer linearizer(robot);
  InverseDynamicsCondenser condenser(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  linearizer.linearizeStageCost(robot, cost, cost_data, t, dtau, s, kkt_residual);
  linearizer.linearizeStateEquation(robot, dtau, q_prev, v_prev, s,
                                    kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau, kkt_residual, kkt_matrix);
  constraints.augmentDualResidual(dtau, kkt_residual.lq(), kkt_residual.lv(), 
                                  kkt_residual.la());
  condenser.setContactStatus(robot);
  condenser.linearizeStageCost(robot, cost, cost_data, t, dtau, s);
  condenser.linearizeInequalityConstraints(robot, constraints, t, dtau, s);
  condenser.linearizeInverseDynamics(robot, dtau, s);
  condenser.augmentInverseDynamicsDerivatives(robot, dtau, s, kkt_residual);
  condenser.augmentInverseDynamicsDerivatives(dtau, s);
  condenser.linearizeFloatingBaseConstraint(dtau, s);
  condenser.augmentFloatingBaseConstraint(dtau, s, kkt_residual);
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_active();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_active();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual.lf().noalias() += kkt_matrix.Cf().transpose() * s.mu_active();
  }
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  kkt_error_ref += condenser.squaredKKTErrorNorm();
  kkt_error_ref += constraints.residualSquaredNrom(dtau, s.q, s.v, s.a, s.u);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, coarseUpdate) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  Eigen::MatrixXd aux_mat_seed = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed.transpose() * aux_mat_seed;
  SplitSolution s_new_coarse(robot);
  s_new_coarse.set(robot);
  parnmpc.coarseUpdate(robot, t, dtau, q_prev, v_prev, s, lmd_next, gmm_next, 
                       q_next, aux_mat_next, d, s_new_coarse);
  pdipm::JointSpaceConstraints constraints(robot);
  constraints.setSlackAndDual(dtau, s.q, s.v, s.a, s.u);
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  ParNMPCLinearizer linearizer(robot);
  InverseDynamicsCondenser condenser(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  linearizer.linearizeStageCost(robot, cost, cost_data, t, dtau, s, kkt_residual);
  linearizer.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, lmd_next,
                                    gmm_next, q_next, kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau, kkt_residual, kkt_matrix);
  constraints.augmentDualResidual(dtau, kkt_residual.lq(), kkt_residual.lv(), 
                                  kkt_residual.la());
  condenser.setContactStatus(robot);
  condenser.linearizeStageCost(robot, cost, cost_data, t, dtau, s);
  condenser.linearizeInequalityConstraints(robot, constraints, t, dtau, s);
  condenser.condenseInequalityConstraints(robot, constraints, t, dtau, s);
  condenser.linearizeInverseDynamics(robot, dtau, s);
  condenser.augmentInverseDynamicsDerivatives(robot, dtau, s, kkt_residual);
  condenser.augmentInverseDynamicsDerivatives(dtau, s);
  condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
  condenser.linearizeFloatingBaseConstraint(dtau, s);
  condenser.condenseFloatingBaseConstraint(dtau, s, kkt_residual, kkt_matrix);
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_active();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_active();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual.lf().noalias() += kkt_matrix.Cf().transpose() * s.mu_active();
  }
  constraints.condenseSlackAndDual(dtau, s.q, s.v, s.a, kkt_matrix.Qqq(), 
                                   kkt_matrix.Qvv(), kkt_matrix.Qaa(), 
                                   kkt_residual.lq(), kkt_residual.lv(),
                                   kkt_residual.la());
  cost->augment_lqq(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qqq());
  cost->augment_lvv(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qvv());
  cost->augment_laa(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qaa());
  cost->augment_lff(robot, cost_data, t, dtau, s.f, kkt_matrix.Qff());
  kkt_matrix.symmetrize();
  kkt_matrix.Qxx().noalias() += aux_mat_next;
  const int dim_kkt = kkt_matrix.dimKKT();
  Eigen::MatrixXd kkt_matrix_inverse = Eigen::MatrixXd::Zero(dim_kkt, dim_kkt);
  kkt_matrix.invert(kkt_matrix_inverse);
  // coarse update of the solution
  d.split_direction() = kkt_matrix_inverse * kkt_residual.KKT_residual();
  SplitSolution s_new_coarse_ref(robot);
  s_new_coarse_ref.set(robot);
  s_new_coarse_ref.lmd = s.lmd - d.dlmd();
  s_new_coarse_ref.gmm = s.gmm - d.dgmm();
  s_new_coarse_ref.mu_active() = s.mu_active() - d.dmu();
  s_new_coarse_ref.a = s.a - d.da();
  s_new_coarse_ref.f_active() = s.f_active() - d.df();
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v = s.v - d.dv();
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  Eigen::MatrixXd aux_mat_ref = Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv());
  parnmpc.getAuxiliaryMatrix(aux_mat_ref);
  EXPECT_TRUE(aux_mat_ref.isApprox(-1*kkt_matrix_inverse.topLeftCorner(2*robot.dimv(), 2*robot.dimv())));
  parnmpc.backwardCollectionSerial(robot, s_old, s_new, s_new_coarse);
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(2*robot.dimv());
  x_res.head(robot.dimv()) = s_new.lmd - s_old.lmd;
  x_res.tail(robot.dimv()) = s_new.gmm - s_old.gmm;
  Eigen::VectorXd dx = kkt_matrix_inverse.topRightCorner(2*robot.dimv(), 2*robot.dimv()) * x_res;
  s_new_coarse_ref.lmd -= dx.head(robot.dimv());
  s_new_coarse_ref.gmm -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  parnmpc.backwardCollectionParallel(robot, d, s_new_coarse);
  const int dimx = 2*robot.dimv();
  d.split_direction().segment(dimx, dim_kkt-dimx) 
      = kkt_matrix_inverse.bottomRightCorner(dim_kkt-dimx, dimx) * x_res;
  s_new_coarse_ref.mu_active().noalias() -= d.dmu();
  s_new_coarse_ref.a.noalias() -= d.da();
  s_new_coarse_ref.f_active().noalias() -= d.df();
  robot.integrateConfiguration(d.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v.noalias() -= d.dv();
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  parnmpc.forwardCollectionSerial(robot, s_old, s_new, s_new_coarse);
  robot.subtractConfiguration(s_new.q, s_old.q, x_res.head(robot.dimv()));
  x_res.tail(robot.dimv()) = s_new.v - s_old.v;
  dx = kkt_matrix_inverse.bottomLeftCorner(2*robot.dimv(), 2*robot.dimv()) * x_res;
  robot.integrateConfiguration(dx.head(robot.dimv()), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  parnmpc.forwardCollectionParallel(robot, d, s_new_coarse);
  d.split_direction().segment(0, dim_kkt-dimx) 
      = kkt_matrix_inverse.topLeftCorner(dim_kkt-dimx, dimx) * x_res;
  s_new_coarse_ref.lmd -= d.dlmd();
  s_new_coarse_ref.gmm -= d.dgmm();
  s_new_coarse_ref.mu_active() -= d.dmu();
  s_new_coarse_ref.a -= d.da();
  s_new_coarse_ref.f_active() -= d.df();
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
}


TEST_F(FixedBaseSplitParNMPCTest, coarseUpdateTerminal) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  Eigen::MatrixXd aux_mat_seed = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed.transpose() * aux_mat_seed;
  SplitSolution s_new_coarse(robot);
  s_new_coarse.set(robot);
  parnmpc.coarseUpdate(robot, t, dtau, q_prev, v_prev, s, d, s_new_coarse);
  pdipm::JointSpaceConstraints constraints(robot);
  constraints.setSlackAndDual(dtau, s.q, s.v, s.a, s.u);
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  ParNMPCLinearizer linearizer(robot);
  InverseDynamicsCondenser condenser(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  linearizer.linearizeStageCost(robot, cost, cost_data, t, dtau, s, kkt_residual);
  linearizer.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau, kkt_residual, kkt_matrix);
  constraints.augmentDualResidual(dtau, kkt_residual.lq(), kkt_residual.lv(), 
                                  kkt_residual.la());
  condenser.setContactStatus(robot);
  condenser.linearizeStageCost(robot, cost, cost_data, t, dtau, s);
  condenser.linearizeInequalityConstraints(robot, constraints, t, dtau, s);
  condenser.condenseInequalityConstraints(robot, constraints, t, dtau, s);
  condenser.linearizeInverseDynamics(robot, dtau, s);
  condenser.augmentInverseDynamicsDerivatives(robot, dtau, s, kkt_residual);
  condenser.augmentInverseDynamicsDerivatives(dtau, s);
  condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
  condenser.linearizeFloatingBaseConstraint(dtau, s);
  condenser.condenseFloatingBaseConstraint(dtau, s, kkt_residual, kkt_matrix);
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_active();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_active();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual.lf().noalias() += kkt_matrix.Cf().transpose() * s.mu_active();
  }
  constraints.condenseSlackAndDual(dtau, s.q, s.v, s.a, kkt_matrix.Qqq(), 
                                   kkt_matrix.Qvv(), kkt_matrix.Qaa(), 
                                   kkt_residual.lq(), kkt_residual.lv(),
                                   kkt_residual.la());
  cost->augment_lqq(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qqq());
  cost->augment_lvv(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qvv());
  cost->augment_laa(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qaa());
  cost->augment_lff(robot, cost_data, t, dtau, s.f, kkt_matrix.Qff());
  cost->augment_phiqq(robot, cost_data, t, s.q, s.v, kkt_matrix.Qqq());
  cost->augment_phivv(robot, cost_data, t, s.q, s.v, kkt_matrix.Qvv());
  kkt_matrix.symmetrize();
  const int dim_kkt = kkt_matrix.dimKKT();
  Eigen::MatrixXd kkt_matrix_inverse = Eigen::MatrixXd::Zero(dim_kkt, dim_kkt);
  kkt_matrix.invert(kkt_matrix_inverse);
  // coarse update of the solution
  d.split_direction() = kkt_matrix_inverse * kkt_residual.KKT_residual();
  SplitSolution s_new_coarse_ref(robot);
  s_new_coarse_ref.set(robot);
  s_new_coarse_ref.lmd = s.lmd - d.dlmd();
  s_new_coarse_ref.gmm = s.gmm - d.dgmm();
  s_new_coarse_ref.mu_active() = s.mu_active() - d.dmu();
  s_new_coarse_ref.a = s.a - d.da();
  s_new_coarse_ref.f_active() = s.f_active() - d.df();
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v = s.v - d.dv();
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  Eigen::MatrixXd aux_mat_ref = Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv());
  parnmpc.getAuxiliaryMatrix(aux_mat_ref);
  EXPECT_TRUE(aux_mat_ref.isApprox(-1*kkt_matrix_inverse.topLeftCorner(2*robot.dimv(), 2*robot.dimv())));
  parnmpc.backwardCollectionSerial(robot, s_old, s_new, s_new_coarse);
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(2*robot.dimv());
  x_res.head(robot.dimv()) = s_new.lmd - s_old.lmd;
  x_res.tail(robot.dimv()) = s_new.gmm - s_old.gmm;
  Eigen::VectorXd dx = kkt_matrix_inverse.topRightCorner(2*robot.dimv(), 2*robot.dimv()) * x_res;
  s_new_coarse_ref.lmd -= dx.head(robot.dimv());
  s_new_coarse_ref.gmm -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  parnmpc.backwardCollectionParallel(robot, d, s_new_coarse);
  const int dimx = 2*robot.dimv();
  d.split_direction().segment(dimx, dim_kkt-dimx) 
      = kkt_matrix_inverse.bottomRightCorner(dim_kkt-dimx, dimx) * x_res;
  s_new_coarse_ref.mu_active().noalias() -= d.dmu();
  s_new_coarse_ref.a.noalias() -= d.da();
  s_new_coarse_ref.f_active().noalias() -= d.df();
  robot.integrateConfiguration(d.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v.noalias() -= d.dv();
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  parnmpc.forwardCollectionSerial(robot, s_old, s_new, s_new_coarse);
  robot.subtractConfiguration(s_new.q, s_old.q, x_res.head(robot.dimv()));
  x_res.tail(robot.dimv()) = s_new.v - s_old.v;
  dx = kkt_matrix_inverse.bottomLeftCorner(2*robot.dimv(), 2*robot.dimv()) * x_res;
  robot.integrateConfiguration(dx.head(robot.dimv()), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  parnmpc.forwardCollectionParallel(robot, d, s_new_coarse);
  d.split_direction().segment(0, dim_kkt-dimx) 
      = kkt_matrix_inverse.topLeftCorner(dim_kkt-dimx, dimx) * x_res;
  s_new_coarse_ref.lmd -= d.dlmd();
  s_new_coarse_ref.gmm -= d.dgmm();
  s_new_coarse_ref.mu_active() -= d.dmu();
  s_new_coarse_ref.a -= d.da();
  s_new_coarse_ref.f_active() -= d.df();
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
}

TEST_F(FixedBaseSplitParNMPCTest, computePrimalDualDirection) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  parnmpc.computePrimalAndDualDirection(robot, dtau, s, s_new, d);
  EXPECT_TRUE(d.dlmd().isApprox(s_new.lmd-s.lmd));
  EXPECT_TRUE(d.dgmm().isApprox(s_new.gmm-s.gmm));
  EXPECT_TRUE(d.dmu().isApprox(s_new.mu_active()-s.mu_active()));
  EXPECT_TRUE(d.da().isApprox(s_new.a-s.a));
  EXPECT_TRUE(d.df().isApprox(s_new.f_active()-s.f_active()));
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s_new.q, s.q, qdiff);
  EXPECT_TRUE(d.dq().isApprox(qdiff));
  EXPECT_TRUE(d.dv().isApprox(s_new.v-s.v));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
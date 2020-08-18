#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"


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
    s.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s.q);
    s.v = Eigen::VectorXd::Random(robot.dimv());
    s.a = Eigen::VectorXd::Random(robot.dimv());
    s.f = Eigen::VectorXd::Random(robot.max_dimf());
    s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s.lmd = Eigen::VectorXd::Random(robot.dimv());
    s.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_tmp = SplitSolution(robot);
    s_old = SplitSolution(robot);
    s_old.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s_old.q);
    s_old.v = Eigen::VectorXd::Random(robot.dimv());
    s_old.a = Eigen::VectorXd::Random(robot.dimv());
    s_old.f = Eigen::VectorXd::Random(robot.max_dimf());
    s_old.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s_old.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_old.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_new = SplitSolution(robot);
    s_new.setContactStatus(robot);
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
    auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
    auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
    auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
    auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
    constraints->push_back(joint_upper_limit); 
    constraints->push_back(joint_lower_limit);
    constraints->push_back(velocity_lower_limit); 
    constraints->push_back(velocity_upper_limit);
    constraints_data = constraints->createConstraintsData(robot);
    kkt_matrix = KKTMatrix(robot);
    kkt_residual = KKTResidual(robot);
    state_equation = StateEquation(robot);
    inverse_dynamics = InverseDynamics(robot);
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  std::shared_ptr<CostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<Constraints> constraints;
  ConstraintsData constraints_data;
  SplitSolution s, s_tmp, s_old, s_new;
  SplitDirection d, d_prev;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  StateEquation state_equation;
  InverseDynamics inverse_dynamics;
  Eigen::VectorXd q_prev, v_prev, lmd_next, gmm_next, q_next;
};


TEST_F(FixedBaseSplitParNMPCTest, isFeasible) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  EXPECT_EQ(parnmpc.isFeasible(robot, s), 
            constraints->isFeasible(robot, constraints_data, s));
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNorm) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  kkt_error_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNormEmptyCostAndEmptyConstraints) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<CostFunction>();
  SplitParNMPC parnmpc(robot, empty_cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<Constraints>();
  SplitParNMPC parnmpc(robot, cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNormOnlyStateEquation) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false}));
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  robot.RNEA(s.q, s.v, s.a, s.u);
  s.beta.setZero();
  s.lmd.setZero();
  s.gmm.setZero();
  lmd_next.setZero();
  gmm_next.setZero();
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()+kkt_residual.Fv().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNormTerminal) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error 
      = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, q_prev, v_prev, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  state_equation.linearizeStateEquationTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  kkt_error_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitParNMPCTest, KKTErrorNormTerminalEmptyCostAndEmptyConstraints) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error 
      = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, q_prev, v_prev, s);
  robot.updateKinematics(s.q, s.v, s.a);
  state_equation.linearizeStateEquationTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


// TEST_F(FixedBaseSplitParNMPCTest, coarseUpdate) {
//   SplitParNMPC parnmpc(robot, cost, constraints);
//   parnmpc.initConstraints(robot, 2, dtau, s);
//   constraints->setSlackAndDual(robot, constraints_data, dtau, 
//                                s.a, s.f, s.q, s.v, s.u);
//   Eigen::MatrixXd aux_mat_seed = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
//   const Eigen::MatrixXd aux_mat_next = aux_mat_seed.transpose() * aux_mat_seed;
//   SplitSolution s_new_coarse(robot);
//   s_new_coarse.setContactStatus(robot);
//   parnmpc.coarseUpdate(robot, t, dtau, q_prev, v_prev, s, lmd_next, gmm_next, 
//                        q_next, aux_mat_next, d, s_new_coarse);
//   KKTResidual kkt_residual(robot);
//   KKTMatrix kkt_matrix(robot);
//   ParNMPCLinearizer linearizer(robot);
//   InverseDynamicsCondenser condenser(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   linearizer.linearizeCostAndConstraints(robot, cost, cost_data, constraints, 
//                                          constraints_data, t, dtau, s, kkt_residual);
//   linearizer.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, lmd_next,
//                                     gmm_next, q_next, kkt_residual, kkt_matrix);
//   linearizer.linearizeContactConstraints(robot, dtau, s, kkt_residual, kkt_matrix);
//   condenser.setContactStatus(robot);
//   condenser.linearizeCostAndConstraints(robot, cost, cost_data, constraints,
//                                         constraints_data, t, dtau, s);
//   condenser.linearizeInverseDynamics(robot, dtau, s);
//   condenser.linearizeFloatingBaseConstraint(dtau, s, kkt_residual);
//   condenser.augmentInverseDynamicsDerivatives(dtau, s, kkt_residual);
//   condenser.condenseInequalityConstraints(robot, constraints, constraints_data, 
//                                           t, dtau, s);
//   condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
//   condenser.condenseFloatingBaseConstraint(dtau, kkt_residual, kkt_matrix);
//   linearizer.condenseInequalityConstraints(robot, constraints, constraints_data,
//                                            t, dtau, s, kkt_residual, kkt_matrix);
//   cost->augment_laa(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qaa());
//   if (robot.dimf() > 0) {
//     cost->augment_lff(robot, cost_data, t, dtau, s.f, kkt_matrix.Qff());
//   }
//   cost->augment_lqq(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qqq());
//   cost->augment_lvv(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qvv());
//   kkt_matrix.Qxx().noalias() += aux_mat_next;
//   kkt_matrix.symmetrize();
//   KKTMatrixInverse kkt_matrix_inverse(robot);
//   kkt_matrix_inverse.setContactStatus(robot);
//   kkt_matrix.invert(kkt_matrix_inverse.KKT_matrix_inverse());
//   d.split_direction() 
//       = kkt_matrix_inverse.KKT_matrix_inverse() * kkt_residual.KKT_residual();
//   std::cout << kkt_matrix.KKT_matrix() << std::endl;
//   std::cout << std::endl;
//   std::cout << kkt_matrix_inverse.KKT_matrix_inverse() << std::endl;
//   EXPECT_TRUE((kkt_matrix.KKT_matrix()*kkt_matrix_inverse.KKT_matrix_inverse())
//               .isApprox(Eigen::MatrixXd::Identity(kkt_matrix.dimKKT(), kkt_matrix.dimKKT())));
//   SplitSolution s_new_coarse_ref(robot);
//   s_new_coarse_ref.setContactStatus(robot);
//   s_new_coarse_ref.lmd = s.lmd - d.dlmd();
//   s_new_coarse_ref.gmm = s.gmm - d.dgmm();
//   s_new_coarse_ref.mu_active() = s.mu_active() - d.dmu();
//   s_new_coarse_ref.a = s.a - d.da();
//   s_new_coarse_ref.f_active() = s.f_active() - d.df();
//   robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse_ref.q);
//   s_new_coarse_ref.v = s.v - d.dv();
//   EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
//   EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
//   EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
//   EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
//   EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
//   Eigen::MatrixXd aux_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv());
//   parnmpc.getAuxiliaryMatrix(aux_mat);
//   Eigen::MatrixXd aux_mat_ref = kkt_matrix_inverse.auxiliaryMatrix();
//   EXPECT_TRUE(aux_mat.isApprox(-1*aux_mat_ref));
//   parnmpc.backwardCorrectionSerial(robot, s_old, s_new, s_new_coarse);
//   Eigen::VectorXd x_res = Eigen::VectorXd::Zero(2*robot.dimv());
//   x_res.head(robot.dimv()) = s_new.lmd - s_old.lmd;
//   x_res.tail(robot.dimv()) = s_new.gmm - s_old.gmm;
//   Eigen::VectorXd dx = kkt_matrix_inverse.backwardCorrectionSerialCoeff() * x_res;
//   s_new_coarse_ref.lmd -= dx.head(robot.dimv());
//   s_new_coarse_ref.gmm -= dx.tail(robot.dimv());
//   EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
//   EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
//   parnmpc.backwardCorrectionParallel(robot, d, s_new_coarse);
//   const int dimx = 2*robot.dimv();
//   d.backwardCorrectionParallelDirection()
//       = kkt_matrix_inverse.backwardCorrectionParallelCoeff()* x_res;
//   s_new_coarse_ref.mu_active().noalias() -= d.dmu();
//   s_new_coarse_ref.a.noalias() -= d.da();
//   s_new_coarse_ref.f_active().noalias() -= d.df();
//   robot.integrateConfiguration(d.dq(), -1, s_new_coarse_ref.q);
//   s_new_coarse_ref.v.noalias() -= d.dv();
//   EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
//   EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
//   EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
//   parnmpc.forwardCorrectionSerial(robot, s_old, s_new, s_new_coarse);
//   robot.subtractConfiguration(s_new.q, s_old.q, x_res.head(robot.dimv()));
//   x_res.tail(robot.dimv()) = s_new.v - s_old.v;
//   dx = kkt_matrix_inverse.forwardCorrectionSerialCoeff() * x_res;
//   robot.integrateConfiguration(dx.head(robot.dimv()), -1, s_new_coarse_ref.q);
//   s_new_coarse_ref.v -= dx.tail(robot.dimv());
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
//   parnmpc.forwardCorrectionParallel(robot, d, s_new_coarse);
//   d.forwardCorrectionParallelDirection()
//       = kkt_matrix_inverse.forwardCorrectionParallelCoeff() * x_res;
//   s_new_coarse_ref.lmd -= d.dlmd();
//   s_new_coarse_ref.gmm -= d.dgmm();
//   s_new_coarse_ref.mu_active() -= d.dmu();
//   s_new_coarse_ref.a -= d.da();
//   s_new_coarse_ref.f_active() -= d.df();
//   EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
//   EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
//   EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
//   EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
//   EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
// }


// TEST_F(FixedBaseSplitParNMPCTest, coarseUpdateTerminal) {
//   SplitParNMPC parnmpc(robot, cost, constraints);
//   parnmpc.initConstraints(robot, 2, dtau, s);
//   constraints->setSlackAndDual(robot, constraints_data, dtau, 
//                                s.a, s.f, s.q, s.v, s.u);
//   Eigen::MatrixXd aux_mat_seed = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
//   const Eigen::MatrixXd aux_mat_next = aux_mat_seed.transpose() * aux_mat_seed;
//   SplitSolution s_new_coarse(robot);
//   s_new_coarse.setContactStatus(robot);
//   parnmpc.coarseUpdate(robot, t, dtau, q_prev, v_prev, s, d, s_new_coarse);
//   KKTResidual kkt_residual(robot);
//   KKTMatrix kkt_matrix(robot);
//   ParNMPCLinearizer linearizer(robot);
//   InverseDynamicsCondenser condenser(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   linearizer.linearizeCostAndConstraints(robot, cost, cost_data, constraints, 
//                                          constraints_data, t, dtau, s, kkt_residual);
//   linearizer.linearizeStateEquation(robot, dtau, q_prev, v_prev, s,
//                                     kkt_residual, kkt_matrix);
//   linearizer.linearizeContactConstraints(robot, dtau, s, kkt_residual, kkt_matrix);
//   condenser.setContactStatus(robot);
//   condenser.linearizeCostAndConstraints(robot, cost, cost_data, constraints,
//                                         constraints_data, t, dtau, s);
//   condenser.linearizeInverseDynamics(robot, dtau, s);
//   condenser.linearizeFloatingBaseConstraint(dtau, s, kkt_residual);
//   condenser.augmentInverseDynamicsDerivatives(dtau, s, kkt_residual);
//   condenser.condenseInequalityConstraints(robot, constraints, constraints_data, 
//                                           t, dtau, s);
//   condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
//   condenser.condenseFloatingBaseConstraint(dtau, kkt_residual, kkt_matrix);
//   linearizer.condenseInequalityConstraints(robot, constraints, constraints_data,
//                                            t, dtau, s, kkt_residual, kkt_matrix);
//   cost->augment_laa(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qaa());
//   if (robot.dimf() > 0) {
//     cost->augment_lff(robot, cost_data, t, dtau, s.f, kkt_matrix.Qff());
//   }
//   cost->augment_lqq(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qqq());
//   cost->augment_lvv(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_matrix.Qvv());
//   cost->augment_phiqq(robot, cost_data, t, s.q, s.v, kkt_matrix.Qqq());
//   cost->augment_phivv(robot, cost_data, t, s.q, s.v, kkt_matrix.Qvv());
//   kkt_matrix.symmetrize();
//   KKTMatrixInverse kkt_matrix_inverse(robot);
//   kkt_matrix_inverse.setContactStatus(robot);
//   kkt_matrix.invert(kkt_matrix_inverse.KKT_matrix_inverse());
//   d.split_direction() 
//       = kkt_matrix_inverse.KKT_matrix_inverse() * kkt_residual.KKT_residual();
//   SplitSolution s_new_coarse_ref(robot);
//   s_new_coarse_ref.setContactStatus(robot);
//   s_new_coarse_ref.lmd = s.lmd - d.dlmd();
//   s_new_coarse_ref.gmm = s.gmm - d.dgmm();
//   s_new_coarse_ref.mu_active() = s.mu_active() - d.dmu();
//   s_new_coarse_ref.a = s.a - d.da();
//   s_new_coarse_ref.f_active() = s.f_active() - d.df();
//   robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse_ref.q);
//   s_new_coarse_ref.v = s.v - d.dv();
//   EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
//   EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
//   EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
//   EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
//   EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
//   Eigen::MatrixXd aux_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv());
//   parnmpc.getAuxiliaryMatrix(aux_mat);
//   Eigen::MatrixXd aux_mat_ref = kkt_matrix_inverse.auxiliaryMatrix();
//   EXPECT_TRUE(aux_mat.isApprox(-1*aux_mat_ref));
//   parnmpc.backwardCorrectionSerial(robot, s_old, s_new, s_new_coarse);
//   Eigen::VectorXd x_res = Eigen::VectorXd::Zero(2*robot.dimv());
//   x_res.head(robot.dimv()) = s_new.lmd - s_old.lmd;
//   x_res.tail(robot.dimv()) = s_new.gmm - s_old.gmm;
//   Eigen::VectorXd dx = kkt_matrix_inverse.backwardCorrectionSerialCoeff() * x_res;
//   s_new_coarse_ref.lmd -= dx.head(robot.dimv());
//   s_new_coarse_ref.gmm -= dx.tail(robot.dimv());
//   EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
//   EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
//   parnmpc.backwardCorrectionParallel(robot, d, s_new_coarse);
//   const int dimx = 2*robot.dimv();
//   d.backwardCorrectionParallelDirection()
//       = kkt_matrix_inverse.backwardCorrectionParallelCoeff()* x_res;
//   s_new_coarse_ref.mu_active().noalias() -= d.dmu();
//   s_new_coarse_ref.a.noalias() -= d.da();
//   s_new_coarse_ref.f_active().noalias() -= d.df();
//   robot.integrateConfiguration(d.dq(), -1, s_new_coarse_ref.q);
//   s_new_coarse_ref.v.noalias() -= d.dv();
//   EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
//   EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
//   EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
//   parnmpc.forwardCorrectionSerial(robot, s_old, s_new, s_new_coarse);
//   robot.subtractConfiguration(s_new.q, s_old.q, x_res.head(robot.dimv()));
//   x_res.tail(robot.dimv()) = s_new.v - s_old.v;
//   dx = kkt_matrix_inverse.forwardCorrectionSerialCoeff() * x_res;
//   robot.integrateConfiguration(dx.head(robot.dimv()), -1, s_new_coarse_ref.q);
//   s_new_coarse_ref.v -= dx.tail(robot.dimv());
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
//   parnmpc.forwardCorrectionParallel(robot, d, s_new_coarse);
//   d.forwardCorrectionParallelDirection()
//       = kkt_matrix_inverse.forwardCorrectionParallelCoeff() * x_res;
//   s_new_coarse_ref.lmd -= d.dlmd();
//   s_new_coarse_ref.gmm -= d.dgmm();
//   s_new_coarse_ref.mu_active() -= d.dmu();
//   s_new_coarse_ref.a -= d.da();
//   s_new_coarse_ref.f_active() -= d.df();
//   EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
//   EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
//   EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
//   EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
//   EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
//   EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
//   EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
// }


// TEST_F(FixedBaseSplitParNMPCTest, computePrimalDualDirection) {
//   SplitParNMPC parnmpc(robot, cost, constraints);
//   parnmpc.initConstraints(robot, 2, dtau, s);
//   parnmpc.computePrimalAndDualDirection(robot, dtau, s, s_new, d);
//   EXPECT_TRUE(d.dlmd().isApprox(s_new.lmd-s.lmd));
//   EXPECT_TRUE(d.dgmm().isApprox(s_new.gmm-s.gmm));
//   EXPECT_TRUE(d.dmu().isApprox(s_new.mu_active()-s.mu_active()));
//   EXPECT_TRUE(d.da().isApprox(s_new.a-s.a));
//   EXPECT_TRUE(d.df().isApprox(s_new.f_active()-s.f_active()));
//   Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
//   robot.subtractConfiguration(s_new.q, s.q, qdiff);
//   EXPECT_TRUE(d.dq().isApprox(qdiff));
//   EXPECT_TRUE(d.dv().isApprox(s_new.v-s.v));
// }

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
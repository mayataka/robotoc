#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"

#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_gain.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"


namespace idocp {

class FixedBaseSplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    std::vector<int> contact_frames = {18};
    std::vector<double> mu;
    for (int i=0; i<contact_frames.size(); ++i) {
      mu.push_back(std::abs(Eigen::VectorXd::Random(1)[0]));
    }
    const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
    const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
    robot = Robot(urdf, contact_frames, mu, baum_a, baum_b);
    std::random_device rnd;
    contact_status.push_back(rnd()%2==0);
    robot.setContactStatus(contact_status);
    q_prev = Eigen::VectorXd::Zero(robot.dimq());
    robot.generateFeasibleConfiguration(q_prev);
    s = SplitSolution::Random(robot);
    s_next = SplitSolution::Random(robot);
    s_tmp = SplitSolution::Random(robot);
    d = SplitDirection::Random(robot);
    d_next = SplitDirection::Random(robot);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    auto joint_cost = std::make_shared<JointSpaceCost>(robot);
    auto contact_cost = std::make_shared<ContactCost>(robot);
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
    const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    std::vector<Eigen::Vector3d> f_weight, f_ref;
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      const Eigen::Vector3d f_weight_positive = Eigen::Vector3d::Random().array().abs();
      f_weight.push_back(f_weight_positive);
      f_ref.push_back(Eigen::Vector3d::Random());
    }
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
    robotdynamics = RobotDynamics(robot);
    gain = RiccatiGain(robot);
    factorizer = RiccatiMatrixFactorizer(robot);
    inverter = RiccatiMatrixInverter(robot);
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  std::vector<bool> contact_status;
  std::shared_ptr<CostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<Constraints> constraints;
  ConstraintsData constraints_data;
  Eigen::VectorXd q_prev;
  SplitSolution s, s_next, s_tmp;
  SplitDirection d, d_next;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  StateEquation state_equation;
  RobotDynamics robotdynamics;
  RiccatiGain gain;
  RiccatiMatrixFactorizer factorizer;
  RiccatiMatrixInverter inverter;
};


TEST_F(FixedBaseSplitOCPTest, isFeasible) {
  s.setContactStatus(robot);
  SplitOCP ocp(robot, cost, constraints);
  EXPECT_EQ(constraints->isFeasible(robot, constraints_data, s),
            ocp.isFeasible(robot, s));
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormOnlyStateEquation) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  robot.RNEA(s.q, s.v, s.a, s.u);
  s.beta.setZero();
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.computeSquaredKKTErrorNorm(robot, t, dtau, q_prev, s, s_next);
  kkt_residual.Fq() = s.q + dtau * s.v - s_next.q;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormStateEquationAndInverseDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.computeSquaredKKTErrorNorm(robot, t, dtau, q_prev, s, s_next);
  kkt_residual.Fq() = s.q + dtau * s.v - s_next.q;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  kkt_residual.lq() += dtau * du_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * du_dv.transpose() * s.beta;
  kkt_residual.la() += dtau * du_da.transpose() * s.beta;
  kkt_residual.lu -= dtau * s.beta;
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la().squaredNorm()
                         + kkt_residual.lu.squaredNorm()
                         + dtau*dtau*kkt_residual.u_res.squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormStateEquationAndRobotDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({true}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.computeSquaredKKTErrorNorm(robot, t, dtau, q_prev, s, s_next);
  kkt_residual.Fq() = s.q + dtau * s.v - s_next.q;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimf());
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  kkt_residual.lq() += dtau * du_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * du_dv.transpose() * s.beta;
  kkt_residual.la() += dtau * du_da.transpose() * s.beta;
  kkt_residual.lf() += dtau * du_df.transpose() * s.beta;
  kkt_residual.lu -= dtau * s.beta;
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                    kkt_matrix.Ca());
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_stack();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_stack();
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la().squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.lu.squaredNorm()
                         + dtau*dtau*kkt_residual.u_res.squaredNorm()
                         + kkt_residual.C().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + dtau * kkt_residual.u_res.lpNorm<1>()
                                + kkt_residual.C().lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<CostFunction>();
  SplitOCP ocp(robot, empty_cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.computeSquaredKKTErrorNorm(robot, t, dtau, q_prev, s, s_next);
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s.u, kkt_residual.lu);
  state_equation.linearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  robotdynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<Constraints>();
  SplitOCP ocp(robot, cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.computeSquaredKKTErrorNorm(robot, t, dtau, q_prev, s, s_next);
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  state_equation.linearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  robotdynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNorm) {
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.computeSquaredKKTErrorNorm(robot, t, dtau, q_prev, s, s_next);
  kkt_matrix.setContactStatus(robot);
  kkt_residual.setContactStatus(robot);
  kkt_residual.setZero();
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, 
                                     kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s.u, kkt_residual.lu);
  state_equation.linearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  robotdynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  kkt_error_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, costAndViolation) {
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.computeSquaredKKTErrorNorm(robot, t, dtau, q_prev, s, s_next);
  const auto pair = ocp.costAndViolation(robot, t, dtau, s); 
  const double cost_ref 
      = cost->l(robot, cost_data, t, dtau, s) 
          + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C().head(robot.dimf()).lpNorm<1>()
          + constraints->residualL1Nrom(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseSplitOCPTest, costAndViolationWithStepSize) {
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double step_size = 0.3;
  const auto pair = ocp.costAndViolation(robot, step_size, t, dtau, s, d, s_next, d_next); 
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.a = s.a + step_size * d.da();
  s_tmp.f_stack() = s.f_stack() + step_size * d.df();
  s_tmp.set_f();
  s_tmp.u = s.u + step_size * d.du;
  const double cost_ref 
      = cost->l(robot, cost_data, t, dtau, s_tmp) 
          + constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s_tmp.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s_tmp.v - step_size * d_next.dq();
  kkt_residual.Fv() = s_tmp.v + dtau * s_tmp.a - s_next.v - step_size * d_next.dv();
  robot.setContactForces(s_tmp.f);
  robot.RNEA(s_tmp.q, s_tmp.v, s_tmp.a, kkt_residual.u_res);
  kkt_residual.u_res -= s_tmp.u;
  robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C().head(robot.dimf()).lpNorm<1>()
          + constraints->residualL1Nrom(robot, constraints_data, dtau, s_tmp);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseSplitOCPTest, riccatiRecursion) {
  const int dimv = robot.dimv();
  const int dimf = robot.dimf();
  const int dimc = robot.dimf() + robot.dim_passive();
  std::cout << "dimf = " << dimf << std::endl;
  std::cout << "dimc = " << dimc << std::endl;
  s.setContactStatus(robot);
  SplitOCP ocp(robot, cost, constraints);
  ASSERT_FALSE(robot.has_floating_base());
  ASSERT_TRUE(robot.dim_passive() == 0);
  while (!ocp.isFeasible(robot, s)) {
    s.v = Eigen::VectorXd::Random(dimv);
    s.u = Eigen::VectorXd::Random(dimv);
  }
  ASSERT_TRUE(ocp.isFeasible(robot, s));
  ocp.initConstraints(robot, 5, dtau, s);
  ocp.linearizeOCP(robot, t, dtau, q_prev, s, s_next);
  RiccatiFactorization riccati_next(robot);
  Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pqq = seed_mat * seed_mat.transpose();
  seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pqv = seed_mat * seed_mat.transpose();
  seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pvq = seed_mat * seed_mat.transpose();
  seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pvv = seed_mat * seed_mat.transpose();
  RiccatiFactorization riccati(robot);
  ocp.backwardRiccatiRecursion(dtau, riccati_next, riccati);
  ocp.forwardRiccatiRecursion(dtau, d, d_next);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->lq(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lv(robot, cost_data, t, dtau, s, kkt_residual);
  cost->la(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lf(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  Eigen::VectorXd q_res = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd v_res = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv, dimf);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.dRNEAPartialdFext(du_df);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                    kkt_matrix.Ca());
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s.u,
                                   kkt_residual.lu);
  cost->luu(robot, cost_data, t, dtau, s.u, kkt_matrix.Quu);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s.u, 
                                    kkt_matrix.Quu, kkt_residual.lu);
  const Eigen::VectorXd lu_condensed = kkt_residual.lu + kkt_matrix.Quu * kkt_residual.u_res;
  kkt_residual.lq() += du_dq.transpose() * lu_condensed;
  kkt_residual.lv() += du_dv.transpose() * lu_condensed;
  kkt_residual.la() += du_da.transpose() * lu_condensed;
  kkt_residual.lf() += du_df.transpose() * lu_condensed;
  kkt_residual.lq() += s_next.lmd - s.lmd;
  kkt_residual.lv() += dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() += dtau * s_next.gmm;
  constraints->augmentDualResidual(robot, constraints_data, dtau, s,
                                   kkt_residual);
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_stack();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_stack();
  kkt_residual.lf() += kkt_matrix.Cf().transpose() * s.mu_stack();
  kkt_matrix.Qqq() = du_dq.transpose() * kkt_matrix.Quu * du_dq;
  kkt_matrix.Qqv() = du_dq.transpose() * kkt_matrix.Quu * du_dv;
  kkt_matrix.Qqa() = du_dq.transpose() * kkt_matrix.Quu * du_da;
  kkt_matrix.Qvv() = du_dv.transpose() * kkt_matrix.Quu * du_dv;
  kkt_matrix.Qva() = du_dv.transpose() * kkt_matrix.Quu * du_da;
  kkt_matrix.Qaa() = du_da.transpose() * kkt_matrix.Quu * du_da;
  kkt_matrix.Qqf() = du_dq.transpose() * kkt_matrix.Quu * du_df;
  kkt_matrix.Qvf() = du_dv.transpose() * kkt_matrix.Quu * du_df;
  kkt_matrix.Qaf() = du_da.transpose() * kkt_matrix.Quu * du_df;
  kkt_matrix.Qff() = du_df.transpose() * kkt_matrix.Quu * du_df;
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s, 
                                    kkt_matrix, kkt_residual);
  cost->computeStageCostHessian(robot, cost_data, t, dtau, s, kkt_matrix);
  inverter.setContactStatus(robot);
  gain.setContactStatus(robot);
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qaq() = kkt_matrix.Qqa().transpose();
  kkt_matrix.Qav() = kkt_matrix.Qva().transpose();
  if (dimf > 0) {
    kkt_matrix.Qfq() = kkt_matrix.Qqf().transpose();
    kkt_matrix.Qfv() = kkt_matrix.Qvf().transpose();
    kkt_matrix.Qfa() = kkt_matrix.Qaf().transpose();
  }
  factorizer.factorizeF(dtau, riccati_next.Pqq, riccati_next.Pqv, 
                        riccati_next.Pvq, riccati_next.Pvv, 
                        kkt_matrix.Qqq(), kkt_matrix.Qqv(), 
                        kkt_matrix.Qvq(), kkt_matrix.Qvv());
  factorizer.factorizeH(dtau, riccati_next.Pqv, riccati_next.Pvv, 
                        kkt_matrix.Qqa(), kkt_matrix.Qva());
  factorizer.factorizeG(dtau, riccati_next.Pvv, kkt_matrix.Qaa());
  kkt_matrix.Qaq() = kkt_matrix.Qqa().transpose();
  kkt_matrix.Qav() = kkt_matrix.Qva().transpose();
  kkt_residual.la() += dtau * riccati_next.Pvq * kkt_residual.Fq();
  kkt_residual.la() += dtau * riccati_next.Pvv * kkt_residual.Fv();
  kkt_residual.la() -= dtau * riccati_next.sv;
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(dimv+dimf+dimc, dimv+dimf+dimc);
  G.topLeftCorner(dimv+dimf, dimv+dimf) = kkt_matrix.Qafaf();
  G.topRightCorner(dimv+dimf, dimc) = kkt_matrix.Caf().transpose();
  G.bottomLeftCorner(dimc, dimv+dimf) = kkt_matrix.Caf();
  const Eigen::MatrixXd Ginv = G.inverse();
  Eigen::MatrixXd Qqvaf = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  Qqvaf.topLeftCorner(dimv, dimv) = kkt_matrix.Qqa();
  Qqvaf.topRightCorner(dimv, dimf) = kkt_matrix.Qqf();
  Qqvaf.bottomLeftCorner(dimv, dimv) = kkt_matrix.Qva();
  Qqvaf.bottomRightCorner(dimv, dimf) = kkt_matrix.Qvf();
  gain.computeFeedbackGain(Ginv, Qqvaf.transpose(), kkt_matrix.Cqv());
  gain.computeFeedforward(Ginv, kkt_residual.laf(), kkt_residual.C());
  Eigen::MatrixXd Pqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Pqv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Pvq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Pvv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd sq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd sv_ref = Eigen::VectorXd::Zero(dimv);
  Pqq_ref = kkt_matrix.Qqq();
  Pqq_ref += gain.Kaq().transpose() * kkt_matrix.Qqa().transpose();
  Pqv_ref = kkt_matrix.Qqv();
  Pqv_ref += gain.Kaq().transpose() * kkt_matrix.Qva().transpose();
  Pvq_ref = kkt_matrix.Qvq();
  Pvq_ref += gain.Kav().transpose() * kkt_matrix.Qqa().transpose();
  Pvv_ref = kkt_matrix.Qvv();
  Pvv_ref += gain.Kav().transpose() * kkt_matrix.Qva().transpose();
  Pqq_ref += gain.Kfq().transpose() * kkt_matrix.Qqf().transpose();
  Pqv_ref += gain.Kfq().transpose() * kkt_matrix.Qvf().transpose();
  Pvq_ref += gain.Kfv().transpose() * kkt_matrix.Qqf().transpose();
  Pvv_ref += gain.Kfv().transpose() * kkt_matrix.Qvf().transpose();
  Pqq_ref += gain.Kmuq().transpose() * kkt_matrix.Cq();
  Pqv_ref += gain.Kmuq().transpose() * kkt_matrix.Cv();
  Pvq_ref += gain.Kmuv().transpose() * kkt_matrix.Cq();
  Pvv_ref += gain.Kmuv().transpose() * kkt_matrix.Cv();
  sq_ref = riccati_next.sq - kkt_residual.lq();
  sq_ref -= riccati_next.Pqq * kkt_residual.Fq();
  sq_ref -= riccati_next.Pqv * kkt_residual.Fv();
  sq_ref -= kkt_matrix.Qqa() * gain.ka();
  sv_ref = dtau * riccati_next.sq + riccati_next.sv - kkt_residual.lv();
  sv_ref -= dtau * riccati_next.Pqq * kkt_residual.Fq();
  sv_ref -= riccati_next.Pvq * kkt_residual.Fq();
  sv_ref -= dtau * riccati_next.Pqv * kkt_residual.Fv();
  sv_ref -= riccati_next.Pvv * kkt_residual.Fv();
  sv_ref -= kkt_matrix.Qva() * gain.ka();
  sq_ref -= kkt_matrix.Qqf() * gain.kf();
  sv_ref -= kkt_matrix.Qvf() * gain.kf();
  sq_ref -= kkt_matrix.Cq().transpose() * gain.kmu();
  sv_ref -= kkt_matrix.Cv().transpose() * gain.kmu();
  EXPECT_TRUE(riccati.Pqq.isApprox(Pqq_ref));
  EXPECT_TRUE(riccati.Pqv.isApprox(Pqv_ref));
  EXPECT_TRUE(riccati.Pvq.isApprox(Pvq_ref));;
  EXPECT_TRUE(riccati.Pvv.isApprox(Pvv_ref));
  EXPECT_TRUE(riccati.sq.isApprox(sq_ref));
  EXPECT_TRUE(riccati.sv.isApprox(sv_ref));
  const Eigen::VectorXd da_ref = gain.ka() + gain.Kaq() * d.dq() + gain.Kav() * d.dv(); 
  const Eigen::VectorXd dq_next_ref = d.dq() + dtau * d.dv() + kkt_residual.Fq();
  const Eigen::VectorXd dv_next_ref = d.dv() + dtau * d.da() + kkt_residual.Fv();
  EXPECT_TRUE(d_next.dq().isApprox(dq_next_ref));
  EXPECT_TRUE(d_next.dv().isApprox(dv_next_ref));
  gain.computeFeedbackGain(Ginv, Qqvaf.transpose(), kkt_matrix.Cqv());
  gain.computeFeedforward(Ginv, kkt_residual.laf(), kkt_residual.C());

  double condensed_KKT_ref = kkt_residual.KKT_residual().squaredNorm();
  condensed_KKT_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(condensed_KKT_ref, ocp.condensedSquaredKKTErrorNorm(robot, t, dtau, s));

  ocp.computeCondensedDirection(robot, dtau, s, d);
  EXPECT_TRUE(d.df().isApprox(gain.kf()+gain.Kfq()*d.dq()+gain.Kfv()*d.dv()));
  EXPECT_TRUE(d.dmu().isApprox(gain.kmu()+gain.Kmuq()*d.dq()+gain.Kmuv()*d.dv()));

  const double step_size = 0.3;
  const Eigen::VectorXd lmd_ref 
      = s.lmd + step_size * (riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq);
  const Eigen::VectorXd gmm_ref 
      = s.gmm + step_size * (riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv);
  Eigen::VectorXd q_ref = s.q;
  robot.integrateConfiguration(d.dq(), step_size, q_ref);
  const Eigen::VectorXd v_ref = s.v + step_size * d.dv();
  const Eigen::VectorXd a_ref = s.a + step_size * d.da();
  const Eigen::VectorXd f_ref = s.f_stack() + step_size * d.df();
  const Eigen::VectorXd mu_ref = s.mu_stack() + step_size * d.dmu();
  const Eigen::VectorXd u_ref = s.u + step_size * d.du;
  const Eigen::VectorXd beta_ref = s.beta + step_size * d.dbeta;
  ocp.updatePrimal(robot, step_size, dtau, riccati, d, s);
  EXPECT_TRUE(lmd_ref.isApprox(s.lmd));
  EXPECT_TRUE(gmm_ref.isApprox(s.gmm));
  EXPECT_TRUE(q_ref.isApprox(s.q));
  EXPECT_TRUE(v_ref.isApprox(s.v));
  EXPECT_TRUE(a_ref.isApprox(s.a));
  EXPECT_TRUE(f_ref.isApprox(s.f_stack()));
  EXPECT_TRUE(mu_ref.isApprox(s.mu_stack()));
  EXPECT_TRUE(u_ref.isApprox(s.u));
  EXPECT_TRUE(beta_ref.isApprox(s.beta));

  Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  ocp.getStateFeedbackGain(Kuq, Kuv);
  Eigen::MatrixXd Kuq_ref = du_dq + du_da * gain.Kaq();
  if (dimf > 0) {
    Kuq_ref += du_df * gain.Kfq();
  }
  Eigen::MatrixXd Kuv_ref = du_dv + du_da * gain.Kav();
  if (dimf > 0) {
    Kuv_ref += du_df * gain.Kfv();
  }
  EXPECT_TRUE(Kuq.isApprox(Kuq_ref));
  EXPECT_TRUE(Kuv.isApprox(Kuv_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
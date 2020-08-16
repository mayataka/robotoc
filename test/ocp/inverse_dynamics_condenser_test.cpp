#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/inverse_dynamics_condenser.hpp"


namespace idocp {

class InverseDynamicsCondenserTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(InverseDynamicsCondenserTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  Robot robot(fixed_base_urdf_, contact_frames, baum_a, baum_b);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
  s.lmd = Eigen::VectorXd::Random(robot.dimv());
  s.gmm = Eigen::VectorXd::Random(robot.dimv());
  s.u = Eigen::VectorXd::Random(robot.dimv());
  s.beta = Eigen::VectorXd::Random(robot.dimv());
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  cost->push_back(joint_cost);
  CostFunctionData cost_data(robot);
  auto constraints = std::make_shared<Constraints>();
  auto joint_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
  auto joint_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
  constraints->push_back(joint_upper_limit);
  constraints->push_back(joint_lower_limit);
  ConstraintsData constraints_data(constraints->createConstraintsData(robot));
  constraints->setSlackAndDual(robot, constraints_data, dtau_, s.a, s.f, s.q, s.v, s.u);
  InverseDynamicsCondenser id_condenser(robot);
  id_condenser.setContactStatus(robot);
  id_condenser.linearizeCostAndConstraints(robot, cost, cost_data, constraints, 
                                           constraints_data, t_, dtau_, s);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  cost->lu(robot, cost_data, t_, dtau_, s.u, lu);
  cost->luu(robot, cost_data, t_, dtau_, s.u, luu);
  constraints->augmentDualResidual(robot, constraints_data, dtau_, lu);
  EXPECT_DOUBLE_EQ(lu.squaredNorm(), id_condenser.squaredKKTErrorNorm(dtau_));
  lu.head(robot.dim_passive()) += dtau_ * s.mu_active().tail(robot.dim_passive());
  lu -= dtau_ * s.beta;
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  id_condenser.linearizeFloatingBaseConstraint(dtau_, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.C().isZero());
  id_condenser.augmentInverseDynamicsDerivatives(dtau_, s, kkt_residual);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());  
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta));
  EXPECT_DOUBLE_EQ(lu.squaredNorm()+dtau_*dtau_*u_res.squaredNorm(), 
                   id_condenser.squaredKKTErrorNorm(dtau_));
  id_condenser.condenseInequalityConstraints(robot, constraints, 
                                             constraints_data, t_, dtau_, s);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau_, s.u, luu, lu);
  EXPECT_DOUBLE_EQ(lu.squaredNorm()+dtau_*dtau_*u_res.squaredNorm(), 
                   id_condenser.squaredKKTErrorNorm(dtau_));
  id_condenser.condenseFloatingBaseConstraint(dtau_, kkt_residual, kkt_matrix);
  EXPECT_TRUE(kkt_residual.C().isZero());
  EXPECT_TRUE(kkt_matrix.Cq().isZero());
  EXPECT_TRUE(kkt_matrix.Cv().isZero());
  EXPECT_TRUE(kkt_matrix.Ca().isZero());
  EXPECT_TRUE(kkt_matrix.Cf().isZero());
  kkt_residual.setZero();
  id_condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
  Eigen::VectorXd lu_condensed = lu + luu * u_res;
  EXPECT_TRUE(kkt_residual.lq().isApprox(du_dq.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lv().isApprox(du_dv.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.la().isApprox(du_da.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lf().isApprox((du_df.leftCols(robot.dimf()).transpose()*lu_condensed).head(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*luu*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*luu*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*luu*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*luu*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(robot.dimf()).transpose()*luu*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(robot.dimf()).transpose()*luu*du_dq));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(robot.dimf()).transpose()*luu*du_dv));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*luu*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*luu*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*luu*du_dv));
  SplitDirection d(robot);
  d.split_direction() = Eigen::VectorXd::Random(d.split_direction().size());
  id_condenser.computeCondensedDirection(dtau_, d);
  Eigen::VectorXd du_ref = u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  if (robot.dimf() > 0) {
    du_ref += du_df.leftCols(robot.dimf()) * d.df();
  }
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (lu + luu * d.du) / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsResidualL1Norm(dtau_), 
                   dtau_*u_res.lpNorm<1>());
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  if (robot.dimf() > 0) {
    s.f = Eigen::VectorXd::Random(robot.max_dimf());
  }
  s.u = Eigen::VectorXd::Random(robot.dimv());
  if (robot.dimf() > 0) {
    robot.setContactForces(s.f);
  }
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsResidualL1Norm(robot, dtau_, s), 
                   dtau_*u_res.lpNorm<1>());
}


TEST_F(InverseDynamicsCondenserTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  Robot robot(floating_base_urdf_, contact_frames, baum_a, baum_b);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
  s.lmd = Eigen::VectorXd::Random(robot.dimv());
  s.gmm = Eigen::VectorXd::Random(robot.dimv());
  s.u = Eigen::VectorXd::Random(robot.dimv());
  s.beta = Eigen::VectorXd::Random(robot.dimv());
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  cost->push_back(joint_cost);
  CostFunctionData cost_data(robot);
  auto constraints = std::make_shared<Constraints>();
  auto joint_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
  auto joint_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
  constraints->push_back(joint_upper_limit);
  constraints->push_back(joint_lower_limit);
  ConstraintsData constraints_data(constraints->createConstraintsData(robot));
  constraints->setSlackAndDual(robot, constraints_data, dtau_, s.a, s.f, s.q, s.v, s.u);
  InverseDynamicsCondenser id_condenser(robot);
  id_condenser.setContactStatus(robot);
  id_condenser.linearizeCostAndConstraints(robot, cost, cost_data, constraints, 
                                           constraints_data, t_, dtau_, s);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  cost->lu(robot, cost_data, t_, dtau_, s.u, lu);
  cost->luu(robot, cost_data, t_, dtau_, s.u, luu);
  constraints->augmentDualResidual(robot, constraints_data, dtau_, lu);
  EXPECT_DOUBLE_EQ(lu.squaredNorm(), id_condenser.squaredKKTErrorNorm(dtau_));
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  id_condenser.linearizeFloatingBaseConstraint(dtau_, s, kkt_residual);
  lu.head(robot.dim_passive()) += dtau_ * s.mu_active().tail(robot.dim_passive());
  lu -= dtau_ * s.beta;
  const int dimc = robot.dim_passive() + robot.dimf();
  Eigen::VectorXd C = Eigen::VectorXd::Zero(dimc);
  C.tail(robot.dim_passive()) = dtau_ * s.u.head(robot.dim_passive());
  EXPECT_TRUE(kkt_residual.C().isApprox(C));
  id_condenser.augmentInverseDynamicsDerivatives(dtau_, s, kkt_residual);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());  
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta));
  EXPECT_DOUBLE_EQ(lu.squaredNorm()+dtau_*dtau_*u_res.squaredNorm(), 
                   id_condenser.squaredKKTErrorNorm(dtau_));
  id_condenser.condenseInequalityConstraints(robot, constraints, 
                                             constraints_data, t_, dtau_, s);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau_, s.u, luu, lu);
  EXPECT_DOUBLE_EQ(lu.squaredNorm()+dtau_*dtau_*u_res.squaredNorm(), 
                   id_condenser.squaredKKTErrorNorm(dtau_));
  id_condenser.condenseFloatingBaseConstraint(dtau_, kkt_residual, kkt_matrix);
  C.tail(robot.dim_passive()) 
      = dtau_ * (s.u.head(robot.dim_passive())+u_res.head(robot.dim_passive()));
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(dimc, robot.dimv());
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(dimc, robot.dimv());
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(dimc, robot.dimv());
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(dimc, robot.dimf());
  Cq.bottomRows(robot.dim_passive()) = dtau_ * du_dq.topRows(robot.dim_passive());
  Cv.bottomRows(robot.dim_passive()) = dtau_ * du_dv.topRows(robot.dim_passive());
  Ca.bottomRows(robot.dim_passive()) = dtau_ * du_da.topRows(robot.dim_passive());
  Cf.bottomRows(robot.dim_passive()) = dtau_ * du_df.leftCols(robot.dimf()).topRows(robot.dim_passive());
  EXPECT_TRUE(kkt_residual.C().isApprox(C));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(Cq));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(Cv));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(Ca));
  EXPECT_TRUE(kkt_matrix.Cf().isApprox(Cf));
  kkt_residual.setZero();
  id_condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
  Eigen::VectorXd lu_condensed = lu + luu * u_res;
  EXPECT_TRUE(kkt_residual.lq().isApprox(du_dq.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lv().isApprox(du_dv.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.la().isApprox(du_da.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lf().isApprox((du_df.leftCols(robot.dimf()).transpose()*lu_condensed).head(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*luu*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*luu*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*luu*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*luu*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(robot.dimf()).transpose()*luu*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(robot.dimf()).transpose()*luu*du_dq));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(robot.dimf()).transpose()*luu*du_dv));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*luu*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*luu*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*luu*du_dv));
  SplitDirection d(robot);
  d.split_direction() = Eigen::VectorXd::Random(d.split_direction().size());
  id_condenser.computeCondensedDirection(dtau_, d);
  Eigen::VectorXd du_ref = u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  du_ref += du_df.leftCols(robot.dimf()) * d.df();
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (lu + luu * d.du) / dtau_;
  dbeta_ref.head(robot.dim_passive()) += d.dmu().tail(robot.dim_passive());
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsResidualL1Norm(dtau_), 
                   dtau_*u_res.lpNorm<1>());
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.u = Eigen::VectorXd::Random(robot.dimv());
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsResidualL1Norm(robot, dtau_, s), 
                   dtau_*u_res.lpNorm<1>());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
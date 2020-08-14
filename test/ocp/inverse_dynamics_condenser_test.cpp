#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
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
  s.set(robot);
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
  pdipm::JointSpaceConstraints constraints(robot);
  constraints.setSlackAndDual(dtau_, s.q, s.v, s.a, s.u);
  CostFunctionData data(robot);
  InverseDynamicsCondenser id_condenser(robot);
  id_condenser.setContactStatus(robot);
  id_condenser.linearizeStageCost(robot, cost, data, t_, dtau_, s);
  id_condenser.linearizeInequalityConstraints(robot, constraints, t_, dtau_, s);
  id_condenser.condenseInequalityConstraints(robot, constraints, t_, dtau_, s);
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  id_condenser.linearizeFloatingBaseConstraint(dtau_, s);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  id_condenser.condenseFloatingBaseConstraint(dtau_, s, kkt_residual, kkt_matrix);
  EXPECT_TRUE(kkt_residual.C().isZero());
  EXPECT_TRUE(kkt_matrix.Cq().isZero());
  EXPECT_TRUE(kkt_matrix.Cv().isZero());
  EXPECT_TRUE(kkt_matrix.Ca().isZero());
  EXPECT_TRUE(kkt_matrix.Cf().isZero());
  id_condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  cost->lu(robot, data, t_, dtau_, s.u, lu);
  cost->luu(robot, data, t_, dtau_, s.u, luu);
  constraints.augmentDualResidual(dtau_, lu);
  constraints.condenseSlackAndDual(dtau_, s.u, luu, lu);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());  
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  Eigen::VectorXd lu_condensed = lu + luu * u_res;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  EXPECT_TRUE(kkt_residual.lq().isApprox(du_dq.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lv().isApprox(du_dv.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.la().isApprox(du_da.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lf().isApprox((du_df.transpose()*lu_condensed).head(robot.dimf())));
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
  id_condenser.computeCondensedDirection(s, dtau_, d);
  Eigen::VectorXd du_ref = u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  du_ref += du_df.leftCols(robot.dimf()) * d.df();
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (lu + luu * d.du) / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsConstraintViolation(dtau_), 
                   dtau_*u_res.lpNorm<1>());
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.u = Eigen::VectorXd::Random(robot.dimv());
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsConstraintViolation(robot, dtau_, s), 
                   dtau_*u_res.lpNorm<1>());
  kkt_residual.setZero();
  kkt_matrix.setZero();
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  id_condenser.augmentInverseDynamicsDerivatives(robot, dtau_, s, kkt_residual);
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta));
  id_condenser.augmentFloatingBaseConstraint(dtau_, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.C().isZero());
  cost->lu(robot, data, t_, dtau_, s.u, lu);
  constraints.augmentDualResidual(dtau_, lu);
  lu -= dtau_ * s.beta;
  id_condenser.linearizeStageCost(robot, cost, data, t_, dtau_, s);
  id_condenser.linearizeInequalityConstraints(robot, constraints, t_, dtau_, s);
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  id_condenser.augmentInverseDynamicsDerivatives(dtau_, s);
  id_condenser.linearizeFloatingBaseConstraint(dtau_, s);
  EXPECT_DOUBLE_EQ(lu.squaredNorm()+u_res.squaredNorm(), 
                   id_condenser.squaredKKTErrorNorm());
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
  s.set(robot);
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
  pdipm::JointSpaceConstraints constraints(robot);
  constraints.setSlackAndDual(dtau_, s.q, s.v, s.a, s.u);
  CostFunctionData data(robot);
  InverseDynamicsCondenser id_condenser(robot);
  id_condenser.setContactStatus(robot);
  id_condenser.linearizeStageCost(robot, cost, data, t_, dtau_, s);
  id_condenser.linearizeInequalityConstraints(robot, constraints, t_, dtau_, s);
  id_condenser.condenseInequalityConstraints(robot, constraints, t_, dtau_, s);
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  id_condenser.linearizeFloatingBaseConstraint(dtau_, s);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  id_condenser.condenseFloatingBaseConstraint(dtau_, s, kkt_residual, kkt_matrix);
  id_condenser.condenseInverseDynamics(kkt_residual, kkt_matrix);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  cost->lu(robot, data, t_, dtau_, s.u, lu);
  cost->luu(robot, data, t_, dtau_, s.u, luu);
  lu.head(robot.dim_passive()) += dtau_ * s.mu_active().tail(robot.dim_passive());
  constraints.augmentDualResidual(dtau_, lu);
  constraints.condenseSlackAndDual(dtau_, s.u, luu, lu);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());  
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  Eigen::VectorXd lu_condensed = lu + luu * u_res;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  EXPECT_TRUE(kkt_residual.lq().isApprox(du_dq.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lv().isApprox(du_dv.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.la().isApprox(du_da.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lf().isApprox((du_df.transpose()*lu_condensed).head(robot.dimf())));
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
  Eigen::VectorXd C = Eigen::VectorXd::Zero(robot.dimf()+robot.dim_passive());
  C.tail(robot.dim_passive()) = dtau_ * (s.u.head(robot.dim_passive())+u_res.head(robot.dim_passive()));
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimf());
  Cq.bottomRows(robot.dim_passive()) = dtau_ * du_dq.topRows(robot.dim_passive());
  Cv.bottomRows(robot.dim_passive()) = dtau_ * du_dv.topRows(robot.dim_passive());
  Ca.bottomRows(robot.dim_passive()) = dtau_ * du_da.topRows(robot.dim_passive());
  Cf.bottomRows(robot.dim_passive()) = dtau_ * du_df.topLeftCorner(robot.dim_passive(), robot.dimf());
  EXPECT_TRUE(kkt_residual.C().isApprox(C));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(Cq));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(Cv));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(Ca));
  EXPECT_TRUE(kkt_matrix.Cf().isApprox(Cf));
  SplitDirection d(robot);
  d.split_direction() = Eigen::VectorXd::Random(d.split_direction().size());
  id_condenser.computeCondensedDirection(s, dtau_, d);
  Eigen::VectorXd du_ref = u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  du_ref += du_df.leftCols(robot.dimf()) * d.df();
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (lu + luu * d.du) / dtau_;
  dbeta_ref.head(robot.dim_passive()) += d.dmu().tail(robot.dim_passive()) / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsConstraintViolation(dtau_), 
                   dtau_*u_res.lpNorm<1>());
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.u = Eigen::VectorXd::Random(robot.dimv());
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res -= s.u;
  EXPECT_DOUBLE_EQ(id_condenser.inverseDynamicsConstraintViolation(robot, dtau_, s), 
                   dtau_*u_res.lpNorm<1>());
  kkt_residual.setZero();
  kkt_matrix.setZero();
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  id_condenser.augmentInverseDynamicsDerivatives(robot, dtau_, s, kkt_residual);
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta));
  id_condenser.augmentFloatingBaseConstraint(dtau_, s, kkt_residual);
  C.setZero();
  C.tail(robot.dim_passive()) = dtau_ * s.u.head(robot.dim_passive());
  EXPECT_TRUE(kkt_residual.C().isApprox(C));
  cost->lu(robot, data, t_, dtau_, s.u, lu);
  constraints.augmentDualResidual(dtau_, lu);
  lu.head(robot.dim_passive()) += dtau_ * s.mu_active().tail(robot.dim_passive());
  lu -= dtau_ * s.beta;
  id_condenser.linearizeStageCost(robot, cost, data, t_, dtau_, s);
  id_condenser.linearizeInequalityConstraints(robot, constraints, t_, dtau_, s);
  id_condenser.linearizeInverseDynamics(robot, dtau_, s);
  id_condenser.augmentInverseDynamicsDerivatives(dtau_, s);
  id_condenser.linearizeFloatingBaseConstraint(dtau_, s);
  EXPECT_DOUBLE_EQ(lu.squaredNorm()+u_res.squaredNorm(), 
                   id_condenser.squaredKKTErrorNorm());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
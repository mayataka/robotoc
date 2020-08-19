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
#include "idocp/ocp/inverse_dynamics.hpp"


namespace idocp {

class InverseDynamicsTest : public ::testing::Test {
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


TEST_F(InverseDynamicsTest, fixed_base) {
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
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  InverseDynamics id(robot);
  Eigen::VectorXd lu_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::MatrixXd Quu_ref = Eigen::MatrixXd::Random(robot.dimv(), 
                                                          robot.dimv());
  kkt_residual.lu = lu_ref;
  kkt_matrix.Quu = Quu_ref;
  id.linearizeInverseDynamics(robot, dtau_, s, kkt_residual);
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
  EXPECT_TRUE(kkt_residual.u_res.isApprox(u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta));
  lu_ref -= dtau_ * s.beta;
  Eigen::VectorXd lu_condensed = lu_ref + Quu_ref * u_res; 
  id.condenseInverseDynamics(kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta+du_dq.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta+du_dv.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta+du_da.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta+(du_df.leftCols(robot.dimf()).transpose()*lu_condensed).head(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
  id.condenseEqualityConstraint(dtau_, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.C().isZero());
  EXPECT_TRUE(kkt_matrix.Cq().isZero());
  EXPECT_TRUE(kkt_matrix.Cv().isZero());
  EXPECT_TRUE(kkt_matrix.Ca().isZero());
  EXPECT_TRUE(kkt_matrix.Cf().isZero());
  SplitDirection d(robot);
  d.setContactStatus(robot);
  d.dv() = Eigen::VectorXd::Random(robot.dimv());
  d.da() = Eigen::VectorXd::Random(robot.dimv());
  d.df() = Eigen::VectorXd::Random(robot.dimf());
  d.dmu() = Eigen::VectorXd::Random(robot.dim_passive()+robot.dimf());
  id.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
  Eigen::VectorXd du_ref = u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  if (robot.dimf() > 0) {
    du_ref += du_df.leftCols(robot.dimf()) * d.df();
  }
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (lu_ref + Quu_ref * du_ref) / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  std::cout << "du_dq" << std::endl;
  std::cout << du_dq << std::endl;
  std::cout << "du_dv" << std::endl;
  std::cout << du_dv << std::endl;
  std::cout << "du_da" << std::endl;
  std::cout << du_da << std::endl;
  std::cout << "du_df" << std::endl;
  std::cout << du_df << std::endl;
}



TEST_F(InverseDynamicsTest, floating_base) {
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
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  InverseDynamics id(robot);
  Eigen::VectorXd lu_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::MatrixXd Quu_ref = Eigen::MatrixXd::Random(robot.dimv(), 
                                                          robot.dimv());
  kkt_residual.lu = lu_ref;
  kkt_matrix.Quu = Quu_ref;
  id.linearizeInverseDynamics(robot, dtau_, s, kkt_residual);
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
  EXPECT_TRUE(kkt_residual.u_res.isApprox(u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta));
  lu_ref -= dtau_ * s.beta;
  Eigen::VectorXd lu_condensed = lu_ref + Quu_ref * u_res; 
  id.condenseInverseDynamics(kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*du_dq.transpose()*s.beta+du_dq.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*du_dv.transpose()*s.beta+du_dv.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*du_da.transpose()*s.beta+du_da.transpose()*lu_condensed));
  EXPECT_TRUE(kkt_residual.lf().isApprox(dtau_*du_df.leftCols(robot.dimf()).transpose()*s.beta+(du_df.leftCols(robot.dimf()).transpose()*lu_condensed).head(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
  id.condenseEqualityConstraint(dtau_, kkt_matrix, kkt_residual);
  Eigen::VectorXd C = Eigen::VectorXd::Zero(robot.dimf()+robot.dim_passive());
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimf());
  C.tail(robot.dim_passive()) = dtau_ * kkt_residual.u_res.head(robot.dim_passive());
  Cq.bottomRows(robot.dim_passive()) = dtau_ * du_dq.topRows(robot.dim_passive());
  Cv.bottomRows(robot.dim_passive()) = dtau_ * du_dv.topRows(robot.dim_passive());
  Ca.bottomRows(robot.dim_passive()) = dtau_ * du_da.topRows(robot.dim_passive());
  Cf.bottomRows(robot.dim_passive()) = dtau_ * du_df.leftCols(robot.dimf()).topRows(robot.dim_passive());
  EXPECT_TRUE(kkt_residual.C().isApprox(C));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(Cq));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(Cv));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(Ca));
  EXPECT_TRUE(kkt_matrix.Cf().isApprox(Cf));
  SplitDirection d(robot);
  d.setContactStatus(robot);
  d.dv() = Eigen::VectorXd::Random(robot.dimv());
  d.da() = Eigen::VectorXd::Random(robot.dimv());
  d.df() = Eigen::VectorXd::Random(robot.dimf());
  d.dmu() = Eigen::VectorXd::Random(robot.dim_passive()+robot.dimf());
  id.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
  Eigen::VectorXd du_ref = u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  if (robot.dimf() > 0) {
    du_ref += du_df.leftCols(robot.dimf()) * d.df();
  }
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (lu_ref + Quu_ref * d.du) / dtau_;
  dbeta_ref.head(robot.dim_passive()) += d.dmu().tail(robot.dim_passive());
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
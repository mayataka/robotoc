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
#include "idocp/ocp/split_parnmpc.hpp"


namespace idocp {

class SplitParNMPCTest : public ::testing::Test {
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


TEST_F(SplitParNMPCTest, coarseUpdate_fixed_base) {
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
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Zero(robot.dimv()).array().abs();
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
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
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  std::shared_ptr<Constraints> constraints = std::make_shared<Constraints>();
  SplitParNMPC split_parnmpc(robot, cost, constraints);
  split_parnmpc.initConstraints(robot, 2, dtau_, s);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_next);
  const Eigen::MatrixXd aux_mat_old_seed = Eigen::MatrixXd::Random(2*robot.dimv(), 
                                                                   2*robot.dimv());
  const Eigen::MatrixXd aux_mat_old = aux_mat_old_seed.transpose() * aux_mat_old_seed;
  SplitDirection d(robot);
  SplitSolution s_new_coarse(robot);
  d.set(robot);
  s_new_coarse.set(robot);
  split_parnmpc.coarseUpdate(robot, t_, dtau_, q_prev, v_prev, s, lmd_next, gmm_next, q_next, aux_mat_old, d, s_new_coarse);
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
  CostFunctionData data(robot);
  kkt_matrix.setZero();
  kkt_matrix.setContactStatus(robot);
  kkt_residual.setZero();
  kkt_residual.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  parnmpclinearizer::linearizeStageCost(robot, cost, data, t_, dtau_, s, 
                                        kkt_residual, lu)

  // KKTResidual kkt_residual(robot);
  // kkt_residual.setContactStatus(robot);
  // Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  // robot.updateKinematics(s.q, s.v, s.a);
  // parnmpclinearizer::linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual, lu);
  // EXPECT_TRUE(kkt_residual.lq().isApprox(-dtau_*q_weight.asDiagonal()*q_ref));
  // EXPECT_TRUE(kkt_residual.lv().isApprox(-dtau_*v_weight.asDiagonal()*v_ref));
  // EXPECT_TRUE(kkt_residual.la().isApprox(-dtau_*a_weight.asDiagonal()*a_ref));
  // EXPECT_TRUE(kkt_residual.lf().isApprox((-dtau_*f_weight.asDiagonal()*f_ref).head(robot.dimf())));
  // EXPECT_TRUE(lu.isApprox(-dtau_*u_weight.asDiagonal()*u_ref));
  // const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  // Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  // Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());  
  // parnmpclinearizer::linearizeDynamics(robot, dtau_, q_prev, v_prev, s, kkt_residual, u_res, du_dq, du_dv, du_da, du_df);
  // EXPECT_TRUE(kkt_residual.Fq().isApprox(q_prev));
  // EXPECT_TRUE(kkt_residual.Fv().isApprox(v_prev));
  // Eigen::VectorXd u_res_ref = Eigen::VectorXd::Zero(robot.dimv());
  // robot.setContactForces(s.f);
  // robot.RNEA(s.q, s.v, s.a, u_res_ref);
  // u_res_ref -= s.u;
  // EXPECT_TRUE(u_res_ref.isApprox(u_res));
  // Eigen::MatrixXd du_dq_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_dv_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_da_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_df_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());
  // robot.RNEADerivatives(s.q, s.v, s.a, du_dq_ref, du_dv_ref, du_da_ref);
  // robot.dRNEAPartialdFext(du_df_ref);
  // EXPECT_TRUE(du_dq_ref.isApprox(du_dq));
  // EXPECT_TRUE(du_dv_ref.isApprox(du_dv));
  // EXPECT_TRUE(du_da_ref.isApprox(du_da));
  // EXPECT_TRUE(du_df_ref.isApprox(du_df));
  // KKTMatrix kkt_matrix(robot);
  // kkt_matrix.setContactStatus(robot);
  // parnmpclinearizer::linearizeConstraints(robot, dtau_, s, u_res, du_dq, du_dv, du_da, du_df, kkt_residual, kkt_matrix);
  // Eigen::VectorXd C_ref = Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf());
  // robot.computeBaumgarteResidual(dtau_, C_ref);
  // EXPECT_TRUE(C_ref.head(robot.dim_passive()+robot.dimf()).isApprox(kkt_residual.C()));
  // Eigen::MatrixXd Cq_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  // Eigen::MatrixXd Cv_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  // Eigen::MatrixXd Ca_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  // Eigen::MatrixXd Cf_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.max_dimf());
  // robot.computeBaumgarteDerivatives(dtau_, Cq_ref, Cv_ref, Ca_ref);
  // EXPECT_TRUE(Cq_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cq()));
  // EXPECT_TRUE(Cv_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cv()));
  // EXPECT_TRUE(Ca_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Ca()));
  // EXPECT_TRUE(Cf_ref.topLeftCorner(robot.dim_passive()+robot.dimf(), robot.dimf()).isApprox(kkt_matrix.Cf()));
}


TEST_F(SplitParNMPCTest, floating_base) {
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
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Zero(robot.dimv()).array().abs();
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
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
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  std::shared_ptr<Constraints> constraints = std::make_shared<Constraints>();
  SplitParNMPC split_parnmpc(robot, cost, constraints);
  split_parnmpc.initConstraints(robot, 2, dtau_, s);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_next);
  const Eigen::MatrixXd aux_mat_old_seed = Eigen::MatrixXd::Random(2*robot.dimv(), 
                                                                   2*robot.dimv());
  const Eigen::MatrixXd aux_mat_old = aux_mat_old_seed.transpose() * aux_mat_old_seed;
  SplitDirection d(robot);
  SplitSolution s_new_coarse(robot);
  d.set(robot);
  s_new_coarse.set(robot);
  split_parnmpc.coarseUpdate(robot, t_, dtau_, q_prev, v_prev, s, lmd_next, gmm_next, q_next, aux_mat_old, d, s_new_coarse);
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
  CostFunctionData data(robot);
  kkt_matrix.setZero();
  kkt_matrix.setContactStatus(robot);
  kkt_residual.setZero();
  kkt_residual.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);



  // CostFunctionData data(robot);
  // KKTResidual kkt_residual(robot);
  // kkt_residual.setContactStatus(robot);
  // Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  // Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(robot.dimv());
  // robot.subtractConfiguration(s.q, q_ref, q_diff);
  // Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
  // parnmpclinearizer::linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual, lu);
  // EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*Jq_diff.transpose()*q_weight.asDiagonal()*q_diff));
  // EXPECT_TRUE(kkt_residual.lv().isApprox(-dtau_*v_weight.asDiagonal()*v_ref));
  // EXPECT_TRUE(kkt_residual.la().isApprox(-dtau_*a_weight.asDiagonal()*a_ref));
  // Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  // robot.normalizeConfiguration(q_prev);
  // const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  // Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  // Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());
  // parnmpclinearizer::linearizeDynamics(robot, dtau_, q_prev, v_prev, s, kkt_residual, u_res, du_dq, du_dv, du_da, du_df);
  // robot.subtractConfiguration(q_prev, s.q, q_diff);
  // EXPECT_TRUE(kkt_residual.Fq().isApprox(q_diff));
  // EXPECT_TRUE(kkt_residual.Fv().isApprox(v_prev));
  // Eigen::VectorXd u_res_ref = Eigen::VectorXd::Zero(robot.dimv());
  // robot.setContactForces(s.f);
  // robot.RNEA(s.q, s.v, s.a, u_res_ref);
  // u_res_ref -= s.u;
  // EXPECT_TRUE(u_res_ref.isApprox(u_res));
  // Eigen::MatrixXd du_dq_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_dv_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_da_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd du_df_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf());
  // robot.RNEADerivatives(s.q, s.v, s.a, du_dq_ref, du_dv_ref, du_da_ref);
  // robot.dRNEAPartialdFext(du_df_ref);
  // EXPECT_TRUE(du_dq_ref.isApprox(du_dq));
  // EXPECT_TRUE(du_dv_ref.isApprox(du_dv));
  // EXPECT_TRUE(du_da_ref.isApprox(du_da));
  // EXPECT_TRUE(du_df_ref.isApprox(du_df));
  // EXPECT_TRUE(lu.isApprox(-dtau_*u_weight.asDiagonal()*u_ref));
  // KKTMatrix kkt_matrix(robot);
  // kkt_matrix.setContactStatus(robot);
  // parnmpclinearizer::linearizeConstraints(robot, dtau_, s, u_res, du_dq, du_dv, du_da, du_df, kkt_residual, kkt_matrix);
  // Eigen::VectorXd C_ref = Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf());
  // robot.computeBaumgarteResidual(dtau_, C_ref);
  // C_ref.head(robot.dim_passive()+robot.dimf()).tail(robot.dim_passive()) = dtau_ * (s.u.head(robot.dim_passive())+u_res.head(robot.dim_passive()));
  // EXPECT_TRUE(C_ref.head(robot.dim_passive()+robot.dimf()).isApprox(kkt_residual.C()));
  // Eigen::MatrixXd Cq_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  // Eigen::MatrixXd Cv_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  // Eigen::MatrixXd Ca_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  // Eigen::MatrixXd Cf_ref = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.max_dimf());
  // robot.computeBaumgarteDerivatives(dtau_, Cq_ref.topRows(robot.dimf()), Cv_ref.topRows(robot.dimf()), Ca_ref.topRows(robot.dimf()));
  // Cq_ref.topRows(robot.dim_passive()+robot.dimf()).bottomRows(robot.dim_passive()) = dtau_ * du_dq.topRows(robot.dim_passive());
  // Cv_ref.topRows(robot.dim_passive()+robot.dimf()).bottomRows(robot.dim_passive()) = dtau_ * du_dv.topRows(robot.dim_passive());
  // Ca_ref.topRows(robot.dim_passive()+robot.dimf()).bottomRows(robot.dim_passive()) = dtau_ * du_da.topRows(robot.dim_passive());
  // Cf_ref.topLeftCorner(robot.dim_passive()+robot.dimf(), robot.dimf()).bottomRows(robot.dim_passive()) = dtau_ * du_df.topLeftCorner(robot.dim_passive()+robot.dimf(), robot.dimf()).topRows(robot.dim_passive());
  // EXPECT_TRUE(Cq_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cq()));
  // EXPECT_TRUE(Cv_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cv()));
  // EXPECT_TRUE(Ca_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Ca()));
  // EXPECT_TRUE(Cf_ref.topLeftCorner(robot.dim_passive()+robot.dimf(), robot.dimf()).isApprox(kkt_matrix.Cf()));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
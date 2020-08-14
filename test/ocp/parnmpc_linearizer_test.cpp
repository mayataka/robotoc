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
#include "idocp/ocp/parnmpc_linearizer.hpp"


namespace idocp {

class ParNMPCLinearizerTest : public ::testing::Test {
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


TEST_F(ParNMPCLinearizerTest, fixed_base) {
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf());
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
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  CostFunctionData data(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  ParNMPCLinearizer linearizer(robot);
  linearizer.linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*q_weight.asDiagonal()*(s.q-q_ref)));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*v_weight.asDiagonal()*(s.v-v_ref)));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*a_weight.asDiagonal()*(s.a-a_ref)));
  EXPECT_TRUE(kkt_residual.lf().isApprox((dtau_*f_weight.asDiagonal()*(s.f-f_ref)).head(robot.dimf())));
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_next);
  kkt_residual.setZero();
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, lmd_next, 
                                    gmm_next, q_next, kkt_residual, kkt_matrix);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((lmd_next-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+gmm_next)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq()
              .isApprox(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  linearizer.linearizeTerminalCost(robot, cost, data, t_, s);
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, 
                                    kkt_residual, kkt_matrix);
  Eigen::VectorXd phiq = qf_weight.asDiagonal()*(s.q-q_ref);
  Eigen::VectorXd phiv = vf_weight.asDiagonal()*(s.v-v_ref);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((phiq-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+phiv)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq()
              .isApprox(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  robot.updateKinematics(s.q, s.v, s.a);
  linearizer.linearizeContactConstraints(robot, dtau_, kkt_residual, kkt_matrix);
  Eigen::VectorXd C_ref = Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf());
  robot.computeBaumgarteResidual(dtau_, C_ref);
  EXPECT_TRUE(C_ref.head(robot.dim_passive()+robot.dimf()).isApprox(kkt_residual.C()));
  Eigen::MatrixXd Cq_ref 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  Eigen::MatrixXd Cv_ref 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  Eigen::MatrixXd Ca_ref 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  robot.computeBaumgarteDerivatives(dtau_, Cq_ref, Cv_ref, Ca_ref);
  EXPECT_TRUE(Cq_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cq()));
  EXPECT_TRUE(Cv_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cv()));
  EXPECT_TRUE(Ca_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Ca()));
}


TEST_F(ParNMPCLinearizerTest, fixed_base_KKT_error) {
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf());
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
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  CostFunctionData data(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  ParNMPCLinearizer linearizer(robot);
  linearizer.linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, lmd_next, 
                                    gmm_next, q_next, kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau_, kkt_residual, kkt_matrix);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  kkt_residual_ref.lq() = dtau_*q_weight.asDiagonal()*(s.q-q_ref);
  kkt_residual_ref.lv() = dtau_*v_weight.asDiagonal()*(s.v-v_ref);
  kkt_residual_ref.la() = dtau_*a_weight.asDiagonal()*(s.a-a_ref);
  kkt_residual_ref.lf() = (dtau_*f_weight.asDiagonal()*(s.f-f_ref)).head(robot.dimf());
  kkt_residual_ref.Fq() = q_prev - s.q + dtau_ * s.v;
  kkt_residual_ref.Fv() = v_prev - s.v + dtau_ * s.a;
  kkt_residual_ref.lq() += lmd_next - s.lmd;
  kkt_residual_ref.lv() += dtau_ * s.lmd - s.gmm + gmm_next;
  kkt_residual_ref.la() += dtau_ * s.gmm;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.Fq().isApprox(kkt_residual_ref.Fq()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  kkt_residual_ref.setZero();
  linearizer.linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual);
  linearizer.linearizeTerminalCost(robot, cost, data, t_, s);
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, 
                                    kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau_, kkt_residual, kkt_matrix);
  Eigen::VectorXd phiq = qf_weight.asDiagonal()*(s.q-q_ref);
  Eigen::VectorXd phiv = vf_weight.asDiagonal()*(s.v-v_ref);
  kkt_residual_ref.lq() = dtau_*q_weight.asDiagonal()*(s.q-q_ref);
  kkt_residual_ref.lv() = dtau_*v_weight.asDiagonal()*(s.v-v_ref);
  kkt_residual_ref.la() = dtau_*a_weight.asDiagonal()*(s.a-a_ref);
  kkt_residual_ref.lf() = (dtau_*f_weight.asDiagonal()*(s.f-f_ref)).head(robot.dimf());
  kkt_residual_ref.Fq() = q_prev - s.q + dtau_ * s.v;
  kkt_residual_ref.Fv() = v_prev - s.v + dtau_ * s.a;
  kkt_residual_ref.lq() += phiq - s.lmd;
  kkt_residual_ref.lv() += dtau_ * s.lmd - s.gmm + phiv;
  kkt_residual_ref.la() += dtau_ * s.gmm;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.Fq().isApprox(kkt_residual_ref.Fq()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  cost->lf(robot, data, t_, dtau_, s.f, kkt_residual_ref.lf());
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
}


TEST_F(ParNMPCLinearizerTest, floating_base) {
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf());
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
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  CostFunctionData data(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, q_ref, q_diff);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
  ParNMPCLinearizer linearizer(robot);
  linearizer.linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.lq().isApprox(dtau_*Jq_diff.transpose()*q_weight.asDiagonal()*q_diff));
  EXPECT_TRUE(kkt_residual.lv().isApprox(dtau_*v_weight.asDiagonal()*(s.v-v_ref)));
  EXPECT_TRUE(kkt_residual.la().isApprox(dtau_*a_weight.asDiagonal()*(s.a-a_ref)));
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_next);
  kkt_residual.setZero();
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, lmd_next, 
                                    gmm_next, q_next, kkt_residual, kkt_matrix);
  Eigen::MatrixXd Jsub_minus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Jsub_plus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, q_diff);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, Jsub_minus);
  robot.dSubtractdConfigurationPlus(s.q, q_next, Jsub_plus);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_diff+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq()
              .isApprox((Jsub_minus.transpose()*s.lmd+Jsub_plus.transpose()*lmd_next)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+gmm_next)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Jsub_minus));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  linearizer.linearizeTerminalCost(robot, cost, data, t_, s);
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, 
                                    kkt_residual, kkt_matrix);
  robot.subtractConfiguration(q_prev, s.q, q_diff);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_diff+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  robot.subtractConfiguration(s.q, q_ref, q_diff);
  robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
  Eigen::VectorXd phiq = Jq_diff.transpose()*qf_weight.asDiagonal()*q_diff;
  Eigen::VectorXd phiv = vf_weight.asDiagonal()*(s.v-v_ref);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, Jsub_minus);
  EXPECT_TRUE(kkt_residual.lq().isApprox((phiq+Jsub_minus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+phiv)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Jsub_minus));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  robot.updateKinematics(s.q, s.v, s.a);
  linearizer.linearizeContactConstraints(robot, dtau_, kkt_residual, kkt_matrix);
  Eigen::VectorXd C_ref = Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf());
  robot.computeBaumgarteResidual(dtau_, C_ref);
  EXPECT_TRUE(C_ref.head(robot.dim_passive()+robot.dimf()).isApprox(kkt_residual.C()));
  Eigen::MatrixXd Cq_ref 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  Eigen::MatrixXd Cv_ref 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  Eigen::MatrixXd Ca_ref 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), robot.dimv());
  robot.computeBaumgarteDerivatives(dtau_, Cq_ref, Cv_ref, Ca_ref);
  EXPECT_TRUE(Cq_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cq()));
  EXPECT_TRUE(Cv_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Cv()));
  EXPECT_TRUE(Ca_ref.topRows(robot.dim_passive()+robot.dimf()).isApprox(kkt_matrix.Ca()));
}


TEST_F(ParNMPCLinearizerTest, floating_base_KKT_error) {
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf());
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
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  CostFunctionData data(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  ParNMPCLinearizer linearizer(robot);
  linearizer.linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_next);
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, lmd_next, 
                                    gmm_next, q_next, kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau_, kkt_residual, kkt_matrix);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, q_ref, q_diff);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
  kkt_residual_ref.lq() = dtau_*Jq_diff.transpose()*q_weight.asDiagonal()*q_diff;
  kkt_residual_ref.lv() = dtau_*v_weight.asDiagonal()*(s.v-v_ref);
  kkt_residual_ref.la() = dtau_*a_weight.asDiagonal()*(s.a-a_ref);
  Eigen::MatrixXd Jsub_minus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Jsub_plus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, q_diff);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, Jsub_minus);
  robot.dSubtractdConfigurationPlus(s.q, q_next, Jsub_plus);
  kkt_residual_ref.Fq() = q_diff + dtau_ * s.v;
  kkt_residual_ref.Fv() = v_prev - s.v + dtau_ * s.a;
  kkt_residual_ref.lq() += Jsub_minus.transpose() * s.lmd + Jsub_plus.transpose() * lmd_next;
  kkt_residual_ref.lv() += dtau_ * s.lmd - s.gmm + gmm_next;
  kkt_residual_ref.la() += dtau_ * s.gmm;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.Fq().isApprox(kkt_residual_ref.Fq()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  kkt_residual_ref.setZero();
  linearizer.linearizeStageCost(robot, cost, data, t_, dtau_, s, kkt_residual);
  linearizer.linearizeTerminalCost(robot, cost, data, t_, s);
  linearizer.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, 
                                    kkt_residual, kkt_matrix);
  linearizer.linearizeContactConstraints(robot, dtau_, kkt_residual, kkt_matrix);
  robot.subtractConfiguration(s.q, q_ref, q_diff);
  robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
  Eigen::VectorXd phiq = Jq_diff.transpose()*qf_weight.asDiagonal()*q_diff;
  Eigen::VectorXd phiv = vf_weight.asDiagonal()*(s.v-v_ref);
  kkt_residual_ref.lq() = dtau_*Jq_diff.transpose()*q_weight.asDiagonal()*q_diff;
  kkt_residual_ref.lv() = dtau_*v_weight.asDiagonal()*(s.v-v_ref);
  kkt_residual_ref.la() = dtau_*a_weight.asDiagonal()*(s.a-a_ref);
  robot.subtractConfiguration(q_prev, s.q, q_diff);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, Jsub_minus);
  kkt_residual_ref.Fq() = q_diff + dtau_ * s.v;
  kkt_residual_ref.Fv() = v_prev - s.v + dtau_ * s.a;
  kkt_residual_ref.lq() += Jsub_minus.transpose() * s.lmd + phiq;
  kkt_residual_ref.lv() += dtau_ * s.lmd - s.gmm + phiv;
  kkt_residual_ref.la() += dtau_ * s.gmm;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.Fq().isApprox(kkt_residual_ref.Fq()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  cost->lf(robot, data, t_, dtau_, s.f, kkt_residual_ref.lf());
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
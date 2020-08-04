#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"
#include "idocp/constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "idocp/quadruped/cost_function.hpp"
#include "idocp/quadruped/constraints.hpp"


namespace idocp {

class FloatingBaseSplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    contact_frames_ = {14, 24, 34, 44};
    baum_on_velocity_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    baum_on_position_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    robot_ = Robot(urdf_, contact_frames_, baum_on_velocity_, baum_on_position_);
    cost_ = quadruped::CostFunction(robot_);
    constraints_ = quadruped::Constraints(robot_);
    ocp_ = SplitOCP(robot_, cost_, constraints_);
    factorizer_ = RiccatiMatrixFactorizer(robot_);
    inverter_ = RiccatiMatrixInverter(robot_);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    time_step_ = rnd()%10;
    dimq_ = robot_.dimq();
    dimv_ = robot_.dimv();
    max_dimf_ = robot_.max_dimf();
    dim_passive_ = robot_.dim_passive();
    max_dimc_ = max_dimf_ + dim_passive_;
    q_ = Eigen::VectorXd::Random(dimq_);
    v_ = Eigen::VectorXd::Random(dimv_);
    a_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
    f_ = Eigen::VectorXd::Random(max_dimf_);
    lmd_ = Eigen::VectorXd::Random(dimv_);
    gmm_ = Eigen::VectorXd::Random(dimv_);
    mu_ = Eigen::VectorXd::Random(max_dimc_);
    q_next_ = Eigen::VectorXd::Random(dimq_);
    v_next_ = Eigen::VectorXd::Random(dimv_);
    lmd_next_ = Eigen::VectorXd::Random(dimv_);
    gmm_next_ = Eigen::VectorXd::Random(dimv_);
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  Robot robot_;
  SplitOCP ocp_;
  RiccatiMatrixFactorizer factorizer_;
  RiccatiMatrixInverter inverter_;
  quadruped::CostFunction cost_;
  quadruped::Constraints constraints_;
  double t_, dtau_, baum_on_velocity_, baum_on_position_;
  int time_step_, dimq_, dimv_, max_dimf_, dim_passive_, max_dimc_;
  std::vector<int> contact_frames_;
  Eigen::VectorXd q_, v_, a_, u_, f_, lmd_, gmm_, mu_, q_next_, v_next_, 
                  lmd_next_, gmm_next_;
};


TEST_F(FloatingBaseSplitOCPTest, isFeasible) {
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  EXPECT_EQ(ocp_.isFeasible(q_, v_, a_, u_), 
            joint_space_constraints_ref.isFeasible(q_, v_, a_, u_));
}


TEST_F(FloatingBaseSplitOCPTest, withoutActiveContacts) {
  ASSERT_TRUE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.dim_passive() == 6);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  robot_.generateFeasibleConfiguration(q_next_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  std::vector<bool> active_contacts = {false, false, false, false};
  robot_.setActiveContacts(active_contacts);
  ASSERT_TRUE(robot_.max_dimf() == 12);
  ASSERT_TRUE(robot_.dimf() == 0);
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref.setTimeStep(time_step_);
  joint_space_constraints_ref.setSlackAndDual(dtau_, q_, v_, a_, u_);
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
  ocp_.linearizeOCP(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                    lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::VectorXd q_res, v_res;
  q_res = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_, q_next_, q_res);
  q_res.noalias() += dtau_ * v_;
  v_res = v_ + dtau_ * a_ - v_next_;
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv_);
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv_, dimv_);
  cost_.lu(t_, dtau_, u_, lu);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lu);
  cost_.luu(t_, dtau_, u_, luu);
  joint_space_constraints_ref.condenseSlackAndDual(dtau_, u_, luu, lu);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv_);
  robot_.RNEA(q_, v_, a_, u_res);
  u_res -= u_;
  const Eigen::VectorXd lu_condensed = lu + luu * u_res;
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  cost_.setConfigurationJacobian(robot_, q_);
  cost_.lq(t_, dtau_, q_, v_, a_, lq);
  cost_.lv(t_, dtau_, q_, v_, a_, lv);
  cost_.la(t_, dtau_, q_, v_, a_, la);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lq, lv, la);
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  robot_.RNEADerivatives(q_, v_, a_, du_dq, du_dv, du_da);
  lq += du_dq.transpose() * lu_condensed;
  lv += du_dv.transpose() * lu_condensed;
  la += du_da.transpose() * lu_condensed;
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::VectorXd C_res = Eigen::VectorXd::Zero(max_dimc_);
  Cq.topRows(dim_passive_) = dtau_ * du_dq.topRows(dim_passive_);
  Cv.topRows(dim_passive_) = dtau_ * du_dv.topRows(dim_passive_);
  Ca.topRows(dim_passive_) = dtau_ * du_da.topRows(dim_passive_);
  C_res.head(dim_passive_) = dtau_ * u_.head(dim_passive_);
  C_res.head(dim_passive_) += dtau_ * u_res.head(dim_passive_);
  lq += Cq.topRows(dim_passive_).transpose() * mu_.head(dim_passive_);
  lv += Cv.topRows(dim_passive_).transpose() * mu_.head(dim_passive_);
  la += Ca.topRows(dim_passive_).transpose() * mu_.head(dim_passive_);
  Eigen::MatrixXd Qqq = du_dq.transpose() * luu * du_dq;
  Eigen::MatrixXd Qqv = du_dq.transpose() * luu * du_dv;
  Eigen::MatrixXd Qqa = du_dq.transpose() * luu * du_da;
  Eigen::MatrixXd Qvv = du_dv.transpose() * luu * du_dv;
  Eigen::MatrixXd Qva = du_dv.transpose() * luu * du_da;
  Eigen::MatrixXd Qaa = du_da.transpose() * luu * du_da;
  joint_space_constraints_ref.condenseSlackAndDual(dtau_, q_, v_, a_, Qqq, Qvv, 
                                                   Qaa, lq, lv, la);
  cost_.augment_lqq(t_, dtau_, q_, v_, a_, Qqq);
  cost_.augment_lvv(t_, dtau_, q_, v_, a_, Qvv);
  cost_.augment_laa(t_, dtau_, q_, v_, a_, Qaa);
  Eigen::MatrixXd Qvq = Qqv.transpose();
  Eigen::MatrixXd gen_mat = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pqq_next = gen_mat * gen_mat.transpose();
  const Eigen::MatrixXd Pqv_next = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pvq_next = Eigen::MatrixXd::Random(dimv_, dimv_);
  gen_mat = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pvv_next = gen_mat * gen_mat.transpose();
  const Eigen::VectorXd sq_next = Eigen::VectorXd::Random(dimv_);
  const Eigen::VectorXd sv_next = Eigen::VectorXd::Random(dimv_);
  Eigen::MatrixXd Pqq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pqv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::VectorXd sq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd sv = Eigen::VectorXd::Zero(dimv_);
  ocp_.backwardRiccatiRecursion(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next,  
                                sq_next, sv_next, Pqq, Pqv, Pvq, Pvv, sq, sv);
  Eigen::MatrixXd Pqq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pqv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  factorizer_.setIntegrationSensitivities(robot_, dtau_, q_, v_);
  factorizer_.factorize(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next, 
                        Qqq, Qqv, Qvq, Qvv);
  factorizer_.factorize(dtau_, Pqv_next, Pvv_next, Qqa, Qva);
  factorizer_.factorize(dtau_, Pvv_next, Qaa);
  la += dtau_ * Pvq_next * q_res;
  la += dtau_ * Pvv_next * v_res;
  la -= dtau_ * sv_next;
  Eigen::MatrixXd Kaq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Kav = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Kmuq = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Kmuv = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd kmu = Eigen::VectorXd::Zero(max_dimc_);
  inverter_.invert(Qqa, Qva, Qaa, Cq, Cv, Ca, la, C_res, Kaq, Kav, Kmuq, Kmuv, 
                   ka, kmu);
  Pqq_ref = Qqq;
  Pqq_ref += Kaq.transpose() * Qqa.transpose();
  Pqv_ref = Qqv;
  Pqv_ref += Kaq.transpose() * Qva.transpose();
  Pvq_ref = Qvq;
  Pvq_ref += Kav.transpose() * Qqa.transpose();
  Pvv_ref = Qvv;
  Pvv_ref += Kav.transpose() * Qva.transpose();
  Pqq_ref += Kmuq.topRows(dim_passive_).transpose() * Cq.topRows(dim_passive_);
  Pqv_ref += Kmuq.topRows(dim_passive_).transpose() * Cv.topRows(dim_passive_);
  Pvq_ref += Kmuv.topRows(dim_passive_).transpose() * Cq.topRows(dim_passive_);
  Pvv_ref += Kmuv.topRows(dim_passive_).transpose() * Cv.topRows(dim_passive_);
  Eigen::VectorXd sq_ref = sq_next - lq;
  sq_ref -= Pqq_next * q_res;
  sq_ref -= Pqv_next * v_res;
  sq_ref -= Qqa * ka;
  Eigen::VectorXd sv_ref = dtau_ * sq_next + sv_next - lv;
  sv_ref -= dtau_ * Pqq_next * q_res;
  sv_ref -= Pvq_next * q_res;
  sv_ref -= dtau_ * Pqv_next * v_res;
  sv_ref -= Pvv_next * v_res;
  sv_ref -= Qva * ka;
  sq_ref -= Cq.topRows(dim_passive_).transpose() * kmu.head(dim_passive_);
  sv_ref -= Cv.topRows(dim_passive_).transpose() * kmu.head(dim_passive_);
  EXPECT_TRUE(Pqq.isApprox(Pqq_ref));
  EXPECT_TRUE(Pqv.isApprox(Pqv_ref));
  EXPECT_TRUE(Pvq.isApprox(Pvq_ref));
  EXPECT_TRUE(Pvv.isApprox(Pvv_ref));
  EXPECT_TRUE(sq.isApprox(sq_ref));
  EXPECT_TRUE(sv.isApprox(sv_ref));
  std::cout << "Pqq - Pqq_ref" << std::endl;
  std::cout << Pqq - Pqq_ref << "\n" << std::endl;
  std::cout << "Pqv - Pqv_ref" << std::endl;
  std::cout << Pqv - Pqv_ref  << "\n" << std::endl;
  std::cout << "Pvq - Pvq_ref" << std::endl;
  std::cout << Pvq - Pvq_ref  << "\n" << std::endl;
  std::cout << "Pvv - Pvv_ref" << std::endl;
  std::cout << Pvv - Pvv_ref << "\n" << std::endl;
  std::cout << "sq - sq_ref" << std::endl;
  std::cout << sq - sq_ref  << "\n" << std::endl;
  std::cout << "sv - sv_ref" << std::endl;
  std::cout << sv - sv_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv_);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv_);
  Eigen::VectorXd dq_next = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd dv_next = Eigen::VectorXd::Zero(dimv_);
  ocp_.forwardRiccatiRecursion(dtau_, dq, dv, dq_next, dv_next);
  const Eigen::VectorXd da = ka + Kaq * dq + Kav * dv; 
  const Eigen::VectorXd dq_next_ref = dq + dtau_ * dv + q_res;
  const Eigen::VectorXd dv_next_ref = dv + dtau_ * da + v_res;
  EXPECT_TRUE(dq_next.isApprox(dq_next_ref));
  EXPECT_TRUE(dv_next.isApprox(dv_next_ref));
  ocp_.computeCondensedDirection(dtau_, dq, dv);
  Eigen::VectorXd dmu 
      = kmu.head(dim_passive_) + Kmuq.topRows(dim_passive_) * dq 
                               + Kmuv.topRows(dim_passive_) * dv;
  Eigen::VectorXd du = u_res + du_dq * dq + du_dv * dv + du_da * da;
  joint_space_constraints_ref.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  const double max_primal_step_size = ocp_.maxPrimalStepSize();
  const double max_dual_step_size = ocp_.maxDualStepSize();
  const double max_primal_step_size_ref 
      = joint_space_constraints_ref.maxSlackStepSize();
  const double max_dual_step_size_ref 
      = joint_space_constraints_ref.maxDualStepSize();
  EXPECT_DOUBLE_EQ(max_primal_step_size, max_primal_step_size_ref);
  EXPECT_DOUBLE_EQ(max_dual_step_size, max_dual_step_size_ref);
  const std::pair<double, double> cost_and_violation 
      = ocp_.costAndConstraintsViolation(robot_, max_primal_step_size, t_, 
                                         dtau_, q_, v_, a_, u_, q_next_, 
                                         v_next_, dq, dv, dq_next, dv_next);
  Eigen::VectorXd q_tmp = q_;
  robot_.integrateConfiguration(dq, max_primal_step_size, q_tmp);
  const Eigen::VectorXd v_tmp = v_ + max_primal_step_size * dv;
  const Eigen::VectorXd a_tmp = a_ + max_primal_step_size * da;
  const Eigen::VectorXd u_tmp = u_ + max_primal_step_size * du;
  cost_.setConfigurationJacobian(robot_, q_tmp);
  const double cost_ref 
      = cost_.l(t_, dtau_, q_tmp, v_tmp, a_tmp, u_tmp, f_) 
        + joint_space_constraints_ref.costSlackBarrier(max_primal_step_size);
  Eigen::VectorXd q_res_tmp = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_tmp, q_next_, q_res_tmp);
  q_res_tmp += dtau_ * v_tmp - max_primal_step_size * dq_next;
  const Eigen::VectorXd v_res_tmp 
      = v_tmp - v_next_ + dtau_ * a_tmp - max_primal_step_size * dv_next;
  Eigen::VectorXd u_res_tmp = Eigen::VectorXd::Zero(dimv_);
  robot_.RNEA(q_tmp, v_tmp, a_tmp, u_res_tmp);
  u_res_tmp -= u_tmp;
  const double violation_ref
      = q_res_tmp.lpNorm<1>() + v_res_tmp.lpNorm<1>() 
          + dtau_ * u_res_tmp.lpNorm<1>() 
          + joint_space_constraints_ref.residualL1Nrom(dtau_, q_tmp, v_tmp, a_tmp, u_tmp)
          + dtau_ * u_tmp.head(dim_passive_).lpNorm<1>();
  EXPECT_DOUBLE_EQ(cost_and_violation.first, cost_ref);
  EXPECT_DOUBLE_EQ(cost_and_violation.second, violation_ref);
  ocp_.updateDual(max_dual_step_size);
  joint_space_constraints_ref.updateDual(max_dual_step_size);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  Eigen::VectorXd q_ref = q_;
  Eigen::VectorXd v_ref = v_;
  Eigen::VectorXd a_ref = a_;
  Eigen::VectorXd u_ref = u_;
  Eigen::VectorXd beta_ref = beta;
  Eigen::VectorXd lmd_ref = lmd_;
  Eigen::VectorXd gmm_ref = gmm_;
  ocp_.updatePrimal(robot_, max_primal_step_size, dtau_, dq, dv, Pqq, Pqv, Pvq, 
                    Pvv, sq, sv, q_, v_, a_, u_, beta, lmd_, gmm_);
  robot_.integrateConfiguration(dq, max_primal_step_size, q_ref);
  v_ref += max_primal_step_size * dv;
  a_ref += max_primal_step_size * da;
  u_ref += max_primal_step_size * du;
  lu -= dtau_ * beta_ref;
  beta_ref += max_primal_step_size * lu / dtau_;
  beta_ref += max_primal_step_size * luu * du / dtau_;
  lmd_ref += max_primal_step_size * (Pqq * dq + Pqv * dv - sq);
  gmm_ref += max_primal_step_size * (Pvq * dq + Pvv * dv - sv);
  joint_space_constraints_ref.updateSlack(max_primal_step_size);
  EXPECT_TRUE(q_.isApprox(q_ref));
  EXPECT_TRUE(v_.isApprox(v_ref));
  EXPECT_TRUE(a_.isApprox(a_ref));
  EXPECT_TRUE(u_.isApprox(u_ref));
  EXPECT_TRUE(beta.isApprox(beta_ref));
  EXPECT_TRUE(lmd_.isApprox(lmd_ref));
  EXPECT_TRUE(gmm_.isApprox(gmm_ref));
}


TEST_F(FloatingBaseSplitOCPTest, withoutContacts) {
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
  robot_.generateFeasibleConfiguration(q_);
  robot_.generateFeasibleConfiguration(q_next_);
  ASSERT_TRUE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.max_dimf() == 12);
  ASSERT_TRUE(robot_.dim_passive() == 6);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  ocp_.linearizeOCP(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                    lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::MatrixXd gen_mat = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pqq_next = gen_mat * gen_mat.transpose();
  const Eigen::MatrixXd Pqv_next = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pvq_next = Eigen::MatrixXd::Random(dimv_, dimv_);
  gen_mat = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pvv_next = gen_mat * gen_mat.transpose();
  const Eigen::VectorXd sq_next = Eigen::VectorXd::Random(dimv_);
  const Eigen::VectorXd sv_next = Eigen::VectorXd::Random(dimv_);
  Eigen::MatrixXd Pqq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pqv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::VectorXd sq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd sv = Eigen::VectorXd::Zero(dimv_);
  ocp_.backwardRiccatiRecursion(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next,  
                                sq_next, sv_next, Pqq, Pqv, Pvq, Pvv, sq, sv);
  const Eigen::VectorXd dq = Eigen::VectorXd::Zero(dimv_);
  const Eigen::VectorXd dv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd dq_next = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd dv_next = Eigen::VectorXd::Zero(dimv_);
  ocp_.forwardRiccatiRecursion(dtau_, dq, dv, dq_next, dv_next);
  ocp_.computeCondensedDirection(dtau_, dq, dv);
  const double primal_step_size = ocp_.maxPrimalStepSize();
  const double dual_step_size = ocp_.maxPrimalStepSize();
  const std::pair<double, double> cost_and_violation 
      = ocp_.costAndConstraintsViolation(robot_, primal_step_size, t_, 
                                         dtau_, q_, v_, a_, u_, q_next_, 
                                         v_next_, dq, dv, dq_next, dv_next);
  const Eigen::VectorXd beta = Eigen::VectorXd::Zero(dimv_);
  const double KKT_error 
      = ocp_.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_,
                                 beta, lmd_next_, gmm_next_, q_next_, v_next_);
  robot_ = Robot(urdf_);
  ASSERT_TRUE(robot_.has_floating_base());
  ASSERT_EQ(robot_.max_dimf(), 0);
  ASSERT_EQ(robot_.dimf(), 0);
  cost_ = quadruped::CostFunction(robot_);
  constraints_ = quadruped::Constraints(robot_);
  ocp_ = SplitOCP(robot_, cost_, constraints_);
  ocp_.set_f(f_.head(robot_.dimf()));
  ocp_.set_mu(mu_.head(robot_.dim_passive()));
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  ocp_.linearizeOCP(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                    lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::MatrixXd Pqq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pqv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::VectorXd sq_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd sv_ref = Eigen::VectorXd::Zero(dimv_);
  ocp_.backwardRiccatiRecursion(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next,  
                                sq_next, sv_next, Pqq_ref, Pqv_ref, Pvq_ref, 
                                Pvv_ref, sq_ref, sv_ref);
  EXPECT_TRUE(Pqq.isApprox(Pqq_ref));
  EXPECT_TRUE(Pqv.isApprox(Pqv_ref));
  EXPECT_TRUE(Pvq.isApprox(Pvq_ref));
  EXPECT_TRUE(Pvv.isApprox(Pvv_ref));
  EXPECT_TRUE(sq.isApprox(sq_ref));
  EXPECT_TRUE(sv.isApprox(sv_ref));
  Eigen::VectorXd dq_next_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd dv_next_ref = Eigen::VectorXd::Zero(dimv_);
  ocp_.forwardRiccatiRecursion(dtau_, dq, dv, dq_next_ref, dv_next_ref);
  EXPECT_TRUE(dq_next.isApprox(dq_next_ref));
  EXPECT_TRUE(dv_next.isApprox(dv_next_ref));
  ocp_.computeCondensedDirection(dtau_, dq, dv);
  const double primal_step_size_ref = ocp_.maxPrimalStepSize();
  const double dual_step_size_ref = ocp_.maxPrimalStepSize();
  EXPECT_DOUBLE_EQ(primal_step_size, primal_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const std::pair<double, double> cost_and_violation_ref 
      = ocp_.costAndConstraintsViolation(robot_, primal_step_size, t_, 
                                         dtau_, q_, v_, a_, u_, q_next_, 
                                         v_next_, dq, dv, dq_next, dv_next);
  EXPECT_DOUBLE_EQ(cost_and_violation.first, cost_and_violation_ref.first);
  EXPECT_DOUBLE_EQ(cost_and_violation.second, cost_and_violation_ref.second);
  const double KKT_error_ref
      = ocp_.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_,
                                 beta, lmd_next_, gmm_next_, q_next_, v_next_);
  EXPECT_DOUBLE_EQ(KKT_error, KKT_error_ref);
}


TEST_F(FloatingBaseSplitOCPTest, withActiveContacts) {
  ASSERT_TRUE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.dim_passive() == 6);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  robot_.generateFeasibleConfiguration(q_next_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<contact_frames_.size(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setActiveContacts(active_contacts);
  ASSERT_TRUE(robot_.max_dimf() == 12);
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  const int dimf = robot_.dimf();
  const int dimc = robot_.dim_passive() + robot_.dimf();
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref.setTimeStep(time_step_);
  joint_space_constraints_ref.setSlackAndDual(dtau_, q_, v_, a_, u_);
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
  ocp_.linearizeOCP(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                    lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::VectorXd q_res, v_res;
  q_res = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_, q_next_, q_res);
  q_res.noalias() += dtau_ * v_;
  v_res = v_ + dtau_ * a_ - v_next_;
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv_);
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv_, dimv_);
  cost_.lu(t_, dtau_, u_, lu);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lu);
  cost_.luu(t_, dtau_, u_, luu);
  joint_space_constraints_ref.condenseSlackAndDual(dtau_, u_, luu, lu);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv_);
  robot_.setContactForces(f_);
  robot_.RNEA(q_, v_, a_, u_res);
  u_res -= u_;
  const Eigen::VectorXd lu_condensed = lu + luu * u_res;
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(max_dimf_);
  cost_.setConfigurationJacobian(robot_, q_);
  cost_.lq(t_, dtau_, q_, v_, a_, lq);
  cost_.lv(t_, dtau_, q_, v_, a_, lv);
  cost_.la(t_, dtau_, q_, v_, a_, la);
  cost_.setContactStatus(robot_);
  cost_.lf(t_, dtau_, f_, lf);
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  joint_space_constraints_ref.augmentDualResidual(dtau_, lq, lv, la);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  robot_.RNEADerivatives(q_, v_, a_, du_dq, du_dv, du_da);
  robot_.updateKinematics(q_, v_, a_);
  robot_.dRNEAPartialdFext(du_df);
  lq += du_dq.transpose() * lu_condensed;
  lv += du_dv.transpose() * lu_condensed;
  la += du_da.transpose() * lu_condensed;
  lf.head(dimf) += du_df.leftCols(dimf).transpose() * lu_condensed;
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(dim_passive_, max_dimf_);
  Eigen::VectorXd C_res = Eigen::VectorXd::Zero(max_dimc_);
  Cq.topRows(dim_passive_) = dtau_ * du_dq.topRows(dim_passive_);
  Cv.topRows(dim_passive_) = dtau_ * du_dv.topRows(dim_passive_);
  Ca.topRows(dim_passive_) = dtau_ * du_da.topRows(dim_passive_);
  Cf.leftCols(dimf) = dtau_ * du_df.topLeftCorner(dim_passive_, dimf);
  C_res.head(dim_passive_) = dtau_ * u_.head(dim_passive_);
  C_res.head(dim_passive_) += dtau_ * u_res.head(dim_passive_);
  robot_.computeBaumgarteResidual(dim_passive_, dtau_, C_res);
  robot_.computeBaumgarteDerivatives(dim_passive_, dtau_, Cq, Cv, Ca);
  lq += Cq.topRows(dimc).transpose() * mu_.head(dimc);
  lv += Cv.topRows(dimc).transpose() * mu_.head(dimc);
  la += Ca.topRows(dimc).transpose() * mu_.head(dimc);
  lf.head(dimf) += Cf.leftCols(dimf).transpose() * mu_.head(dim_passive_);
  Eigen::MatrixXd Qqq = du_dq.transpose() * luu * du_dq;
  Eigen::MatrixXd Qqv = du_dq.transpose() * luu * du_dv;
  Eigen::MatrixXd Qqa = du_dq.transpose() * luu * du_da;
  Eigen::MatrixXd Qvv = du_dv.transpose() * luu * du_dv;
  Eigen::MatrixXd Qva = du_dv.transpose() * luu * du_da;
  Eigen::MatrixXd Qaa = du_da.transpose() * luu * du_da;
  Eigen::MatrixXd Qvq = Qqv.transpose();
  Eigen::MatrixXd Qqf = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  Eigen::MatrixXd Qvf = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  Eigen::MatrixXd Qaf = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  Eigen::MatrixXd Qff = Eigen::MatrixXd::Zero(max_dimf_, max_dimf_);
  Qqf.leftCols(dimf) = du_dq.transpose() * luu * du_df.leftCols(dimf);
  Qvf.leftCols(dimf) = du_dv.transpose() * luu * du_df.leftCols(dimf);
  Qaf.leftCols(dimf) = du_da.transpose() * luu * du_df.leftCols(dimf);
  Qff.topLeftCorner(dimf, dimf) 
      = du_df.leftCols(dimf).transpose() * luu * du_df.leftCols(dimf);
  joint_space_constraints_ref.condenseSlackAndDual(dtau_, q_, v_, a_, Qqq, Qvv, 
                                                   Qaa, lq, lv, la);
  cost_.augment_lqq(t_, dtau_, q_, v_, a_, Qqq);
  cost_.augment_lvv(t_, dtau_, q_, v_, a_, Qvv);
  cost_.augment_laa(t_, dtau_, q_, v_, a_, Qaa);
  cost_.augment_lff(t_, dtau_, f_, Qff);
  inverter_.setContactStatus(robot_);
  inverter_.precompute(Qaf, Qff);
  Eigen::MatrixXd gen_mat = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pqq_next = gen_mat * gen_mat.transpose();
  const Eigen::MatrixXd Pqv_next = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pvq_next = Eigen::MatrixXd::Random(dimv_, dimv_);
  gen_mat = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd Pvv_next = gen_mat * gen_mat.transpose();
  const Eigen::VectorXd sq_next = Eigen::VectorXd::Random(dimv_);
  const Eigen::VectorXd sv_next = Eigen::VectorXd::Random(dimv_);
  Eigen::MatrixXd Pqq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pqv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::VectorXd sq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd sv = Eigen::VectorXd::Zero(dimv_);
  ocp_.backwardRiccatiRecursion(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next,  
                                sq_next, sv_next, Pqq, Pqv, Pvq, Pvv, sq, sv);
  Eigen::MatrixXd Pqq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pqv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  factorizer_.setIntegrationSensitivities(robot_, dtau_, q_, v_);
  factorizer_.factorize(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next, 
                        Qqq, Qqv, Qvq, Qvv);
  factorizer_.factorize(dtau_, Pqv_next, Pvv_next, Qqa, Qva);
  factorizer_.factorize(dtau_, Pvv_next, Qaa);
  la += dtau_ * Pvq_next * q_res;
  la += dtau_ * Pvv_next * v_res;
  la -= dtau_ * sv_next;
  Eigen::MatrixXd Kaq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Kav = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Kfq = Eigen::MatrixXd::Zero(max_dimf_, dimv_);
  Eigen::MatrixXd Kfv = Eigen::MatrixXd::Zero(max_dimf_, dimv_);
  Eigen::MatrixXd Kmuq = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Kmuv = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd kf = Eigen::VectorXd::Zero(max_dimf_);
  Eigen::VectorXd kmu = Eigen::VectorXd::Zero(max_dimc_);
  if (dimf > 0) {
    inverter_.invert(Qqa, Qva, Qaa, Qqf, Qvf, Cq, Cv, Ca, Cf, la, lf, C_res, 
                     Kaq, Kav, Kfq, Kfv, Kmuq, Kmuv, ka, kf, kmu);
  }
  else {
    inverter_.invert(Qqa, Qva, Qaa, Cq, Cv, Ca, la, C_res, Kaq, Kav, Kmuq, Kmuv, 
                    ka, kmu);
  }
  Pqq_ref = Qqq;
  Pqq_ref += Kaq.transpose() * Qqa.transpose();
  Pqv_ref = Qqv;
  Pqv_ref += Kaq.transpose() * Qva.transpose();
  Pvq_ref = Qvq;
  Pvq_ref += Kav.transpose() * Qqa.transpose();
  Pvv_ref = Qvv;
  Pvv_ref += Kav.transpose() * Qva.transpose();
  Pqq_ref += Kfq.topRows(dimf).transpose() * Qqf.leftCols(dimf).transpose();
  Pqv_ref += Kfq.topRows(dimf).transpose() * Qvf.leftCols(dimf).transpose();
  Pvq_ref += Kfv.topRows(dimf).transpose() * Qqf.leftCols(dimf).transpose();
  Pvv_ref += Kfv.topRows(dimf).transpose() * Qvf.leftCols(dimf).transpose();
  Pqq_ref += Kmuq.topRows(dimc).transpose() * Cq.topRows(dimc);
  Pqv_ref += Kmuq.topRows(dimc).transpose() * Cv.topRows(dimc);
  Pvq_ref += Kmuv.topRows(dimc).transpose() * Cq.topRows(dimc);
  Pvv_ref += Kmuv.topRows(dimc).transpose() * Cv.topRows(dimc);
  Eigen::VectorXd sq_ref = sq_next - lq;
  sq_ref -= Pqq_next * q_res;
  sq_ref -= Pqv_next * v_res;
  sq_ref -= Qqa * ka;
  Eigen::VectorXd sv_ref = dtau_ * sq_next + sv_next - lv;
  sv_ref -= dtau_ * Pqq_next * q_res;
  sv_ref -= Pvq_next * q_res;
  sv_ref -= dtau_ * Pqv_next * v_res;
  sv_ref -= Pvv_next * v_res;
  sv_ref -= Qva * ka;
  sq_ref -= Qqf.leftCols(dimf) * kf.head(dimf);
  sv_ref -= Qvf.leftCols(dimf) * kf.head(dimf);
  sq_ref -= Cq.topRows(dimc).transpose() * kmu.head(dimc);
  sv_ref -= Cv.topRows(dimc).transpose() * kmu.head(dimc);
  EXPECT_TRUE(Pqq.isApprox(Pqq_ref));
  EXPECT_TRUE(Pqv.isApprox(Pqv_ref));
  EXPECT_TRUE(Pvq.isApprox(Pvq_ref));
  EXPECT_TRUE(Pvv.isApprox(Pvv_ref));
  EXPECT_TRUE(sq.isApprox(sq_ref));
  EXPECT_TRUE(sv.isApprox(sv_ref));
  std::cout << "Pqq - Pqq_ref" << std::endl;
  std::cout << Pqq - Pqq_ref << "\n" << std::endl;
  std::cout << "Pqv - Pqv_ref" << std::endl;
  std::cout << Pqv - Pqv_ref  << "\n" << std::endl;
  std::cout << "Pvq - Pvq_ref" << std::endl;
  std::cout << Pvq - Pvq_ref  << "\n" << std::endl;
  std::cout << "Pvv - Pvv_ref" << std::endl;
  std::cout << Pvv - Pvv_ref << "\n" << std::endl;
  std::cout << "sq - sq_ref" << std::endl;
  std::cout << sq - sq_ref  << "\n" << std::endl;
  std::cout << "sv - sv_ref" << std::endl;
  std::cout << sv - sv_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv_);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv_);
  Eigen::VectorXd dq_next = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd dv_next = Eigen::VectorXd::Zero(dimv_);
  ocp_.forwardRiccatiRecursion(dtau_, dq, dv, dq_next, dv_next);
  const Eigen::VectorXd da = ka + Kaq * dq + Kav * dv; 
  const Eigen::VectorXd dq_next_ref = dq + dtau_ * dv + q_res;
  const Eigen::VectorXd dv_next_ref = dv + dtau_ * da + v_res;
  EXPECT_TRUE(dq_next.isApprox(dq_next_ref));
  EXPECT_TRUE(dv_next.isApprox(dv_next_ref));
  ocp_.computeCondensedDirection(dtau_, dq, dv);
  const Eigen::VectorXd df 
      = kf.head(dimf) + Kfq.topRows(dimf) * dq + Kfv.topRows(dimf) * dv;
  const Eigen::VectorXd dmu 
      = kmu.head(dim_passive_) + Kmuq.topRows(dim_passive_) * dq 
                               + Kmuv.topRows(dim_passive_) * dv;
  const Eigen::VectorXd du = u_res + du_dq * dq + du_dv * dv + du_da * da 
                                   + du_df.leftCols(dimf) * df;
  joint_space_constraints_ref.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  const double max_primal_step_size = ocp_.maxPrimalStepSize();
  const double max_dual_step_size = ocp_.maxDualStepSize();
  const double max_primal_step_size_ref 
      = joint_space_constraints_ref.maxSlackStepSize();
  const double max_dual_step_size_ref 
      = joint_space_constraints_ref.maxDualStepSize();
  EXPECT_DOUBLE_EQ(max_primal_step_size, max_primal_step_size_ref);
  EXPECT_DOUBLE_EQ(max_dual_step_size, max_dual_step_size_ref);
  const std::pair<double, double> cost_and_violation 
      = ocp_.costAndConstraintsViolation(robot_, max_primal_step_size, t_, 
                                         dtau_, q_, v_, a_, u_, q_next_, 
                                         v_next_, dq, dv, dq_next, dv_next);
  Eigen::VectorXd q_tmp = q_;
  robot_.integrateConfiguration(dq, max_primal_step_size, q_tmp);
  cost_.setConfigurationJacobian(robot_, q_tmp);
  const Eigen::VectorXd v_tmp = v_ + max_primal_step_size * dv;
  const Eigen::VectorXd a_tmp = a_ + max_primal_step_size * da;
  const Eigen::VectorXd u_tmp = u_ + max_primal_step_size * du;
  Eigen::VectorXd f_tmp = Eigen::VectorXd::Zero(max_dimf_);
  f_tmp.head(dimf) = f_.head(dimf) + max_primal_step_size * df;
  robot_.setContactForces(f_tmp);
  cost_.setContactStatus(robot_);
  double cost_ref = 0;
  cost_ref += cost_.l(t_, dtau_, q_tmp, v_tmp, a_tmp, u_tmp, f_tmp);
  cost_ref += joint_space_constraints_ref.costSlackBarrier(max_primal_step_size);
  Eigen::VectorXd q_res_tmp = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_tmp, q_next_, q_res_tmp);
  q_res_tmp += dtau_ * v_tmp - max_primal_step_size * dq_next;
  const Eigen::VectorXd v_res_tmp 
      = v_tmp - v_next_ + dtau_ * a_tmp - max_primal_step_size * dv_next;
  Eigen::VectorXd u_res_tmp = Eigen::VectorXd::Zero(dimv_);
  robot_.RNEA(q_tmp, v_tmp, a_tmp, u_res_tmp);
  u_res_tmp -= u_tmp;
  C_res.head(dim_passive_) = dtau_ * u_tmp.head(dim_passive_);
  robot_.updateKinematics(q_tmp, v_tmp, a_tmp);
  robot_.computeBaumgarteResidual(dim_passive_, dtau_, C_res);
  const double violation_ref
      = q_res_tmp.lpNorm<1>() + v_res_tmp.lpNorm<1>() 
          + dtau_ * u_res_tmp.lpNorm<1>() 
          + joint_space_constraints_ref.residualL1Nrom(dtau_, q_tmp, v_tmp, a_tmp, u_tmp)
          + C_res.head(dimc).lpNorm<1>();
  EXPECT_DOUBLE_EQ(cost_and_violation.first, cost_ref);
  EXPECT_DOUBLE_EQ(cost_and_violation.second, violation_ref);
  ocp_.updateDual(max_dual_step_size);
  joint_space_constraints_ref.updateDual(max_dual_step_size);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  Eigen::VectorXd q_ref = q_;
  Eigen::VectorXd v_ref = v_;
  Eigen::VectorXd a_ref = a_;
  Eigen::VectorXd u_ref = u_;
  Eigen::VectorXd beta_ref = beta;
  Eigen::VectorXd lmd_ref = lmd_;
  Eigen::VectorXd gmm_ref = gmm_;
  ocp_.updatePrimal(robot_, max_primal_step_size, dtau_, dq, dv, Pqq, Pqv, Pvq, 
                    Pvv, sq, sv, q_, v_, a_, u_, beta, lmd_, gmm_);
  robot_.integrateConfiguration(dq, max_primal_step_size, q_ref);
  v_ref += max_primal_step_size * dv;
  a_ref += max_primal_step_size * da;
  u_ref += max_primal_step_size * du;
  lu -= dtau_ * beta_ref;
  beta_ref += max_primal_step_size * lu / dtau_;
  beta_ref += max_primal_step_size * luu * du / dtau_;
  lmd_ref += max_primal_step_size * (Pqq * dq + Pqv * dv - sq);
  gmm_ref += max_primal_step_size * (Pvq * dq + Pvv * dv - sv);
  joint_space_constraints_ref.updateSlack(max_primal_step_size);
  EXPECT_TRUE(q_.isApprox(q_ref));
  EXPECT_TRUE(v_.isApprox(v_ref));
  EXPECT_TRUE(a_.isApprox(a_ref));
  EXPECT_TRUE(u_.isApprox(u_ref));
  EXPECT_TRUE(beta.isApprox(beta_ref));
  EXPECT_TRUE(lmd_.isApprox(lmd_ref));
  EXPECT_TRUE(gmm_.isApprox(gmm_ref));
}


TEST_F(FloatingBaseSplitOCPTest, KKTError) {
  ASSERT_TRUE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.dim_passive() == 6);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  robot_.generateFeasibleConfiguration(q_next_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<contact_frames_.size(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setActiveContacts(active_contacts);
  ASSERT_TRUE(robot_.max_dimf() == 12);
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  const int dimf = robot_.dimf();
  const int dimc = robot_.dim_passive() + robot_.dimf();
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref.setTimeStep(time_step_);
  joint_space_constraints_ref.setSlackAndDual(dtau_, q_, v_, a_, u_);
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  const double KKT_error 
      = ocp_.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                 beta, lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::VectorXd q_res, v_res;
  q_res = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_, q_next_, q_res);
  q_res.noalias() += dtau_ * v_;
  v_res = v_ + dtau_ * a_ - v_next_;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv_);
  robot_.setContactForces(f_);
  robot_.RNEA(q_, v_, a_, u_res);
  u_res -= u_;
  u_res.array() *= dtau_;
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(max_dimf_);
  cost_.setContactStatus(robot_);
  cost_.setConfigurationJacobian(robot_, q_);
  cost_.lq(t_, dtau_, q_, v_, a_, lq);
  cost_.lv(t_, dtau_, q_, v_, a_, lv);
  cost_.la(t_, dtau_, q_, v_, a_, la);
  cost_.lu(t_, dtau_, u_, lu);
  cost_.lf(t_, dtau_, f_, lf);
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  joint_space_constraints_ref.augmentDualResidual(dtau_, lq, lv, la);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lu);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  robot_.RNEADerivatives(q_, v_, a_, du_dq, du_dv, du_da);
  robot_.updateKinematics(q_, v_, a_);
  robot_.dRNEAPartialdFext(du_df);
  lq += dtau_ * du_dq.transpose() * beta;
  lv += dtau_ * du_dv.transpose() * beta;
  la += dtau_ * du_da.transpose() * beta;
  lu -= dtau_ * beta;
  lf.head(dimf) += dtau_ * du_df.leftCols(dimf).transpose() * beta;
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(dim_passive_, max_dimf_);
  Eigen::VectorXd C_res = Eigen::VectorXd::Zero(max_dimc_);
  Cq.topRows(dim_passive_) = dtau_ * du_dq.topRows(dim_passive_);
  Cv.topRows(dim_passive_) = dtau_ * du_dv.topRows(dim_passive_);
  Ca.topRows(dim_passive_) = dtau_ * du_da.topRows(dim_passive_);
  Cf.leftCols(dimf) = dtau_ * du_df.topLeftCorner(dim_passive_, dimf);
  C_res.head(dim_passive_) = dtau_ * u_.head(dim_passive_);
  robot_.computeBaumgarteResidual(dim_passive_, dtau_, C_res);
  robot_.computeBaumgarteDerivatives(dim_passive_, dtau_, Cq, Cv, Ca);
  lq += Cq.topRows(dimc).transpose() * mu_.head(dimc);
  lv += Cv.topRows(dimc).transpose() * mu_.head(dimc);
  la += Ca.topRows(dimc).transpose() * mu_.head(dimc);
  lf.head(dimf) += Cf.leftCols(dimf).transpose() * mu_.head(dim_passive_);
  double KKT_error_ref = 0;
  KKT_error_ref += q_res.squaredNorm();
  KKT_error_ref += v_res.squaredNorm();
  KKT_error_ref += u_res.squaredNorm();
  KKT_error_ref += lq.squaredNorm();
  KKT_error_ref += lv.squaredNorm();
  KKT_error_ref += la.squaredNorm();
  KKT_error_ref += lu.squaredNorm();
  KKT_error_ref += lf.head(dimf).squaredNorm();
  KKT_error_ref += joint_space_constraints_ref.residualSquaredNrom(dtau_, q_, v_, a_, u_);
  KKT_error_ref += C_res.head(dimc).squaredNorm();
  EXPECT_DOUBLE_EQ(KKT_error, KKT_error_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
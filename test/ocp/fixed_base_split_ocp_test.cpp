#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"
#include "idocp/constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "idocp/manipulator/cost_function.hpp"
#include "idocp/manipulator/constraints.hpp"


namespace idocp {

class FixedBaseSplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    contact_frames_ = {18};
    baum_on_velocity_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    baum_on_position_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    if (rnd()%2==0) {
      robot_ = Robot(urdf_, contact_frames_, baum_on_velocity_, baum_on_position_);
      if (rnd()%2==0) { 
        std::vector<bool> contact_status = {true};
        robot_.setContactStatus(contact_status);
      }
      else {
        std::vector<bool> contact_status = {false};
        robot_.setContactStatus(contact_status);
      }
    }
    else {
      robot_ = Robot(urdf_);
    }
    cost_ = std::make_unique<manipulator::CostFunction>(robot_);
    cost_ref_ = std::make_unique<manipulator::CostFunction>(robot_);
    constraints_ = std::make_unique<manipulator::Constraints>(robot_);
    constraints_ref_ = std::make_unique<manipulator::Constraints>(robot_);
    joint_space_constraints_ref_ = pdipm::JointSpaceConstraints(robot_);
    factorizer_ = RiccatiMatrixFactorizer(robot_);
    inverter_ = RiccatiMatrixInverter(robot_);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    time_step_ = rnd()%10;
    dimq_ = robot_.dimq();
    dimv_ = robot_.dimv();
    dim_passive_ = robot_.dim_passive();
    max_dimf_ = robot_.max_dimf();
    dimf_ = robot_.dimf();
    max_dimc_ = robot_.dim_passive() + robot_.max_dimf();
    dimc_ = robot_.dim_passive() + robot_.dimf();
    q_ = Eigen::VectorXd::Zero(dimq_);
    robot_.generateFeasibleConfiguration(q_);
    v_ = Eigen::VectorXd::Random(dimv_);
    a_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
    f_ = Eigen::VectorXd::Random(max_dimf_);
    lmd_ = Eigen::VectorXd::Random(dimv_);
    gmm_ = Eigen::VectorXd::Random(dimv_);
    mu_ = Eigen::VectorXd::Random(max_dimc_);
    q_next_ = Eigen::VectorXd::Zero(dimq_);
    robot_.generateFeasibleConfiguration(q_next_);
    v_next_ = Eigen::VectorXd::Random(dimv_);
    lmd_next_ = Eigen::VectorXd::Random(dimv_);
    gmm_next_ = Eigen::VectorXd::Random(dimv_);
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  Robot robot_;
  RiccatiMatrixFactorizer factorizer_;
  RiccatiMatrixInverter inverter_;
  std::unique_ptr<CostFunctionInterface> cost_;
  std::unique_ptr<ConstraintsInterface> constraints_;
  std::unique_ptr<CostFunctionInterface> cost_ref_;
  std::unique_ptr<ConstraintsInterface> constraints_ref_;
  pdipm::JointSpaceConstraints joint_space_constraints_ref_;
  double t_, dtau_, baum_on_velocity_, baum_on_position_;
  int time_step_, dimq_, dimv_, dim_passive_, max_dimf_, dimf_, max_dimc_, dimc_;
  std::vector<int> contact_frames_;
  Eigen::VectorXd q_, v_, a_, u_, f_, lmd_, gmm_, mu_, q_next_, v_next_, 
                  lmd_next_, gmm_next_;
};


TEST_F(FixedBaseSplitOCPTest, isFeasible) {
  SplitOCP ocp(robot_, cost_, constraints_);
  EXPECT_EQ(ocp.isFeasible(robot_, q_, v_, a_, u_), 
            joint_space_constraints_ref_.isFeasible(q_, v_, a_, u_));
}


TEST_F(FixedBaseSplitOCPTest, moveConstructor) {
  SplitOCP ocp(robot_, cost_, constraints_);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  const double KKT_error 
      = ocp.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                beta, f_, mu_, lmd_next_, gmm_next_, 
                                q_next_, v_next_);
  SplitOCP ocp_ref(std::move(ocp));
  const double KKT_error_ref
      = ocp_ref.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                    beta, f_, mu_, lmd_next_, gmm_next_, 
                                    q_next_, v_next_);
  EXPECT_DOUBLE_EQ(KKT_error, KKT_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, moveAssignOperator) {
  SplitOCP ocp(robot_, cost_, constraints_);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  const double KKT_error 
      = ocp.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                beta, f_, mu_, lmd_next_, gmm_next_, 
                                q_next_, v_next_);
  SplitOCP ocp_ref;
  ocp_ref = std::move(ocp);
  const double KKT_error_ref
      = ocp_ref.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                    beta, f_, mu_, lmd_next_, gmm_next_, 
                                    q_next_, v_next_);
  EXPECT_DOUBLE_EQ(KKT_error, KKT_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, solveOCP) {
  SplitOCP ocp(robot_, cost_, constraints_);
  ASSERT_FALSE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.dim_passive() == 0);
  while (!ocp.isFeasible(robot_, q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  ASSERT_TRUE(ocp.isFeasible(robot_, q_, v_, a_, u_));
  ocp.initConstraints(robot_, time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref_.setTimeStep(time_step_);
  joint_space_constraints_ref_.setSlackAndDual(dtau_, q_, v_, a_, u_);
  ocp.linearizeOCP(robot_, t_, dtau_, lmd_, gmm_, q_, v_, 
                    a_, u_, f_, mu_, lmd_next_, gmm_next_, q_next_, v_next_);
  if (dimf_ > 0) {
    robot_.updateKinematics(q_, v_, a_);
  }
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(max_dimf_);
  ocplinearizer::linearizeStageCost(robot_, cost_ref_, t_, dtau_, 
                                    q_, v_, a_, u_, f_, 
                                    lq, lv, la, lu, lf);
  Eigen::VectorXd q_res = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd v_res = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv_);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  ocplinearizer::linearizeDynamics(robot_, dtau_, q_, v_, a_, u_, f_, 
                                   q_next_, v_next_, q_res, v_res, u_res, 
                                   du_dq, du_dv, du_da, du_df);
  Eigen::VectorXd C_res = Eigen::VectorXd::Zero(max_dimc_);
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(dim_passive_, max_dimf_);
  ocplinearizer::linearizeConstraints(robot_, dtau_, q_, v_, a_, u_, u_res, 
                                      du_dq, du_dv, du_da, du_df, 
                                      C_res, Cq, Cv, Ca, Cf);
  joint_space_constraints_ref_.augmentDualResidual(dtau_, lu);
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv_, dimv_);
  cost_ref_->luu(robot_, t_, dtau_, u_, luu);
  joint_space_constraints_ref_.condenseSlackAndDual(dtau_, u_, luu, lu);
  Eigen::VectorXd lu_condensed = lu + luu * u_res;
  lq += du_dq.transpose() * lu_condensed;
  lv += du_dv.transpose() * lu_condensed;;
  la += du_da.transpose() * lu_condensed;;
  lf.head(dimf_) += du_df.leftCols(dimf_).transpose() * lu_condensed;
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  joint_space_constraints_ref_.augmentDualResidual(dtau_, lq, lv, la);
  lq += Cq.topRows(dimc_).transpose() * mu_.head(dimc_);
  lv += Cv.topRows(dimc_).transpose() * mu_.head(dimc_);
  la += Ca.topRows(dimc_).transpose() * mu_.head(dimc_);
  lf.head(dimf_) += Cf.leftCols(dimf_).transpose() * mu_.head(dim_passive_);
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
  Qqf = du_dq.transpose() * luu * du_df;
  Qvf = du_dv.transpose() * luu * du_df;
  Qaf = du_da.transpose() * luu * du_df;
  Qff = du_df.transpose() * luu * du_df;
  joint_space_constraints_ref_.condenseSlackAndDual(dtau_, q_, v_, a_, Qqq, Qvv, 
                                                    Qaa, lq, lv, la);
  cost_ref_->augment_lqq(robot_, t_, dtau_, q_, v_, a_, Qqq);
  cost_ref_->augment_lvv(robot_, t_, dtau_, q_, v_, a_, Qvv);
  cost_ref_->augment_laa(robot_, t_, dtau_, q_, v_, a_, Qaa);
  cost_ref_->augment_lff(robot_, t_, dtau_, f_, Qff);
  factorizer_.setIntegrationSensitivities(robot_, dtau_, q_, v_);
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
  ocp.backwardRiccatiRecursion(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next,  
                               sq_next, sv_next, Pqq, Pqv, Pvq, Pvv, sq, sv);
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
  if (dimf_ == 0) {
    inverter_.invert(Qqa, Qva, Qaa, la, Kaq, Kav, ka);
  } 
  else {
    inverter_.invert(Qqa, Qva, Qaa, Qqf, Qvf, Cq, Cv, Ca, la, lf, C_res, 
                     Kaq, Kav, Kfq, Kfv, Kmuq, Kmuv, ka, kf, kmu);
  }
  Eigen::MatrixXd Pqq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pqv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Pvv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::VectorXd sq_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd sv_ref = Eigen::VectorXd::Zero(dimv_);
  Pqq_ref = Qqq;
  Pqq_ref += Kaq.transpose() * Qqa.transpose();
  Pqv_ref = Qqv;
  Pqv_ref += Kaq.transpose() * Qva.transpose();
  Pvq_ref = Qvq;
  Pvq_ref += Kav.transpose() * Qqa.transpose();
  Pvv_ref = Qvv;
  Pvv_ref += Kav.transpose() * Qva.transpose();
  Pqq_ref += Kfq.topRows(dimf_).transpose() * Qqf.leftCols(dimf_).transpose();
  Pqv_ref += Kfq.topRows(dimf_).transpose() * Qvf.leftCols(dimf_).transpose();
  Pvq_ref += Kfv.topRows(dimf_).transpose() * Qqf.leftCols(dimf_).transpose();
  Pvv_ref += Kfv.topRows(dimf_).transpose() * Qvf.leftCols(dimf_).transpose();
  Pqq_ref += Kmuq.topRows(dimc_).transpose() * Cq.topRows(dimc_);
  Pqv_ref += Kmuq.topRows(dimc_).transpose() * Cv.topRows(dimc_);
  Pvq_ref += Kmuv.topRows(dimc_).transpose() * Cq.topRows(dimc_);
  Pvv_ref += Kmuv.topRows(dimc_).transpose() * Cv.topRows(dimc_);
  sq_ref = sq_next - lq;
  sq_ref -= Pqq_next * q_res;
  sq_ref -= Pqv_next * v_res;
  sq_ref -= Qqa * ka;
  sv_ref = dtau_ * sq_next + sv_next - lv;
  sv_ref -= dtau_ * Pqq_next * q_res;
  sv_ref -= Pvq_next * q_res;
  sv_ref -= dtau_ * Pqv_next * v_res;
  sv_ref -= Pvv_next * v_res;
  sv_ref -= Qva * ka;
  sq_ref -= Qqf.leftCols(dimf_) * kf.head(dimf_);
  sv_ref -= Qvf.leftCols(dimf_) * kf.head(dimf_);
  sq_ref -= Cq.topRows(dimc_).transpose() * kmu.head(dimc_);
  sv_ref -= Cv.topRows(dimc_).transpose() * kmu.head(dimc_);
  EXPECT_TRUE(Pqq.isApprox(Pqq_ref));
  EXPECT_TRUE(Pqv.isApprox(Pqv_ref));
  EXPECT_TRUE(Pvq.isApprox(Pvq_ref));;
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
  ocp.forwardRiccatiRecursion(dtau_, dq, dv, dq_next, dv_next);
  const Eigen::VectorXd da = ka + Kaq * dq + Kav * dv; 
  const Eigen::VectorXd dq_next_ref = dq + dtau_ * dv + q_res;
  const Eigen::VectorXd dv_next_ref = dv + dtau_ * da + v_res;
  EXPECT_TRUE(dq_next.isApprox(dq_next_ref));
  EXPECT_TRUE(dv_next.isApprox(dv_next_ref));
  ocp.computeCondensedDirection(dtau_, dq, dv);
  Eigen::VectorXd df = Eigen::VectorXd::Zero(max_dimf_);
  Eigen::VectorXd dmu = Eigen::VectorXd::Zero(max_dimc_);
  df.head(dimf_) = kf.head(dimf_) + Kfq.topRows(dimf_) * dq + Kfv.topRows(dimf_) * dv;
  dmu.head(dimc_) = kmu.head(dimc_) + Kmuq.topRows(dimc_) * dq + Kmuv.topRows(dimc_) * dv;
  Eigen::VectorXd du = u_res + du_dq * dq + du_dv * dv + du_da * da;
  du += du_df.leftCols(dimf_) * df.head(dimf_);
  joint_space_constraints_ref_.computeSlackAndDualDirection(dtau_, dq, dv, da, 
                                                            du);
  const double max_primal_step_size = ocp.maxPrimalStepSize();
  const double max_dual_step_size = ocp.maxDualStepSize();
  const double max_primal_step_size_ref 
      = joint_space_constraints_ref_.maxSlackStepSize();
  const double max_dual_step_size_ref 
      = joint_space_constraints_ref_.maxDualStepSize();
  EXPECT_DOUBLE_EQ(max_primal_step_size, max_primal_step_size_ref);
  EXPECT_DOUBLE_EQ(max_dual_step_size, max_dual_step_size_ref);
  const std::pair<double, double> cost_and_violation 
      = ocp.costAndConstraintsViolation(robot_, max_primal_step_size, t_, 
                                        dtau_, q_, v_, a_, u_, f_, q_next_, 
                                        v_next_, dq, dv, dq_next, dv_next);
  Eigen::VectorXd q_tmp = q_;
  robot_.integrateConfiguration(dq, max_primal_step_size, q_tmp);
  const Eigen::VectorXd v_tmp = v_ + max_primal_step_size * dv;
  const Eigen::VectorXd a_tmp = a_ + max_primal_step_size * da;
  const Eigen::VectorXd u_tmp = u_ + max_primal_step_size * du;
  Eigen::VectorXd f_tmp = f_;
  f_tmp.head(dimf_) += max_primal_step_size * df;
  robot_.setContactForces(f_tmp);
  double cost_ref = 0;
  cost_ref += cost_ref_->l(robot_, t_, dtau_, q_tmp, v_tmp, a_tmp, u_tmp, f_tmp);
  cost_ref += joint_space_constraints_ref_.costSlackBarrier(max_primal_step_size);
  robot_.subtractConfiguration(q_tmp, q_next_, q_res);
  q_res.noalias() += dtau_ * v_tmp - max_primal_step_size * dq_next;
  v_res = v_tmp + dtau_ * a_tmp - v_next_ - max_primal_step_size * dv_next;
  robot_.RNEA(q_tmp, v_tmp, a_tmp, u_res);
  u_res.noalias() -= u_tmp;
  C_res.head(dim_passive_) = dtau_ * u_tmp.head(dim_passive_);
  if (dimf_ > 0) {
    robot_.updateKinematics(q_tmp, v_tmp, a_tmp);
    robot_.computeBaumgarteResidual(dim_passive_, dtau_, C_res);
  }
  double violation_ref = 0;
  violation_ref += q_res.lpNorm<1>();
  violation_ref += v_res.lpNorm<1>();
  violation_ref += dtau_ * u_res.lpNorm<1>();
  violation_ref 
      += joint_space_constraints_ref_.residualL1Nrom(dtau_, q_tmp, v_tmp, 
                                                     a_tmp, u_tmp);
  violation_ref += C_res.head(dimc_).lpNorm<1>();
  EXPECT_DOUBLE_EQ(cost_and_violation.first, cost_ref);
  EXPECT_DOUBLE_EQ(cost_and_violation.second, violation_ref);
  ocp.updateDual(max_dual_step_size);
  joint_space_constraints_ref_.updateDual(max_dual_step_size);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  Eigen::VectorXd q_ref = q_;
  Eigen::VectorXd v_ref = v_;
  Eigen::VectorXd a_ref = a_;
  Eigen::VectorXd u_ref = u_;
  Eigen::VectorXd beta_ref = beta;
  Eigen::VectorXd f_ref = f_;
  Eigen::VectorXd mu_ref = mu_;
  Eigen::VectorXd lmd_ref = lmd_;
  Eigen::VectorXd gmm_ref = gmm_;
  ocp.updatePrimal(robot_, max_primal_step_size, dtau_, dq, dv, Pqq, Pqv, Pvq, 
                   Pvv, sq, sv, q_, v_, a_, u_, beta, f_, mu_, lmd_, gmm_);
  robot_.integrateConfiguration(dq, max_primal_step_size, q_ref);
  v_ref += max_primal_step_size * dv;
  a_ref += max_primal_step_size * da;
  u_ref += max_primal_step_size * du;
  f_ref.head(dimf_) += max_primal_step_size * df.head(dimf_);
  mu_ref.head(dimc_) += max_primal_step_size * dmu.head(dimc_);
  lu -= dtau_ * beta_ref;
  beta_ref += max_primal_step_size * lu / dtau_;
  beta_ref += max_primal_step_size * luu * du / dtau_;
  lmd_ref += max_primal_step_size * (Pqq * dq + Pqv * dv - sq);
  gmm_ref += max_primal_step_size * (Pvq * dq + Pvv * dv - sv);
  joint_space_constraints_ref_.updateSlack(max_primal_step_size);
  EXPECT_TRUE(q_.isApprox(q_ref));
  EXPECT_TRUE(v_.isApprox(v_ref));
  EXPECT_TRUE(a_.isApprox(a_ref));
  EXPECT_TRUE(u_.isApprox(u_ref));
  EXPECT_TRUE(beta.isApprox(beta_ref));
  EXPECT_TRUE(f_.isApprox(f_ref));
  EXPECT_TRUE(mu_.isApprox(mu_ref));
  EXPECT_TRUE(lmd_.isApprox(lmd_ref));
  EXPECT_TRUE(gmm_.isApprox(gmm_ref));
  const double KKT_error 
      = ocp.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                beta, f_, mu_, lmd_next_, gmm_next_, 
                                q_next_, v_next_);
  if (dimf_ > 0) {
    robot_.updateKinematics(q_, v_, a_);
  }
  ocplinearizer::linearizeStageCost(robot_, cost_ref_, t_, dtau_, 
                                    q_, v_, a_, u_, f_, 
                                    lq, lv, la, lu, lf);
  ocplinearizer::linearizeDynamics(robot_, dtau_, q_, v_, a_, u_, f_, 
                                   q_next_, v_next_, q_res, v_res, u_res, 
                                   du_dq, du_dv, du_da, du_df);
  ocplinearizer::linearizeConstraints(robot_, dtau_, q_, v_, a_, u_, u_res, 
                                      du_dq, du_dv, du_da, du_df, 
                                      C_res, Cq, Cv, Ca, Cf);
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  lq += dtau_ * du_dq.transpose() * beta;
  lv += dtau_ * du_dv.transpose() * beta;
  la += dtau_ * du_da.transpose() * beta;
  lu -= dtau_ * beta;
  lf.head(dimf_) += dtau_ * du_df.leftCols(dimf_).transpose() * beta;
  joint_space_constraints_ref_.augmentDualResidual(dtau_, lu);
  joint_space_constraints_ref_.augmentDualResidual(dtau_, lq, lv, la);
  lq += Cq.topRows(dimc_).transpose() * mu_.head(dimc_);
  lv += Cv.topRows(dimc_).transpose() * mu_.head(dimc_);
  la += Ca.topRows(dimc_).transpose() * mu_.head(dimc_);
  lf.head(dimf_) += Cf.leftCols(dimf_).transpose() * mu_.head(dim_passive_);
  double KKT_error_ref = 0;
  KKT_error_ref += q_res.squaredNorm();
  KKT_error_ref += v_res.squaredNorm();
  KKT_error_ref += u_res.squaredNorm();
  KKT_error_ref += lq.squaredNorm();
  KKT_error_ref += lv.squaredNorm();
  KKT_error_ref += la.squaredNorm();
  KKT_error_ref += lu.squaredNorm();
  KKT_error_ref += lf.head(dimf_).squaredNorm();
  KKT_error_ref += joint_space_constraints_ref_.residualSquaredNrom(dtau_, q_, v_, a_, u_);
  KKT_error_ref += C_res.head(dimc_).squaredNorm();
  EXPECT_DOUBLE_EQ(KKT_error, KKT_error_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
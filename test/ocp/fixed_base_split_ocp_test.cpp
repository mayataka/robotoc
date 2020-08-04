#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
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
    robot_ = Robot(urdf_);
    cost_ = manipulator::CostFunction(robot_);
    constraints_ = manipulator::Constraints(robot_);
    ocp_ = SplitOCP(robot_, &cost_, &constraints_);
    factorizer_ = RiccatiMatrixFactorizer(robot_);
    inverter_ = RiccatiMatrixInverter(robot_);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    baum_on_velocity_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    baum_on_position_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    contact_frame_ = 18;
    time_step_ = rnd()%10;
    dimq_ = robot_.dimq();
    dimv_ = robot_.dimv();
    max_dimf_ = robot_.max_dimf();
    dimf_ = 0;
    contact_frames_ = {contact_frame_};
    q_ = Eigen::VectorXd::Random(dimq_);
    v_ = Eigen::VectorXd::Random(dimv_);
    a_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
    f_ = Eigen::VectorXd::Random(max_dimf_);
    lmd_ = Eigen::VectorXd::Random(dimv_);
    gmm_ = Eigen::VectorXd::Random(dimv_);
    mu_ = Eigen::VectorXd::Random(max_dimf_);
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
  manipulator::CostFunction cost_;
  manipulator::Constraints constraints_;
  double t_, dtau_, baum_on_velocity_, baum_on_position_;
  int contact_frame_, time_step_, dimq_, dimv_, max_dimf_, dimf_;
  std::vector<int> contact_frames_;
  Eigen::VectorXd q_, v_, a_, u_, f_, lmd_, gmm_, mu_, q_next_, v_next_, 
                  lmd_next_, gmm_next_;
};


TEST_F(FixedBaseSplitOCPTest, isFeasible) {
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  EXPECT_EQ(ocp_.isFeasible(q_, v_, a_, u_), 
            joint_space_constraints_ref.isFeasible(q_, v_, a_, u_));
}


TEST_F(FixedBaseSplitOCPTest, withoutContacts) {
  ASSERT_FALSE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.max_dimf() == 0);
  ASSERT_TRUE(robot_.dim_passive() == 0);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref.setTimeStep(time_step_);
  joint_space_constraints_ref.setSlackAndDual(dtau_, q_, v_, a_, u_);
  ocp_.linearizeOCP(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                    lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::VectorXd q_res, v_res;
  q_res = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_, q_next_, q_res);
  // q_res = q_ + dtau_ * v_ - q_next_
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
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimq_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
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
  lv += du_dv.transpose() * lu_condensed;;
  la += du_da.transpose() * lu_condensed;;
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
  Qqq += Pqq_next;
  Qqv += dtau_ * Pqq_next;
  Qqv += Pqv_next;
  Qvq += dtau_ * Pqq_next;
  Qvq += Pvq_next;
  Qvv += (dtau_*dtau_) * Pqq_next;
  Qvv += dtau_ * (Pqv_next + Pvq_next);
  Qvv += Pvv_next;
  Qqa += dtau_ * Pqv_next;
  Qva += (dtau_*dtau_) * Pqv_next;
  Qva += dtau_ * Pvv_next;
  Qaa += (dtau_*dtau_) * Pvv_next;
  // factorizer_.factorize(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next, 
  //                       Qqq, Qqv, Qvq, Qvv);
  // factorizer_.factorize(dtau_, Pqv_next, Pvv_next, Qqa, Qva);
  // factorizer_.factorize(dtau_, Pvv_next, Qaa);
  la += dtau_ * Pvq_next * q_res;
  la += dtau_ * Pvv_next * v_res;
  la -= dtau_ * sv_next;
  Eigen::MatrixXd Ginv = Qaa.llt().solve(Eigen::MatrixXd::Identity(dimv_, dimv_));
  Eigen::MatrixXd Kq = - Ginv * Qqa.transpose();
  Eigen::MatrixXd Kv = - Ginv * Qva.transpose();
  Eigen::VectorXd ka = - Ginv * la;
  // Eigen::MatrixXd Kq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  // Eigen::MatrixXd Kv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  // Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv_);
  // inverter_.invert(Qqa, Qva, Qaa, la, Kq, Kv, ka);
  Pqq_ref = Qqq;
  Pqq_ref += Kq.transpose() * Qqa.transpose();
  Pqv_ref = Qqv;
  Pqv_ref += Kq.transpose() * Qva.transpose();
  Pvq_ref = Qvq;
  Pvq_ref += Kv.transpose() * Qqa.transpose();
  Pvv_ref = Qvv;
  Pvv_ref += Kv.transpose() * Qva.transpose();
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
  const Eigen::VectorXd da = ka + Kq * dq + Kv * dv; 
  const Eigen::VectorXd dq_next_ref = dq + dtau_ * dv + q_res;
  const Eigen::VectorXd dv_next_ref = dv + dtau_ * da + v_res;
  EXPECT_TRUE(dq_next.isApprox(dq_next_ref));
  EXPECT_TRUE(dv_next.isApprox(dv_next_ref));
  ocp_.computeCondensedDirection(dtau_, dq, dv);
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
  const Eigen::VectorXd q_tmp = q_ + max_primal_step_size * dq;
  const Eigen::VectorXd v_tmp = v_ + max_primal_step_size * dv;
  const Eigen::VectorXd a_tmp = a_ + max_primal_step_size * da;
  const Eigen::VectorXd u_tmp = u_ + max_primal_step_size * du;
  const double cost_ref 
      = cost_.l(t_, dtau_, q_tmp, v_tmp, a_tmp, u_tmp, f_) 
        + joint_space_constraints_ref.costSlackBarrier(max_primal_step_size);
  const Eigen::VectorXd q_res_tmp 
      = q_tmp - q_next_ + dtau_ * v_tmp - max_primal_step_size * dq_next;
  const Eigen::VectorXd v_res_tmp 
      = v_tmp - v_next_ + dtau_ * a_tmp - max_primal_step_size * dv_next;
  Eigen::VectorXd u_res_tmp = Eigen::VectorXd::Zero(dimv_);
  robot_.RNEA(q_tmp, v_tmp, a_tmp, u_res_tmp);
  u_res_tmp -= u_tmp;
  const double violation_ref
      = q_res_tmp.lpNorm<1>() + v_res_tmp.lpNorm<1>() 
          + dtau_ * u_res_tmp.lpNorm<1>() 
          + joint_space_constraints_ref.residualL1Nrom(dtau_, q_tmp, v_tmp, a_tmp, u_tmp);
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
  q_ref += max_primal_step_size * dq;
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


TEST_F(FixedBaseSplitOCPTest, withoutActiveContacts) {
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
  ASSERT_FALSE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.max_dimf() == 0);
  ASSERT_TRUE(robot_.dim_passive() == 0);
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
  robot_ = Robot(urdf_, contact_frames_, baum_on_velocity_, baum_on_position_);
  ASSERT_FALSE(robot_.has_floating_base());
  ASSERT_EQ(robot_.max_dimf(), 3);
  ASSERT_EQ(robot_.dimf(), 0);
  f_ = Eigen::VectorXd::Random(robot_.max_dimf());
  mu_ = Eigen::VectorXd::Random(robot_.max_dimf());
  cost_ = manipulator::CostFunction(robot_);
  constraints_ = manipulator::Constraints(robot_);
  ocp_ = SplitOCP(robot_, &cost_, &constraints_);
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
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


TEST_F(FixedBaseSplitOCPTest, withActiveContacts) {
  robot_ = Robot(urdf_, contact_frames_, baum_on_velocity_, baum_on_position_);
  std::vector<bool> active_constacts = {true};
  robot_.setActiveContacts(active_constacts);
  ASSERT_FALSE(robot_.has_floating_base());
  ASSERT_EQ(robot_.max_dimf(), 3);
  ASSERT_EQ(robot_.dimf(), 3);
  dimf_ = robot_.dimf();
  f_ = Eigen::VectorXd::Random(robot_.max_dimf());
  mu_ = Eigen::VectorXd::Random(robot_.max_dimf());
  cost_ = manipulator::CostFunction(robot_);
  constraints_ = manipulator::Constraints(robot_);
  factorizer_ = RiccatiMatrixFactorizer(robot_);
  inverter_ = RiccatiMatrixInverter(robot_);
  ocp_ = SplitOCP(robot_, &cost_, &constraints_);
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref.setTimeStep(time_step_);
  joint_space_constraints_ref.setSlackAndDual(dtau_, q_, v_, a_, u_);
  ocp_.linearizeOCP(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                    lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::VectorXd q_res, v_res;
  q_res = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_, q_next_, q_res);
  // q_res = q_ + dtau_ * v_ - q_next_
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
  cost_.setContactStatus(robot_);
  robot_.RNEA(q_, v_, a_, u_res);
  u_res -= u_;
  const Eigen::VectorXd lu_condensed = lu + luu * u_res;
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimq_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf_);
  cost_.lq(t_, dtau_, q_, v_, a_, lq);
  cost_.lv(t_, dtau_, q_, v_, a_, lv);
  cost_.la(t_, dtau_, q_, v_, a_, la);
  cost_.lf(t_, dtau_, f_, lf);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lq, lv, la);
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv_, dimf_);
  robot_.RNEADerivatives(q_, v_, a_, du_dq, du_dv, du_da);
  robot_.updateKinematics(q_, v_, a_);
  robot_.dRNEAPartialdFext(du_df);
  lq += du_dq.transpose() * lu_condensed;
  lv += du_dv.transpose() * lu_condensed;
  la += du_da.transpose() * lu_condensed;
  lf += du_df.transpose() * lu_condensed;
  Eigen::VectorXd C_res = Eigen::VectorXd::Zero(dimf_);
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(dimf_, dimv_);
  robot_.computeBaumgarteResidual(0, dtau_, C_res);
  robot_.computeBaumgarteDerivatives(0, dtau_, Cq, Cv, Ca);
  lq += Cq.transpose() * mu_;
  lv += Cv.transpose() * mu_;
  la += Ca.transpose() * mu_;
  Eigen::MatrixXd Qqq = du_dq.transpose() * luu * du_dq;
  Eigen::MatrixXd Qqv = du_dq.transpose() * luu * du_dv;
  Eigen::MatrixXd Qqa = du_dq.transpose() * luu * du_da;
  Eigen::MatrixXd Qvv = du_dv.transpose() * luu * du_dv;
  Eigen::MatrixXd Qva = du_dv.transpose() * luu * du_da;
  Eigen::MatrixXd Qaa = du_da.transpose() * luu * du_da;
  Eigen::MatrixXd Qqf = du_dq.transpose() * luu * du_df;
  Eigen::MatrixXd Qvf = du_dv.transpose() * luu * du_df;
  Eigen::MatrixXd Qaf = du_da.transpose() * luu * du_df;
  Eigen::MatrixXd Qff = du_df.transpose() * luu * du_df;
  joint_space_constraints_ref.condenseSlackAndDual(dtau_, q_, v_, a_, Qqq, Qvv, 
                                                   Qaa, lq, lv, la);
  cost_.augment_lqq(t_, dtau_, q_, v_, a_, Qqq);
  cost_.augment_lvv(t_, dtau_, q_, v_, a_, Qvv);
  cost_.augment_laa(t_, dtau_, q_, v_, a_, Qaa);
  Eigen::MatrixXd Qvq = Qqv.transpose();
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
  Qqq += Pqq_next;
  Qqv += dtau_ * Pqq_next;
  Qqv += Pqv_next;
  Qvq += dtau_ * Pqq_next;
  Qvq += Pvq_next;
  Qvv += (dtau_*dtau_) * Pqq_next;
  Qvv += dtau_ * (Pqv_next + Pvq_next);
  Qvv += Pvv_next;
  Qqa += dtau_ * Pqv_next;
  Qva += (dtau_*dtau_) * Pqv_next;
  Qva += dtau_ * Pvv_next;
  Qaa += (dtau_*dtau_) * Pvv_next;
  // factorizer_.factorize(dtau_, Pqq_next, Pqv_next, Pvq_next, Pvv_next, 
  //                       Qqq, Qqv, Qvq, Qvv);
  // factorizer_.factorize(dtau_, Pqv_next, Pvv_next, Qqa, Qva);
  // factorizer_.factorize(dtau_, Pvv_next, Qaa);
  la += dtau_ * Pvq_next * q_res;
  la += dtau_ * Pvv_next * v_res;
  la -= dtau_ * sv_next;
  Eigen::MatrixXd Kaq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Kav = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd Kfq = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::MatrixXd Kfv = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::MatrixXd Kmuq = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::MatrixXd Kmuv = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd kf = Eigen::VectorXd::Zero(dimf_);
  Eigen::VectorXd kmu = Eigen::VectorXd::Zero(dimf_);
  inverter_.invert(Qqa, Qva, Qaa, Qqf, Qvf, Cq, Cv, Ca, la, lf, C_res, 
                   Kaq, Kav, Kfq, Kfv, Kmuq, Kmuv, ka, kf, kmu);
  Pqq_ref = Qqq;
  Pqq_ref += Kaq.transpose() * Qqa.transpose();
  Pqv_ref = Qqv;
  Pqv_ref += Kaq.transpose() * Qva.transpose();
  Pvq_ref = Qvq;
  Pvq_ref += Kav.transpose() * Qqa.transpose();
  Pvv_ref = Qvv;
  Pvv_ref += Kav.transpose() * Qva.transpose();
  Pqq_ref += Kfq.transpose() * Qqf.transpose();
  Pqv_ref += Kfq.transpose() * Qvf.transpose();
  Pvq_ref += Kfv.transpose() * Qqf.transpose();
  Pvv_ref += Kfv.transpose() * Qvf.transpose();
  Pqq_ref += Kmuq.transpose() * Cq;
  Pqv_ref += Kmuq.transpose() * Cv;
  Pvq_ref += Kmuv.transpose() * Cq;
  Pvv_ref += Kmuv.transpose() * Cv;
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
  sq_ref -= Qqf * kf;
  sv_ref -= Qvf * kf;
  sq_ref -= Cq.transpose() * kmu;
  sv_ref -= Cv.transpose() * kmu;
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
  const Eigen::VectorXd df = kf + Kfq * dq + Kfv * dv; 
  const Eigen::VectorXd dmu = kmu + Kmuq * dq + Kmuv * dv; 
  Eigen::VectorXd du = u_res + du_dq * dq + du_dv * dv + du_da * da + du_df * df;
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
  const Eigen::VectorXd q_tmp = q_ + max_primal_step_size * dq;
  const Eigen::VectorXd v_tmp = v_ + max_primal_step_size * dv;
  const Eigen::VectorXd a_tmp = a_ + max_primal_step_size * da;
  const Eigen::VectorXd u_tmp = u_ + max_primal_step_size * du;
  const Eigen::VectorXd f_tmp = f_ + max_primal_step_size * df;
  robot_.setContactForces(f_tmp);
  const double cost_ref 
      = cost_.l(t_, dtau_, q_tmp, v_tmp, a_tmp, u_tmp, f_tmp) 
        + joint_space_constraints_ref.costSlackBarrier(max_primal_step_size);
  const Eigen::VectorXd q_res_tmp 
      = q_tmp - q_next_ + dtau_ * v_tmp - max_primal_step_size * dq_next;
  const Eigen::VectorXd v_res_tmp 
      = v_tmp - v_next_ + dtau_ * a_tmp - max_primal_step_size * dv_next;
  Eigen::VectorXd u_res_tmp = Eigen::VectorXd::Zero(dimv_);
  robot_.RNEA(q_tmp, v_tmp, a_tmp, u_res_tmp);
  u_res_tmp -= u_tmp;
  robot_.updateKinematics(q_tmp, v_tmp, a_tmp);
  robot_.computeBaumgarteResidual(0, dtau_, C_res);
  const double violation_ref
      = q_res_tmp.lpNorm<1>() + v_res_tmp.lpNorm<1>() 
          + dtau_ * u_res_tmp.lpNorm<1>() 
          + joint_space_constraints_ref.residualL1Nrom(dtau_, q_tmp, v_tmp, a_tmp, u_tmp)
          + C_res.lpNorm<1>();
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
  q_ref += max_primal_step_size * dq;
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


TEST_F(FixedBaseSplitOCPTest, KKTErrorWithoutContacts) {
  ASSERT_FALSE(robot_.has_floating_base());
  ASSERT_TRUE(robot_.max_dimf() == 0);
  ASSERT_TRUE(robot_.dim_passive() == 0);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref.setTimeStep(time_step_);
  joint_space_constraints_ref.setSlackAndDual(dtau_, q_, v_, a_, u_);
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  const double KKT_error 
      = ocp_.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                 beta, lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::VectorXd q_res = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_, q_next_, q_res);
  q_res.noalias() += dtau_ * v_;
  Eigen::VectorXd v_res = v_ + dtau_ * a_ - v_next_;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv_);
  robot_.RNEA(q_, v_, a_, u_res);
  u_res -= u_;
  u_res = dtau_ * u_res;
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimq_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv_);
  cost_.lq(t_, dtau_, q_, v_, a_, lq);
  cost_.lv(t_, dtau_, q_, v_, a_, lv);
  cost_.la(t_, dtau_, q_, v_, a_, la);
  cost_.lu(t_, dtau_, u_, lu);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lu);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lq, lv, la);
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  robot_.RNEADerivatives(q_, v_, a_, du_dq, du_dv, du_da);
  lq += dtau_ * du_dq.transpose() * beta;
  lv += dtau_ * du_dv.transpose() * beta;
  la += dtau_ * du_da.transpose() * beta;
  lu -= dtau_ * beta;
  double KKT_error_ref = 0;
  KKT_error_ref += q_res.squaredNorm();
  KKT_error_ref += v_res.squaredNorm();
  KKT_error_ref += u_res.squaredNorm();
  KKT_error_ref += lq.squaredNorm();
  KKT_error_ref += lv.squaredNorm();
  KKT_error_ref += la.squaredNorm();
  KKT_error_ref += lu.squaredNorm();
  KKT_error_ref 
      += joint_space_constraints_ref.residualSquaredNrom(dtau_, q_, v_, a_, u_);
  EXPECT_DOUBLE_EQ(KKT_error_ref, KKT_error);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorWitContacts) {
  robot_ = Robot(urdf_, contact_frames_, baum_on_velocity_, baum_on_position_);
  std::vector<bool> active_constacts = {true};
  robot_.setActiveContacts(active_constacts);
  ASSERT_FALSE(robot_.has_floating_base());
  ASSERT_EQ(robot_.max_dimf(), 3);
  ASSERT_EQ(robot_.dimf(), 3);
  dimf_ = robot_.dimf();
  f_ = Eigen::VectorXd::Random(robot_.max_dimf());
  mu_ = Eigen::VectorXd::Random(robot_.max_dimf());
  cost_ = manipulator::CostFunction(robot_);
  constraints_ = manipulator::Constraints(robot_);
  factorizer_ = RiccatiMatrixFactorizer(robot_);
  inverter_ = RiccatiMatrixInverter(robot_);
  ocp_ = SplitOCP(robot_, &cost_, &constraints_);
  ocp_.set_f(f_);
  ocp_.set_mu(mu_);
  pdipm::JointSpaceConstraints joint_space_constraints_ref(robot_);
  robot_.generateFeasibleConfiguration(q_);
  while (!ocp_.isFeasible(q_, v_, a_, u_)) {
    v_ = Eigen::VectorXd::Random(dimv_);
    u_ = Eigen::VectorXd::Random(dimv_);
  }
  ASSERT_TRUE(ocp_.isFeasible(q_, v_, a_, u_));
  ocp_.initConstraints(time_step_, dtau_, q_, v_, a_, u_);
  joint_space_constraints_ref.setTimeStep(time_step_);
  joint_space_constraints_ref.setSlackAndDual(dtau_, q_, v_, a_, u_);
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(dimv_);
  robot_.setContactForces(Eigen::VectorXd::Zero(dimf_));
  const double KKT_error 
      = ocp_.squaredKKTErrorNorm(robot_, t_, dtau_, lmd_, gmm_, q_, v_, a_, u_, 
                                 beta, lmd_next_, gmm_next_, q_next_, v_next_);
  Eigen::VectorXd q_res = Eigen::VectorXd::Zero(dimv_);
  robot_.differenceConfiguration(q_, q_next_, q_res);
  q_res.noalias() += dtau_ * v_;
  Eigen::VectorXd v_res = v_ + dtau_ * a_ - v_next_;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv_);
  robot_.setContactForces(f_);
  robot_.RNEA(q_, v_, a_, u_res);
  u_res -= u_;
  u_res = dtau_ * u_res;
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimq_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf_);
  cost_.lq(t_, dtau_, q_, v_, a_, lq);
  cost_.lv(t_, dtau_, q_, v_, a_, lv);
  cost_.la(t_, dtau_, q_, v_, a_, la);
  cost_.lu(t_, dtau_, u_, lu);
  cost_.setContactStatus(robot_);
  cost_.lf(t_, dtau_, f_, lf);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lu);
  joint_space_constraints_ref.augmentDualResidual(dtau_, lq, lv, la);
  lq += lmd_next_ - lmd_;
  lv += dtau_ * lmd_next_ + gmm_next_ - gmm_;
  la += dtau_ * gmm_next_;
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  robot_.RNEADerivatives(q_, v_, a_, du_dq, du_dv, du_da);
  robot_.updateKinematics(q_, v_, a_);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv_, dimf_);
  robot_.dRNEAPartialdFext(du_df);
  lq += dtau_ * du_dq.transpose() * beta;
  lv += dtau_ * du_dv.transpose() * beta;
  la += dtau_ * du_da.transpose() * beta;
  lu -= dtau_ * beta;
  lf += dtau_ * du_df.transpose() * beta;
  Eigen::VectorXd C_res = Eigen::VectorXd::Zero(dimf_);
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(dimf_, dimv_);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(dimf_, dimv_);
  robot_.computeBaumgarteResidual(0, dtau_, C_res);
  robot_.computeBaumgarteDerivatives(0, dtau_, Cq, Cv, Ca);
  lq += Cq.transpose() * mu_;
  lv += Cv.transpose() * mu_;
  la += Ca.transpose() * mu_;
  double KKT_error_ref = 0;
  KKT_error_ref += q_res.squaredNorm();
  KKT_error_ref += v_res.squaredNorm();
  KKT_error_ref += u_res.squaredNorm();
  KKT_error_ref += lq.squaredNorm();
  KKT_error_ref += lv.squaredNorm();
  KKT_error_ref += la.squaredNorm();
  KKT_error_ref += lu.squaredNorm();
  KKT_error_ref += lf.squaredNorm();
  KKT_error_ref += C_res.squaredNorm();
  KKT_error_ref 
      += joint_space_constraints_ref.residualSquaredNrom(dtau_, q_, v_, a_, u_);
  EXPECT_DOUBLE_EQ(KKT_error_ref, KKT_error);
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "idocp/manipulator/cost_function.hpp"
#include "idocp/manipulator/constraints.hpp"


namespace idocp {

class FixedBaseOCPLinearizerTest : public ::testing::Test {
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
    cost_ = std::make_shared<manipulator::CostFunction>(robot_);
    constraints_ = std::make_shared<manipulator::Constraints>(robot_);
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
  std::shared_ptr<CostFunctionInterface> cost_;
  std::shared_ptr<ConstraintsInterface> constraints_;
  double t_, dtau_, baum_on_velocity_, baum_on_position_;
  int time_step_, dimq_, dimv_, dim_passive_, max_dimf_, dimf_, max_dimc_, dimc_;
  std::vector<int> contact_frames_;
  Eigen::VectorXd q_, v_, a_, u_, f_, lmd_, gmm_, mu_, q_next_, v_next_, 
                  lmd_next_, gmm_next_;
};


TEST_F(FixedBaseOCPLinearizerTest, linearizeStageCost) {
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(max_dimf_);
  ocplinearizer::linearizeStageCost(robot_, cost_, t_, dtau_, q_, v_, a_, u_, f_,
                                    lq, lv, la, lu, lf);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd la_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lu_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(max_dimf_);
  cost_->lq(robot_, t_, dtau_, q_, v_, a_, lq_ref);
  cost_->lv(robot_, t_, dtau_, q_, v_, a_, lv_ref);
  cost_->la(robot_, t_, dtau_, q_, v_, a_, la_ref);
  cost_->lu(robot_, t_, dtau_, u_, lu_ref);
  cost_->lf(robot_, t_, dtau_, f_, lf_ref);
  EXPECT_TRUE(lq.isApprox(lq_ref));
  EXPECT_TRUE(lv.isApprox(lv_ref));
  EXPECT_TRUE(la.isApprox(la_ref));
  EXPECT_TRUE(lu.isApprox(lu_ref));
  EXPECT_TRUE(lf.isApprox(lf_ref));
}


TEST_F(FixedBaseOCPLinearizerTest, linearizeDynamics) {
  Eigen::VectorXd q_res = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd v_res = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv_);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  if (dimf_ > 0) {
    robot_.updateKinematics(q_, v_, a_);
  }
  ocplinearizer::linearizeDynamics(robot_, dtau_, q_, v_, a_, u_, f_, 
                                   q_next_, v_next_, q_res, v_res, u_res, 
                                   du_dq, du_dv, du_da, du_df);
  Eigen::VectorXd q_res_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd v_res_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd u_res_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::MatrixXd du_dq_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_dv_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_da_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd du_df_ref = Eigen::MatrixXd::Zero(dimv_, max_dimf_);
  if (dimf_ > 0) {
    robot_.updateKinematics(q_, v_, a_);
  }
  robot_.subtractConfiguration(q_, q_next_, q_res_ref);
  q_res_ref += dtau_ * v_;
  v_res_ref = v_ + dtau_ * a_ - v_next_;
  robot_.setContactForces(f_);
  robot_.RNEA(q_, v_, a_, u_res_ref);
  u_res_ref -= u_;
  robot_.RNEADerivatives(q_, v_, a_, du_dq_ref, du_dv_ref, du_da_ref);
  robot_.dRNEAPartialdFext(du_df_ref);
  EXPECT_TRUE(q_res.isApprox(q_res_ref));
  EXPECT_TRUE(v_res.isApprox(v_res_ref));
  EXPECT_TRUE(u_res.isApprox(u_res_ref));
  EXPECT_TRUE(du_dq.isApprox(du_dq_ref));
  EXPECT_TRUE(du_dv.isApprox(du_dv_ref));
  EXPECT_TRUE(du_da.isApprox(du_da_ref));
  EXPECT_TRUE(du_df.isApprox(du_df_ref));
}


TEST_F(FixedBaseOCPLinearizerTest, linearizeConstraints) {
  const Eigen::VectorXd u_res = Eigen::VectorXd::Random(dimv_);
  const Eigen::MatrixXd du_dq = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd du_dv = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd du_da = Eigen::MatrixXd::Random(dimv_, dimv_);
  const Eigen::MatrixXd du_df = Eigen::MatrixXd::Random(dimv_, max_dimf_);
  Eigen::VectorXd C_res = Eigen::VectorXd::Zero(max_dimc_);
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(dim_passive_, max_dimf_);
  if (dimf_ > 0) {
    robot_.updateKinematics(q_, v_, a_);
  }
  ocplinearizer::linearizeConstraints(robot_, dtau_, q_, v_, a_, u_, u_res,
                                      du_dq, du_dv, du_da, du_df, C_res,
                                      Cq, Cv, Ca, Cf);
  Eigen::VectorXd C_res_ref = Eigen::VectorXd::Zero(max_dimc_);
  Eigen::MatrixXd Cq_ref = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cv_ref = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Ca_ref = Eigen::MatrixXd::Zero(max_dimc_, dimv_);
  Eigen::MatrixXd Cf_ref = Eigen::MatrixXd::Zero(dim_passive_, max_dimf_);
  if (dimf_ > 0) {
    robot_.updateKinematics(q_, v_, a_);
  } 
  C_res_ref.head(dim_passive_) = dtau_ * (u_.head(dim_passive_)+u_res.head(dim_passive_));
  Cq_ref.topRows(dim_passive_) = dtau_ * du_dq.topRows(dim_passive_);
  Cv_ref.topRows(dim_passive_) = dtau_ * du_dv.topRows(dim_passive_);
  Ca_ref.topRows(dim_passive_) = dtau_ * du_da.topRows(dim_passive_);
  Cf_ref.leftCols(dimf_) = dtau_ * du_df.topLeftCorner(dim_passive_, dimf_);
  robot_.computeBaumgarteResidual(dim_passive_, dtau_, C_res_ref);
  robot_.computeBaumgarteDerivatives(dim_passive_, dtau_, Cq_ref, Cv_ref, Ca_ref);
  EXPECT_TRUE(C_res.isApprox(C_res_ref));
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cf.isApprox(Cf_ref));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
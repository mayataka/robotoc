#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "robot/robot.hpp"
#include "ocp/split_ocp.hpp"
#include "manipulator/cost_function.hpp"
#include "manipulator/constraints.hpp"
#include "quadruped/cost_function.hpp"
#include "quadruped/constraints.hpp"


namespace idocp {

class SplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    manipulator_urdf_ = "../../urdf/iiwa14/iiwa14.urdf";
    quadruped_urdf_ = "../../urdf/anymal/anymal.urdf";
    manipulator_ = Robot(manipulator_urdf_);
    quadruped_ = Robot(quadruped_urdf_);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double t_, dtau_;
  std::string manipulator_urdf_, quadruped_urdf_;
  Robot manipulator_, quadruped_;
};


TEST_F(SplitOCPTest, isFeasible) {
  manipulator::CostFunction cost(manipulator_);
  manipulator::Constraints constraints(manipulator_);
  SplitOCP ocp(manipulator_, &cost, &constraints);
  const Eigen::VectorXd q = Eigen::VectorXd::Random(manipulator_.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(manipulator_.dimv());
  EXPECT_TRUE(ocp.isFeasible(manipulator_, q, v));
  std::random_device rnd;
  ocp.initConstraints(manipulator_, rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(manipulator_, q, v));
}


TEST_F(SplitOCPTest, linearizeOCPFixedBase) {
  manipulator::CostFunction cost(manipulator_);
  manipulator::Constraints constraints(manipulator_);
  SplitOCP ocp(manipulator_, &cost, &constraints);
  const int dimq = manipulator_.dimq();
  const int dimv = manipulator_.dimv();
  const Eigen::VectorXd q = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  Eigen::MatrixXd Qqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv = Eigen::VectorXd::Zero(dimv);
  std::random_device rnd;
  ocp.initConstraints(manipulator_, rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(manipulator_, q, v));
  ocp.linearizeOCP(manipulator_, t_, lmd, gmm, q, v, Qqq, Qqv, Qvq, Qvv, Qq, Qv);

  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv_ref = Eigen::VectorXd::Zero(dimv);
  cost.phiq(t_, q, v, Qq_ref);
  cost.phiv(t_, q, v, Qv_ref);
  Qq_ref = - Qq_ref + lmd;
  Qv_ref = - Qv_ref + gmm;
  cost.phiqq(t_, q, v, Qqq_ref);
  cost.phivv(t_, q, v, Qvv_ref);
  EXPECT_TRUE(Qqq.isApprox(Qqq_ref));
  EXPECT_TRUE(Qqv.isApprox(Qqv_ref));
  EXPECT_TRUE(Qvq.isApprox(Qvq_ref));
  EXPECT_TRUE(Qvv.isApprox(Qvv_ref));
  EXPECT_TRUE(Qq.isApprox(Qq_ref));
  EXPECT_TRUE(Qv.isApprox(Qv_ref));
}


TEST_F(SplitOCPTest, linearizeOCPFloatingBase) {
  quadruped::CostFunction cost(quadruped_);
  quadruped::Constraints constraints(quadruped_);
  SplitOCP ocp(quadruped_, &cost, &constraints);
  const int dimq = quadruped_.dimq();
  const int dimv = quadruped_.dimv();
  pinocchio::Model model;
  pinocchio::urdf::buildModel(quadruped_urdf_, model);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(dimq);
  quadruped_.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(dimv);
  Eigen::MatrixXd Qqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv = Eigen::VectorXd::Zero(dimv);
  std::random_device rnd;
  ocp.initConstraints(quadruped_, rnd()%50, dtau_, q, v);
  EXPECT_TRUE(ocp.isFeasible(quadruped_, q, v));
  ocp.linearizeOCP(quadruped_, t_, lmd, gmm, q, v, Qqq, Qqv, Qvq, Qvv, Qq, Qv);

  cost.setConfigurationJacobian(quadruped_, q);
  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qqv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Qvv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd Qq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd Qv_ref = Eigen::VectorXd::Zero(dimv);
  cost.phiq(t_, q, v, Qq_ref);
  cost.phiv(t_, q, v, Qv_ref);
  Qq_ref = - Qq_ref + lmd;
  Qv_ref = - Qv_ref + gmm;
  cost.phiqq(t_, q, v, Qqq_ref);
  cost.phivv(t_, q, v, Qvv_ref);
  EXPECT_TRUE(Qqq.isApprox(Qqq_ref));
  EXPECT_TRUE(Qqv.isApprox(Qqv_ref));
  EXPECT_TRUE(Qvq.isApprox(Qvq_ref));
  EXPECT_TRUE(Qvv.isApprox(Qvv_ref));
  EXPECT_TRUE(Qq.isApprox(Qq_ref));
  EXPECT_TRUE(Qv.isApprox(Qv_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
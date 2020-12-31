#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_state_equation.hpp"


namespace idocp {

class ImpulseStateEquationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(ImpulseStateEquationTest, forwardEulerFixedBase) {
  Robot robot(fixed_base_urdf_);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(ImpulseStateEquationTest, forwardEulerFloatingBase) {
  Robot robot(floating_base_urdf_);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq.transpose()*s_next.lmd+dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dq));
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(ImpulseStateEquationTest, backwardEulerFixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next,
                                               kkt_matrix,  kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+s.dv)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::LinearizeImpulseBackwardEulerTerminal(robot, q_prev, v_prev, s,
                                                       kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+s.dv)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((-1*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(ImpulseStateEquationTest, backwardEulerFloatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next,
                                               kkt_matrix,  kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  Eigen::MatrixXd dsubtract_dqminus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dqplus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dqminus);
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dqplus);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+s.dv)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqplus.transpose()*s_next.lmd+dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::LinearizeImpulseBackwardEulerTerminal(robot, q_prev, v_prev, s,
                                                       kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+s.dv)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/state_equation.hpp"


namespace idocp {

class StateEquationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
};


TEST_F(StateEquationTest, forwardEulerFixedbase) {
  Robot robot(fixed_base_urdf);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q+dtau*s.v-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dtau*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dtau*s_next.gmm)));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(StateEquationTest, forwardEulerFloatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dtau*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq.transpose()*s_next.lmd+dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dtau*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dq));
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
  std::cout << "kkt_matrix.Fqq()" << std::endl;
  std::cout << kkt_matrix.Fqq() << std::endl;
  std::cout << "kkt_matrix.Fqq_prev" << std::endl;
  std::cout << kkt_matrix.Fqq_prev << std::endl;
}


TEST_F(StateEquationTest, backwardEulerFixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, s_next,
                                        kkt_matrix,  kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau*s.lmd-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dtau*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, 
                                                s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((-1*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dtau*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(StateEquationTest, backwardEulerFloatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::LinearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, s_next,
                                        kkt_matrix,  kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  Eigen::MatrixXd dsubtract_dqminus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dqplus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dqminus);
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dqplus);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqplus.transpose()*s_next.lmd+dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau*s.lmd-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dtau*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, 
                                                s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dtau*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  std::cout << "kkt_matrix.Fqq()" << std::endl;
  std::cout << kkt_matrix.Fqq() << std::endl;
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
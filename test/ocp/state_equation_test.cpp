#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/state_equation.hpp"


namespace idocp {

class StateEquationTest : public ::testing::Test {
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


TEST_F(StateEquationTest, forwardEuler_fixed_base) {
  Robot robot(fixed_base_urdf_);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  stateequation::LinearizeForwardEuler(robot, dtau_, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q+dtau_*s.v-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dtau_*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s_next.gmm)));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(StateEquationTest, forwardEuler_floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  stateequation::LinearizeForwardEuler(robot, dtau_, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dtau_*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq.transpose()*s_next.lmd+dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq.isApprox(dsubtract_dq));
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
  std::cout << "kkt_matrix.Fqq" << std::endl;
  std::cout << kkt_matrix.Fqq << std::endl;
  std::cout << "kkt_matrix.Fqq_prev" << std::endl;
  std::cout << kkt_matrix.Fqq_prev << std::endl;
}


TEST_F(StateEquationTest, backwardEuler_fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  stateequation::LinearizeBackwardEuler(robot, dtau_, q_prev, v_prev, s, s_next,
                                        kkt_matrix,  kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq.isZero());
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau_, q_prev, v_prev, 
                                                s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((-1*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq.isZero());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(StateEquationTest, backwardEuler_floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  KKTResidual kkt_residual(robot);
  KKTMatrix kkt_matrix(robot);
  stateequation::LinearizeBackwardEuler(robot, dtau_, q_prev, v_prev, s, s_next,
                                        kkt_matrix,  kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  Eigen::MatrixXd dsubtract_dqminus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dqplus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dqminus);
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dqplus);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqplus.transpose()*s_next.lmd+dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq.isApprox(dsubtract_dqminus));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau_, q_prev, v_prev, 
                                                s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq.isApprox(dsubtract_dqminus));
  std::cout << "kkt_matrix.Fqq" << std::endl;
  std::cout << kkt_matrix.Fqq << std::endl;
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), 
                   stateequation::L1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), 
                   stateequation::SquaredNormStateEuqationResidual(kkt_residual));
}


TEST_F(StateEquationTest, normWithStepSizeFixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  const Eigen::VectorXd dq_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv_next = Eigen::VectorXd::Random(robot.dimv());
  KKTResidual kkt_residual(robot);
  const double step_size = 0.3;
  stateequation::ComputeForwardEulerResidual(robot, step_size, dtau_, s, 
                                             s_next.q, s_next.v, dq_next, 
                                             dv_next, kkt_residual);
  const double forward_l1 = stateequation::L1NormStateEuqationResidual(kkt_residual);
  const double forward_squared = stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_residual.setZero();
  kkt_residual.Fq() = s.q + dtau_ * s.v - s_next.q - step_size * dq_next;
  kkt_residual.Fv() = s.v + dtau_ * s.a - s_next.v - step_size * dv_next;
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), forward_l1);
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), forward_squared);
  stateequation::ComputeBackwardEulerResidual(robot, dtau_, q_prev, v_prev,
                                              s, kkt_residual);
  const double backrward_l1_initial = stateequation::L1NormStateEuqationResidual(kkt_residual);
  const double backrward_squared_initial = stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_residual.setZero();
  kkt_residual.Fq() = q_prev - s.q + dtau_ * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau_ * s.a;
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), backrward_l1_initial);
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), backrward_squared_initial);
  const Eigen::VectorXd dq_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv_prev = Eigen::VectorXd::Random(robot.dimv());
  stateequation::ComputeBackwardEulerResidual(robot, step_size, dtau_, 
                                              q_prev, v_prev, dq_prev, dv_prev, 
                                              s, kkt_residual);
  const double backrward_l1 = stateequation::L1NormStateEuqationResidual(kkt_residual);
  const double backrward_squared = stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_residual.setZero();
  kkt_residual.Fq() = q_prev + step_size * dq_prev - s.q + dtau_ * s.v;
  kkt_residual.Fv() = v_prev + step_size * dv_prev - s.v + dtau_ * s.a;
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), backrward_l1);
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), backrward_squared);
}


TEST_F(StateEquationTest, normWithStepSizeFloatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  const Eigen::VectorXd dq_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv_next = Eigen::VectorXd::Random(robot.dimv());
  KKTResidual kkt_residual(robot);
  const double step_size = 0.3;
  stateequation::ComputeForwardEulerResidual(robot, step_size, dtau_, s, 
                                             s_next.q, s_next.v, dq_next, 
                                             dv_next, kkt_residual);
  const double forward_l1 = stateequation::L1NormStateEuqationResidual(kkt_residual);
  const double forward_squared = stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_residual.setZero();
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  kkt_residual.Fq() = qdiff + dtau_ * s.v - step_size * dq_next;
  kkt_residual.Fv() = s.v + dtau_ * s.a - s_next.v - step_size * dv_next;
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), forward_l1);
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), forward_squared);
  stateequation::ComputeBackwardEulerResidual(robot, dtau_, q_prev, v_prev,
                                              s, kkt_residual);
  const double backrward_l1_initial = stateequation::L1NormStateEuqationResidual(kkt_residual);
  const double backrward_squared_initial = stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_residual.setZero();
  qdiff.setZero();
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  kkt_residual.Fq() = qdiff + dtau_ * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau_ * s.a;
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), backrward_l1_initial);
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), backrward_squared_initial);
  const Eigen::VectorXd dq_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv_prev = Eigen::VectorXd::Random(robot.dimv());
  stateequation::ComputeBackwardEulerResidual(robot, step_size, dtau_, 
                                              q_prev, v_prev, dq_prev, dv_prev, 
                                              s, kkt_residual);
  const double backrward_l1 = stateequation::L1NormStateEuqationResidual(kkt_residual);
  const double backrward_squared = stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_residual.setZero();
  qdiff.setZero();
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  kkt_residual.Fq() = qdiff + step_size * dq_prev + dtau_ * s.v;
  kkt_residual.Fv() = v_prev + step_size * dv_prev - s.v + dtau_ * s.a;
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().lpNorm<1>(), backrward_l1);
  EXPECT_DOUBLE_EQ(kkt_residual.Fx().squaredNorm(), backrward_squared);
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
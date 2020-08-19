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


TEST_F(StateEquationTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  Robot robot(fixed_base_urdf_, contact_frames, baum_a, baum_b);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
  s.lmd = Eigen::VectorXd::Random(robot.dimv());
  s.gmm = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_next);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  StateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, lmd_next, 
                                        gmm_next, q_next, kkt_matrix,  kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((lmd_next-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+gmm_next)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  state_equation.linearizeStateEquationTerminal(robot, dtau_, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((-1*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
}


TEST_F(StateEquationTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  Robot robot(floating_base_urdf_, contact_frames, baum_a, baum_b);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0, rnd()%2==0, rnd()%2==0, rnd()%2==0};
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
  s.lmd = Eigen::VectorXd::Random(robot.dimv());
  s.gmm = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd lmd_next = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm_next = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd q_next = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_next);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  StateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, dtau_, q_prev, v_prev, s, lmd_next, 
                                        gmm_next, q_next, kkt_matrix,  kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  Eigen::MatrixXd dsubtract_dqminus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dqplus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dqminus);
  robot.dSubtractdConfigurationPlus(s.q, q_next, dsubtract_dqplus);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqplus.transpose()*lmd_next+dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm+gmm_next)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  state_equation.linearizeStateEquationTerminal(robot, dtau_, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dtau_*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dtau_*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dtau_*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la().isApprox((dtau_*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  EXPECT_TRUE(kkt_matrix.Fqv()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fvv()
              .isApprox(-1*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fva()
              .isApprox(dtau_*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  std::cout << "kkt_matrix.Fqq()" << std::endl;
  std::cout << kkt_matrix.Fqq() << std::endl;
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
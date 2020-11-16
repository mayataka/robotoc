#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_lp_factorizer.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"


namespace idocp {

class StateConstraintRiccatiLPFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 10;
    max_num_impulse = N;
    T = 1;
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  static ContactStatus createRamdomContactStatus(const Robot& robot);
  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse;
  double T, dtau;
};


ContactStatus StateConstraintRiccatiLPFactorizerTest::createRamdomContactStatus(const Robot& robot) {
  ContactStatus contact_status = robot.createContactStatus();
  contact_status.setRandom();
  return contact_status;
}


void StateConstraintRiccatiLPFactorizerTest::test(const Robot& robot) const {
  RiccatiFactorization impulse_riccati_factorization(robot);
  impulse_riccati_factorization.N.setRandom();
  impulse_riccati_factorization.N = (impulse_riccati_factorization.N * impulse_riccati_factorization.N.transpose()).eval();
  impulse_riccati_factorization.Pi.setRandom();
  impulse_riccati_factorization.pi.setRandom();
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  constraint_factorization.setImpulseStatus(impulse_status);
  constraint_factorization.Eq().setRandom();
  constraint_factorization.e().setRandom();
  RiccatiFactorization impulse_riccati_factorization_ref = impulse_riccati_factorization;
  StateConstraintRiccatiFactorization constraint_factorization_ref = constraint_factorization;
  const Eigen::VectorXd dx0 = Eigen::VectorXd::Random(2*robot.dimv());
  StateConstraintRiccatiLPFactorizer factorizer(robot);
  factorizer.factorizeLinearProblem(impulse_riccati_factorization, constraint_factorization, dx0);
  const int dimf = impulse_status.dimp();
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  Eigen::MatrixXd E = Eigen::MatrixXd::Zero(dimf, dimx);
  E.leftCols(dimv) = constraint_factorization_ref.Eq();
  const Eigen::MatrixXd EN = E * impulse_riccati_factorization_ref.N;
  const Eigen::MatrixXd ENEt = EN * E.transpose();
  const Eigen::VectorXd e = constraint_factorization_ref.e() 
                            + E * impulse_riccati_factorization_ref.Pi * dx0
                            + E * impulse_riccati_factorization_ref.pi;
  EXPECT_TRUE(EN.isApprox(constraint_factorization.EN()));
  EXPECT_TRUE(ENEt.isApprox(constraint_factorization.ENEt()));
  EXPECT_TRUE(e.isApprox(constraint_factorization.e()));
  std::cout << dimf << std::endl;
}


TEST_F(StateConstraintRiccatiLPFactorizerTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(StateConstraintRiccatiLPFactorizerTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
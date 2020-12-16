#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"

#include "test_helper.hpp"

namespace idocp {

class StateConstraintRiccatiFactorizationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    T = 1;
    N = 20;
    max_num_impulse = 5;
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  ContactSequence createContactSequence(const Robot& robot) const;
  void testDimension(const Robot& robot) const;
  void testAssignment(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double t, T, dtau;
  int N, max_num_impulse;
};


ContactSequence StateConstraintRiccatiFactorizationTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dtau);
}


void StateConstraintRiccatiFactorizationTest::testDimension(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  StateConstraintRiccatiFactorization factorization(robot, N, max_num_impulse);
  for (int constraint_index=0; constraint_index<max_num_impulse; ++constraint_index) {
    for (int i=0; i<N; ++i) {
      EXPECT_EQ(factorization.T(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T(constraint_index, i).cols(), 0);
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.ENT(constraint_index, i).rows(), 0);
      EXPECT_EQ(factorization.ENT(constraint_index, i).cols(), 0);
    }
    EXPECT_EQ(factorization.Eq(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.Eq(constraint_index).cols(), dimv);
    EXPECT_EQ(factorization.EN(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.EN(constraint_index).cols(), dimx);
    EXPECT_EQ(factorization.ENq(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.ENq(constraint_index).cols(), dimv);
    EXPECT_EQ(factorization.ENEt(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.ENEt(constraint_index).cols(), 0);
    EXPECT_EQ(factorization.e(constraint_index).size(), 0);
    EXPECT_EQ(factorization.dxi(constraint_index).size(), 0);
  }
  EXPECT_EQ(factorization.ENT().rows(), 0);
  EXPECT_EQ(factorization.ENT().cols(), 0);
  EXPECT_EQ(factorization.e().size(), 0);
  EXPECT_EQ(factorization.dxi().size(), 0);
  const auto contact_sequence = createContactSequence(robot);
  const int num_impulse = contact_sequence.numImpulseEvents();
  factorization.setConstraintStatus(contact_sequence);
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    const int dimc = contact_sequence.impulseStatus(constraint_index).dimf();
    for (int i=0; i<N; ++i) {
      EXPECT_EQ(factorization.T(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T(constraint_index, i).cols(), dimc);
    }
    for (int i=0; i<num_impulse; ++i) {
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).cols(), dimc);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).cols(), dimc);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).cols(), dimc);
      const int dimf = contact_sequence.impulseStatus(i).dimf();
      EXPECT_EQ(factorization.ENT(constraint_index, i).rows(), dimc);
      EXPECT_EQ(factorization.ENT(constraint_index, i).cols(), dimf);
    }
    for (int i=num_impulse; i<max_num_impulse; ++i) {
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).cols(), dimc);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).cols(), dimc);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).cols(), dimc);
      EXPECT_EQ(factorization.ENT(constraint_index, i).rows(), dimc);
      EXPECT_EQ(factorization.ENT(constraint_index, i).cols(), 0);
    }
    EXPECT_EQ(factorization.Eq(constraint_index).rows(), dimc);
    EXPECT_EQ(factorization.Eq(constraint_index).cols(), dimv);
    EXPECT_EQ(factorization.EN(constraint_index).rows(), dimc);
    EXPECT_EQ(factorization.EN(constraint_index).cols(), dimx);
    EXPECT_EQ(factorization.ENq(constraint_index).rows(), dimc);
    EXPECT_EQ(factorization.ENq(constraint_index).cols(), dimv);
    EXPECT_EQ(factorization.ENEt(constraint_index).rows(), dimc);
    EXPECT_EQ(factorization.ENEt(constraint_index).cols(), dimc);
    EXPECT_EQ(factorization.e(constraint_index).size(), dimc);
    EXPECT_EQ(factorization.dxi(constraint_index).size(), dimc);
  }
  for (int constraint_index=num_impulse; constraint_index<max_num_impulse; ++constraint_index) {
    for (int i=0; i<N; ++i) {
      EXPECT_EQ(factorization.T(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T(constraint_index, i).cols(), 0);
    }
    for (int i=0; i<num_impulse; ++i) {
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).cols(), 0);
      const int dimf = contact_sequence.impulseStatus(i).dimf();
      EXPECT_EQ(factorization.ENT(constraint_index, i).rows(), 0);
      EXPECT_EQ(factorization.ENT(constraint_index, i).cols(), dimf);
    }
    for (int i=num_impulse; i<max_num_impulse; ++i) {
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_impulse(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_aux(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).rows(), dimx);
      EXPECT_EQ(factorization.T_lift(constraint_index, i).cols(), 0);
      EXPECT_EQ(factorization.ENT(constraint_index, i).rows(), 0);
      EXPECT_EQ(factorization.ENT(constraint_index, i).cols(), 0);
    }
    EXPECT_EQ(factorization.Eq(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.Eq(constraint_index).cols(), dimv);
    EXPECT_EQ(factorization.EN(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.EN(constraint_index).cols(), dimx);
    EXPECT_EQ(factorization.ENq(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.ENq(constraint_index).cols(), dimv);
    EXPECT_EQ(factorization.ENEt(constraint_index).rows(), 0);
    EXPECT_EQ(factorization.ENEt(constraint_index).cols(), 0);
    EXPECT_EQ(factorization.e(constraint_index).size(), 0);
    EXPECT_EQ(factorization.dxi(constraint_index).size(), 0);
  }
  int dimf_total = 0;
  for (int i=0; i<num_impulse; ++i) {
    dimf_total += contact_sequence.impulseStatus(i).dimf();
  }
  EXPECT_EQ(factorization.ENT().rows(), dimf_total);
  EXPECT_EQ(factorization.ENT().cols(), dimf_total);
  EXPECT_EQ(factorization.e().size(), dimf_total);
  EXPECT_EQ(factorization.dxi().size(), dimf_total);
}


void StateConstraintRiccatiFactorizationTest::testAssignment(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const auto contact_sequence = createContactSequence(robot);
  const int num_impulse = contact_sequence.numImpulseEvents();
  StateConstraintRiccatiFactorization factorization(robot, N, max_num_impulse);
  factorization.setConstraintStatus(contact_sequence);
  std::vector<std::vector<Eigen::MatrixXd>> T, T_impulse, T_aux, T_lift;
  std::vector<Eigen::MatrixXd> Eq, EN;
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    const int dimc = contact_sequence.impulseStatus(constraint_index).dimf();
    std::vector<Eigen::MatrixXd> T_vec;
    for (int i=0; i<N; ++i) {
      T_vec.push_back(Eigen::MatrixXd::Random(dimx, dimc));
      factorization.T(constraint_index, i) = T_vec[i];
    }
    T.push_back(T_vec);
    std::vector<Eigen::MatrixXd> T_aux_vec, T_impulse_vec, T_lift_vec;
    for (int i=0; i<num_impulse; ++i) {
      T_impulse_vec.push_back(Eigen::MatrixXd::Random(dimx, dimc));
      T_aux_vec.push_back(Eigen::MatrixXd::Random(dimx, dimc));
      T_lift_vec.push_back(Eigen::MatrixXd::Random(dimx, dimc));
      factorization.T_impulse(constraint_index, i) = T_impulse_vec[i];
      factorization.T_aux(constraint_index, i) = T_aux_vec[i];
      factorization.T_lift(constraint_index, i) = T_lift_vec[i];
    }
    T_impulse.push_back(T_impulse_vec);
    T_aux.push_back(T_aux_vec);
    T_lift.push_back(T_lift_vec);
    Eq.push_back(Eigen::MatrixXd::Random(dimc, dimv));
    EN.push_back(Eigen::MatrixXd::Random(dimc, dimx));
    factorization.Eq(constraint_index) = Eq[constraint_index];
    factorization.EN(constraint_index) = EN[constraint_index];
  }
  int dimf_total = 0;
  for (int i=0; i<num_impulse; ++i) {
    dimf_total += contact_sequence.impulseStatus(i).dimf();
  }
  const Eigen::MatrixXd ENT = Eigen::MatrixXd::Random(dimf_total, dimf_total);
  const Eigen::MatrixXd e = Eigen::VectorXd::Random(dimf_total);
  const Eigen::MatrixXd dxi = Eigen::VectorXd::Random(dimf_total);
  factorization.ENT() = ENT;
  factorization.e() = e;
  factorization.dxi() = dxi;
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    const auto T_vec = T[constraint_index];
    const auto T_impulse_vec = T_impulse[constraint_index];
    const auto T_aux_vec = T_aux[constraint_index];
    const auto T_lift_vec = T_lift[constraint_index];
    for (int i=0; i<N; ++i) {
      EXPECT_TRUE(factorization.T(constraint_index, i).isApprox(T_vec[i]));
    }
    for (int i=0; i<num_impulse; ++i) {
      EXPECT_TRUE(factorization.T_impulse(constraint_index, i).isApprox(T_impulse_vec[i]));
      EXPECT_TRUE(factorization.T_aux(constraint_index, i).isApprox(T_aux_vec[i]));
      EXPECT_TRUE(factorization.T_lift(constraint_index, i).isApprox(T_lift_vec[i]));
      EXPECT_FALSE(factorization.T_impulse(constraint_index, i).isZero());
      EXPECT_FALSE(factorization.T_aux(constraint_index, i).isZero());
      EXPECT_FALSE(factorization.T_lift(constraint_index, i).isZero());
    }
    EXPECT_TRUE(factorization.Eq(constraint_index).isApprox(Eq[constraint_index]));
    EXPECT_TRUE(factorization.EN(constraint_index).isApprox(EN[constraint_index]));
    EXPECT_TRUE(factorization.ENq(constraint_index).isApprox(EN[constraint_index].leftCols(dimv)));
    EXPECT_FALSE(factorization.Eq(constraint_index).isZero());
    EXPECT_FALSE(factorization.EN(constraint_index).isZero());
    EXPECT_FALSE(factorization.ENq(constraint_index).isZero());
  }
  EXPECT_TRUE(factorization.ENT().isApprox(ENT));
  EXPECT_TRUE(factorization.e().isApprox(e));
  EXPECT_TRUE(factorization.dxi().isApprox(dxi));
  EXPECT_FALSE(factorization.ENT().isZero());
  EXPECT_FALSE(factorization.e().isZero());
  EXPECT_FALSE(factorization.dxi().isZero());
  int dimf_stack = 0;
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    const int dimc = contact_sequence.impulseStatus(constraint_index).dimf();
    int dimf_stack_inner = 0;
    EXPECT_TRUE(factorization.ENEt(constraint_index).isApprox(factorization.ENT().block(dimf_stack, dimf_stack, dimc, dimc)));
    for (int i=0; i<num_impulse; ++i) {
      const int dimf = contact_sequence.impulseStatus(i).dimf();
      EXPECT_TRUE(factorization.ENT(constraint_index, i).isApprox(factorization.ENT().block(dimf_stack, dimf_stack_inner, dimc, dimf)));
      dimf_stack_inner += dimf;
    }
    EXPECT_TRUE(factorization.e(constraint_index).isApprox(factorization.e().segment(dimf_stack, dimc)));
    EXPECT_TRUE(factorization.dxi(constraint_index).isApprox(factorization.dxi().segment(dimf_stack, dimc)));
    dimf_stack += dimc;
  }
}


TEST_F(StateConstraintRiccatiFactorizationTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testDimension(robot);
  testAssignment(robot);
}


TEST_F(StateConstraintRiccatiFactorizationTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testDimension(robot);
  testAssignment(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
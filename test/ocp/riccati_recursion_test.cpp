#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/split_riccati_factorizer.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"

#include "test_helper.hpp"

namespace idocp {

class RiccatiRecursionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
    N = 20;
    max_num_impulse = 5;
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N;
  }

  virtual void TearDown() {
  }

  KKTMatrix createKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence) const;
  KKTResidual createKKTResidual(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void test(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


KKTMatrix RiccatiRecursionTest::createKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateKKTMatrix(robot, contact_sequence, N, max_num_impulse);
}


KKTResidual RiccatiRecursionTest::createKKTResidual(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateKKTResidual(robot, contact_sequence, N, max_num_impulse);
}


ContactSequence RiccatiRecursionTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


void RiccatiRecursionTest::test(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  OCPLinearizer linearizer(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  KKTMatrix kkt_matrix(robot, N, max_num_impulse);
  KKTResidual kkt_residual(robot, N, max_num_impulse);
  std::vector<Robot> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  StateConstraintJacobian jac(robot, max_num_impulse);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto s = testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
  linearizer.initConstraints(ocp, robots, contact_sequence, s);
  linearizer.linearizeOCP(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual, jac);
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  auto jac_ref = jac;
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  auto factorization_ref = factorization;
  RiccatiRecursion riccati_recursion(robot, N, max_num_impulse);
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  const auto ocp_discretizer = ocp.discrete();
  riccati_recursion.backwardRiccatiRecursion(ocp_discretizer, kkt_matrix, kkt_residual, jac, factorization);
  factorization_ref[ocp_discretizer.N()].Pqq = kkt_matrix_ref[ocp_discretizer.N()].Qqq();
  factorization_ref[ocp_discretizer.N()].Pvv = kkt_matrix_ref[ocp_discretizer.N()].Qvv();
  factorization_ref[ocp_discretizer.N()].sq = - kkt_residual_ref[ocp_discretizer.N()].lq();
  factorization_ref[ocp_discretizer.N()].sv = - kkt_residual_ref[ocp_discretizer.N()].lv();
  for (int i=N-1; i>=0; --i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_aux = ocp_discretizer.dt_aux(impulse_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      factorizer.aux[impulse_index].backwardRiccatiRecursion(
          factorization_ref[i+1], dt_aux, kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], factorization_ref.aux[impulse_index]);
      factorizer.impulse[impulse_index].backwardRiccatiRecursion(
          factorization_ref.aux[impulse_index], kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], factorization_ref.impulse[impulse_index]);
      factorizer[i].backwardRiccatiRecursion(
          factorization_ref.impulse[impulse_index], dti, kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_lift = ocp_discretizer.dt_lift(lift_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dt);
      if (ocp_discretizer.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i+1);
        factorizer.lift[lift_index].backwardRiccatiRecursion(
            factorization_ref[i+1], ocp_discretizer.dt_lift(lift_index), 
            kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
            jac_ref[impulse_index], factorization_ref.lift[lift_index],
            factorization_ref.constraint[impulse_index]);
      }
      else {
        factorizer.lift[lift_index].backwardRiccatiRecursion(
            factorization_ref[i+1], ocp_discretizer.dt_lift(lift_index), 
            kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
            factorization_ref.lift[lift_index]);
      }
      factorizer[i].backwardRiccatiRecursion(
          factorization_ref.lift[lift_index], dti, kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i]);
    }
    else {
      if (ocp_discretizer.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i+1);
        factorizer[i].backwardRiccatiRecursion(
            factorization_ref[i+1], ocp_discretizer.dt(i), kkt_matrix_ref[i], 
            kkt_residual_ref[i], jac_ref[impulse_index], factorization_ref[i], 
            factorization_ref.constraint[impulse_index]);
      }
      else {
        factorizer[i].backwardRiccatiRecursion(
            factorization_ref[i+1], dt, kkt_matrix_ref[i], 
            kkt_residual_ref[i], factorization_ref[i]);
      }
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.constraint[i].isApprox(factorization_ref.constraint[i]));
  }
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=ocp_discretizer.N(); ++i) {
    d[i].dx().setRandom();
    d[i].du().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.impulse[i].dx().setRandom();
    d.aux[i].dx().setRandom();
    d.aux[i].du().setRandom();
    d.lift[i].dx().setRandom();
    d.lift[i].du().setRandom();
  }
  auto d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(ocp_discretizer, kkt_matrix, kkt_residual, d);
  for (int i=0; i<ocp_discretizer.N(); ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_aux = ocp_discretizer.dt_aux(impulse_index);
      factorizer[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], dti, d_ref[i], d_ref.impulse[impulse_index]);
      factorizer.impulse[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index], 
          d_ref.impulse[impulse_index], d_ref.aux[impulse_index]);
      factorizer.aux[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index], 
          dt_aux, d_ref.aux[impulse_index], d_ref[i+1]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_lift = ocp_discretizer.dt_lift(lift_index);
      factorizer[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], 
          dti, d_ref[i], d_ref.lift[lift_index]);
      factorizer.lift[lift_index].forwardRiccatiRecursion(
          kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
          dt_lift, d_ref.lift[lift_index], d_ref[i+1]);
    }
    else {
      factorizer[i].forwardRiccatiRecursion(kkt_matrix_ref[i], 
                                            kkt_residual_ref[i], dt, 
                                            d_ref[i], d_ref[i+1]);
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
}


TEST_F(RiccatiRecursionTest, fixedBase) {
  test(fixed_base_robot);
}


TEST_F(RiccatiRecursionTest, floating_base) {
  test(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
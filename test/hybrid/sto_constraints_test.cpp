#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/hybrid/sto_constraints.hpp"

#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"


namespace robotoc {

class STOConstraintsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_impulse = 5;
    T = 1.0;
    t = 0.1; 
    min_dt  = std::abs(Eigen::VectorXd::Random(1)[0]);

    auto robot = testhelper::CreateFloatingBaseRobot();
    kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
    kkt_residual = KKTResidual(robot, N, max_num_impulse);
    d = Direction(robot, N, max_num_impulse);

    const double dt = T / N;
    contact_sequence 
        = testhelper::CreateContactSequenceSharedPtr(robot, N, max_num_impulse, 0, 3*dt);

    discretization = HybridOCPDiscretization(T, N, 2*max_num_impulse);
    discretization.discretize(contact_sequence, t);
  }

  virtual void TearDown() {
  }

  int N, max_num_impulse;
  double T, t, min_dt;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  Direction d;
  std::shared_ptr<ContactSequence> contact_sequence;
  HybridOCPDiscretization discretization;
};


TEST_F(STOConstraintsTest, test) {
  const int max_num_switches = 2 * max_num_impulse;
  auto constraints = STOConstraints(max_num_switches, min_dt);
  constraints.setSlack(discretization);
  constraints.evalConstraint(discretization);
  constraints.linearizeConstraints(discretization, kkt_residual);
  constraints.condenseSlackAndDual(discretization, kkt_matrix, kkt_residual);
  constraints.expandSlackAndDual(discretization, d);
  const double primal_step_size = constraints.maxPrimalStepSize();
  const double dual_step_size = constraints.maxDualStepSize();
  EXPECT_TRUE(0 < primal_step_size && primal_step_size <= 1.0);
  EXPECT_TRUE(0 < dual_step_size && dual_step_size <= 1.0);
  constraints.updateSlack(primal_step_size);
  constraints.updateDual(dual_step_size);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/periodic_com_ref.hpp"


namespace robotoc {

class PeriodicCoMRefTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    CoM_ref0 = Eigen::Vector3d::Random();
    v_CoM_ref = Eigen::Vector3d::Random();

    t0 = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_active = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_inactive = std::abs(Eigen::VectorXd::Random(1)[0]);
    period = period_active + period_inactive;
  }

  virtual void TearDown() {
  }

  Eigen::Vector3d CoM_ref0, v_CoM_ref;
  double t0, period_active, period_inactive, period;
};


TEST_F(PeriodicCoMRefTest, first_mode_half_true) {
  auto preiodic_com_ref = std::make_shared<PeriodicCoMRef>(CoM_ref0, v_CoM_ref,
                                                           t0, period_active,
                                                           period_inactive,
                                                           true);
  Eigen::VectorXd com(3), com_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  EXPECT_FALSE(preiodic_com_ref->isActive(t1));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_com_ref->update_CoM_ref(t2, com);
  if (t2 < t0+period_active) {
    com_ref = CoM_ref0 + 0.5 * (t2-t0) * v_CoM_ref;
    EXPECT_TRUE(com_ref.isApprox(com));
    EXPECT_TRUE(preiodic_com_ref->isActive(t2));
  }
  else if (t2 < t0+period) {
    EXPECT_FALSE(preiodic_com_ref->isActive(t2));
  }
  else {
    const int steps = std::floor((t2-t0)/period);
    const double tau = t2 - t0 - steps*period;
    if (tau < period_active) {
      com_ref = CoM_ref0 + period_active*(steps-0.5)*v_CoM_ref + tau*v_CoM_ref;
      EXPECT_TRUE(com_ref.isApprox(com));
      EXPECT_TRUE(preiodic_com_ref->isActive(t2));
    }
    else {
      EXPECT_FALSE(preiodic_com_ref->isActive(t2));
    }
  }
}


TEST_F(PeriodicCoMRefTest, first_mode_half_false) {
  auto preiodic_com_ref = std::make_shared<PeriodicCoMRef>(CoM_ref0, v_CoM_ref,
                                                           t0, period_active,
                                                           period_inactive,
                                                           false);
  Eigen::VectorXd com(3), com_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  EXPECT_FALSE(preiodic_com_ref->isActive(t1));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_com_ref->update_CoM_ref(t2, com);
  const int steps = std::floor((t2-t0)/period);
  const double tau = t2 - t0 - steps*period;
  if (tau < period_active) {
    com_ref = CoM_ref0 + period_active*steps*v_CoM_ref + tau*v_CoM_ref;
    EXPECT_TRUE(preiodic_com_ref->isActive(t2));
    EXPECT_TRUE(com_ref.isApprox(com));
  }
  else {
    EXPECT_FALSE(preiodic_com_ref->isActive(t2));
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/periodic_com_ref.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

class PeriodicCoMRefTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    com_ref0 = Eigen::Vector3d::Random();
    vcom_ref = Eigen::Vector3d::Random();

    t0 = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_active = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_inactive = std::abs(Eigen::VectorXd::Random(1)[0]);
    period = period_active + period_inactive;
  }

  virtual void TearDown() {
  }

  Eigen::Vector3d com_ref0, vcom_ref;
  double t0, period_active, period_inactive, period;
};


TEST_F(PeriodicCoMRefTest, first_mode_half_true) {
  auto preiodic_com_ref = std::make_shared<PeriodicCoMRef>(com_ref0, vcom_ref,
                                                           t0, period_active,
                                                           period_inactive,
                                                           true);
  Eigen::VectorXd com(3), com_ref(3);
  auto grid_info = GridInfo();
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  grid_info.t = t1;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  grid_info.t = t2;
  preiodic_com_ref->updateRef(grid_info, com);
  if (t2 < t0+period_active) {
    com_ref = com_ref0 + 0.5 * (t2-t0) * vcom_ref;
    EXPECT_TRUE(com_ref.isApprox(com));
    EXPECT_TRUE(preiodic_com_ref->isActive(grid_info));
  }
  else if (t2 < t0+period) {
    EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
  }
  else {
    const int steps = std::floor((t2-t0)/period);
    const double tau = t2 - t0 - steps*period;
    if (tau < period_active) {
      com_ref = com_ref0 + period_active*(steps-0.5)*vcom_ref + tau*vcom_ref;
      EXPECT_TRUE(com_ref.isApprox(com));
      EXPECT_TRUE(preiodic_com_ref->isActive(grid_info));
    }
    else {
      EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
    }
  }
}


TEST_F(PeriodicCoMRefTest, first_mode_half_false) {
  auto preiodic_com_ref = std::make_shared<PeriodicCoMRef>(com_ref0, vcom_ref,
                                                           t0, period_active,
                                                           period_inactive,
                                                           false);
  Eigen::VectorXd com(3), com_ref(3);
  auto grid_info = GridInfo();
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  grid_info.t = t1;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  grid_info.t = t2;
  preiodic_com_ref->updateRef(grid_info, com);
  const int steps = std::floor((t2-t0)/period);
  const double tau = t2 - t0 - steps*period;
  if (tau < period_active) {
    com_ref = com_ref0 + period_active*steps*vcom_ref + tau*vcom_ref;
    EXPECT_TRUE(preiodic_com_ref->isActive(grid_info));
    EXPECT_TRUE(com_ref.isApprox(com));
  }
  else {
    EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/periodic_swing_foot_ref.hpp"
#include "robotoc/ocp/grid_info.hpp"


namespace robotoc {

class PeriodicSwingFootRefTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    p0 = Eigen::Vector3d::Random();
    step_length = Eigen::Vector3d::Random();

    t0 = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_swing = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_stance = std::abs(Eigen::VectorXd::Random(1)[0]);
    period = period_swing + period_stance;
  }

  virtual void TearDown() {
  }

  Eigen::Vector3d p0, step_length;
  double step_height, t0, period_swing, period_stance, period;
};


TEST_F(PeriodicSwingFootRefTest, first_mode_half_true) {
  auto preiodic_foot_ref = std::make_shared<PeriodicSwingFootRef>(p0, step_length,
                                                                  step_height, t0, 
                                                                  period_swing,
                                                                  period_stance, true);
  Eigen::VectorXd p(3), p_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  GridInfo grid_info;
  grid_info.t = t1;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  grid_info.t = t2;
  preiodic_foot_ref->updateRef(grid_info, p);
  if (t2 < t0+period_swing) {
    p_ref = p0;
    p_ref += 0.5 * ((t2-t0)/period_swing) * step_length;
    const double tau = t2-t0;
    const double rate = tau / period_swing;
    if (rate < 0.5) {
      p_ref(2) += 2 * step_height * rate;
    }
    else {
      p_ref(2) += 2 * step_height * (1-rate);
    }
    EXPECT_TRUE(p.isApprox(p_ref));
    EXPECT_TRUE(preiodic_foot_ref->isActive(grid_info));
  }
  else if (t2 < period) {
    EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
  }
  else {
    const int steps = std::floor((t2-t0)/period);
    const double tau = t2 - t0 - steps*period;
    if (tau < period_swing) {
      p_ref = p0;
      p_ref += (steps-0.5)*step_length; 
      p_ref += (tau/period_swing)*step_length; 
      const double rate = tau / period_swing;
      if (rate < 0.5) {
        p_ref(2) += 2 * step_height * rate;
      }
      else {
        p_ref(2) += 2 * step_height * (1-rate);
      }
      EXPECT_TRUE(p.isApprox(p_ref));
      EXPECT_TRUE(preiodic_foot_ref->isActive(grid_info));
    }
    else {
      EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
    }
  }
}


TEST_F(PeriodicSwingFootRefTest, first_mode_half_false) {
  auto preiodic_foot_ref = std::make_shared<PeriodicSwingFootRef>(p0, step_length,
                                                                  step_height, t0, 
                                                                  period_swing,
                                                                  period_stance, false);
  Eigen::VectorXd p(3), p_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  GridInfo grid_info;
  grid_info.t = t1;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  grid_info.t = t2;
  preiodic_foot_ref->updateRef(grid_info, p);
  const int steps = std::floor((t2-t0)/period);
  const double tau = t2 - t0 - steps*period;
  if (tau < period_swing) {
    p_ref = p0;
    p_ref += steps*step_length; 
    p_ref += (tau/period_swing)*step_length; 
    const double rate = tau / period_swing;
    if (rate < 0.5) {
      p_ref(2) += 2 * step_height * rate;
    }
    else {
      p_ref(2) += 2 * step_height * (1-rate);
    }
    EXPECT_TRUE(p.isApprox(p_ref));
    EXPECT_TRUE(preiodic_foot_ref->isActive(grid_info));
  }
  else {
    EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
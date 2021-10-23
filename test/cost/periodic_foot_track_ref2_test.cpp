#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/periodic_foot_track_ref2.hpp"


namespace robotoc {

class PeriodicFootTrackRefTest2 : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    p0 = Eigen::Vector3d::Random();

    t0 = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_swing = std::abs(Eigen::VectorXd::Random(1)[0]);
    period_stance = std::abs(Eigen::VectorXd::Random(1)[0]);
    period = period_swing + period_stance;
  }

  virtual void TearDown() {
  }

  Eigen::Vector3d p0;
  double step_length, step_height, t0, period_swing, period_stance, period;
};


TEST_F(PeriodicFootTrackRefTest2, first_mode_half_true) {
  auto preiodic_foot_ref = std::make_shared<PeriodicFootTrackRef2>(p0, step_length,
                                                                   step_height, t0, 
                                                                   period_swing,
                                                                   period_stance, true);
  Eigen::VectorXd p(3), p_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_foot_ref->update_q_3d_ref(t1, p);
  EXPECT_TRUE(p.isApprox(p0));
  EXPECT_TRUE(preiodic_foot_ref->isActive(t1));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_foot_ref->update_q_3d_ref(t2, p);
  EXPECT_TRUE(preiodic_foot_ref->isActive(t2));
  if (t2 < t0+period_swing) {
    p_ref = p0;
    p_ref(0) += 0.5 * ((t2-t0)/period_swing) * step_length;
    const double tau = t2-t0;
    const double rate = tau / period_swing;
    if (rate < 0.5) {
      p_ref(2) += 2 * step_height * rate;
    }
    else {
      p_ref(2) += 2 * step_height * (1-rate);
    }
  }
  else if (t2 < period) {
    p_ref = p0;
    p_ref(0) += 0.5 * step_length;
  }
  else {
    const int steps = std::floor((t2-t0)/period);
    const double tau = t2 - t0 - steps*period;
    if (tau < period_swing) {
      p_ref = p0;
      p_ref(0) += (steps-0.5)*step_length; 
      p_ref(0) += (tau/period_swing)*step_length; 
      const double rate = tau / period_swing;
      if (rate < 0.5) {
        p_ref(2) += 2 * step_height * rate;
      }
      else {
        p_ref(2) += 2 * step_height * (1-rate);
      }
    }
    else {
      p_ref = p0;
      p_ref(0) += (steps+0.5)*step_length; 
    }
  }
  EXPECT_TRUE(p.isApprox(p_ref));
}


TEST_F(PeriodicFootTrackRefTest2, first_mode_half_false) {
  auto preiodic_foot_ref = std::make_shared<PeriodicFootTrackRef2>(p0, step_length,
                                                                   step_height, t0, 
                                                                   period_swing,
                                                                   period_stance, false);
  Eigen::VectorXd p(3), p_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_foot_ref->update_q_3d_ref(t1, p);
  EXPECT_TRUE(preiodic_foot_ref->isActive(t1));
  EXPECT_TRUE(p.isApprox(p0));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_foot_ref->update_q_3d_ref(t2, p);
  EXPECT_TRUE(preiodic_foot_ref->isActive(t2));
  const int steps = std::floor((t2-t0)/period);
  const double tau = t2 - t0 - steps*period;
  if (tau < period_swing) {
    p_ref = p0;
    p_ref(0) += steps*step_length; 
    p_ref(0) += (tau/period_swing)*step_length; 
    const double rate = tau / period_swing;
    if (rate < 0.5) {
      p_ref(2) += 2 * step_height * rate;
    }
    else {
      p_ref(2) += 2 * step_height * (1-rate);
    }
  }
  else {
    p_ref = p0;
    p_ref(0) += (steps+1.0)*step_length; 
  }
  EXPECT_TRUE(p.isApprox(p_ref));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
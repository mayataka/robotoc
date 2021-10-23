#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/periodic_com_ref2.hpp"


namespace robotoc {

class PeriodicCoMRefTest2 : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
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


TEST_F(PeriodicCoMRefTest2, first_mode_half_true) {
  auto preiodic_com_ref = std::make_shared<PeriodicCoMRef2>(CoM_ref0, v_CoM_ref,
                                                            t0, period_active,
                                                            period_inactive,
                                                            true);
  Eigen::VectorXd com(3), com_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_com_ref->update_CoM_ref(t1, com);
  EXPECT_TRUE(com.isApprox(CoM_ref0));
  EXPECT_TRUE(preiodic_com_ref->isActive(t1));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_com_ref->update_CoM_ref(t2, com);
  EXPECT_TRUE(preiodic_com_ref->isActive(t2));
  if (t2 < t0+period_active) {
    com_ref = CoM_ref0 + 0.5*(t2-t0) * v_CoM_ref;
  }
  else if (t2 < t0+period) {
    com_ref = CoM_ref0 + 0.5*period_active * v_CoM_ref;
  }
  else {
    const int steps = std::floor((t2-t0)/period);
    const double tau = t2 - t0 - steps*period;
    if (tau < period_active) {
      com_ref = CoM_ref0 + period_active*(steps-0.5)*v_CoM_ref + tau*v_CoM_ref;
    }
    else {
      com_ref = CoM_ref0 + period_active*(steps+0.5)*v_CoM_ref;
    }
  }
  EXPECT_TRUE(com.isApprox(com_ref));
}


TEST_F(PeriodicCoMRefTest2, first_mode_half_false) {
  auto preiodic_com_ref = std::make_shared<PeriodicCoMRef2>(CoM_ref0, v_CoM_ref,
                                                            t0, period_active,
                                                            period_inactive,
                                                            false);
  Eigen::VectorXd com(3), com_ref(3);
  const double t1 = t0 - std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_com_ref->update_CoM_ref(t1, com);
  EXPECT_TRUE(preiodic_com_ref->isActive(t1));
  EXPECT_TRUE(com.isApprox(CoM_ref0));
  const double t2 = t0 + std::abs(Eigen::VectorXd::Random(1)[0]);
  preiodic_com_ref->update_CoM_ref(t2, com);
  EXPECT_TRUE(preiodic_com_ref->isActive(t2));
  const int steps = std::floor((t2-t0)/period);
  const double tau = t2 - t0 - steps*period;
  if (tau < period_active) {
    com_ref = CoM_ref0 + period_active*steps*v_CoM_ref + tau*v_CoM_ref;
  }
  else {
    com_ref = CoM_ref0 + period_active*(steps+1.0)*v_CoM_ref;
  }
  EXPECT_TRUE(com.isApprox(com_ref));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
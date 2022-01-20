#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/periodic_foot_track_ref2.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

class PeriodicFootTrackRefTest2 : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    p0 = Eigen::Vector3d::Random();
    step_length = std::abs(Eigen::VectorXd::Random(1)[0]);
    step_height = std::abs(Eigen::VectorXd::Random(1)[0]);

    start_phase = 2;
    end_phase = 21;
    active_phases = 3;
    inactive_phases = 4;

    grid_info = GridInfo::Random();
    grid_info.N_phase = 11;
    grid_info.grid_count_in_phase = 4;
  }

  virtual void TearDown() {
  }

  Eigen::Vector3d p0;
  double step_length, step_height;
  int start_phase, end_phase, active_phases, inactive_phases;
  GridInfo grid_info;
};


TEST_F(PeriodicFootTrackRefTest2, first_mode_half_true) {
  auto preiodic_foot_ref = std::make_shared<PeriodicFootTrackRef2>(p0, step_length,
                                                                   step_height, 
                                                                   start_phase, 
                                                                   end_phase, 
                                                                   active_phases,
                                                                   inactive_phases, true);
  Eigen::VectorXd p(3), p_ref(3);
  grid_info.contact_phase = start_phase - 1;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
  grid_info.contact_phase = start_phase + 1;
  preiodic_foot_ref->update_x3d_ref(grid_info, p);
  const double rate1 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  p_ref = p0;
  p_ref(0) += 0.5 * rate1 * step_length;
  if (rate1 < 0.5) {
    p_ref(2) += 2 * step_height * rate1;
  }
  else {
    p_ref(2) += 2 * step_height * (1-rate1);
  }
  EXPECT_TRUE(p.isApprox(p_ref));
  EXPECT_TRUE(preiodic_foot_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases + inactive_phases;
  preiodic_foot_ref->update_x3d_ref(grid_info, p);
  const int steps = std::floor((grid_info.contact_phase-start_phase)/(active_phases+inactive_phases));
  const double rate2 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  p_ref = p0;
  p_ref(0) += ((steps-0.5) + rate2) * step_length;
  if (rate2 < 0.5) {
    p_ref(2) += 2 * step_height * rate2;
  }
  else {
    p_ref(2) += 2 * step_height * (1-rate2);
  }
  EXPECT_TRUE(p.isApprox(p_ref));
  EXPECT_TRUE(preiodic_foot_ref->isActive(grid_info));

  grid_info.contact_phase = end_phase;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
}


TEST_F(PeriodicFootTrackRefTest2, first_mode_half_false) {
  auto preiodic_foot_ref = std::make_shared<PeriodicFootTrackRef2>(p0, step_length,
                                                                   step_height, 
                                                                   start_phase, 
                                                                   end_phase, 
                                                                   active_phases,
                                                                   inactive_phases, false);
  Eigen::VectorXd p(3), p_ref(3);
  grid_info.contact_phase = start_phase - 1;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
  grid_info.contact_phase = start_phase + 1;
  preiodic_foot_ref->update_x3d_ref(grid_info, p);
  const double rate1 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  p_ref = p0;
  p_ref(0) += rate1 * step_length;
  if (rate1 < 0.5) {
    p_ref(2) += 2 * step_height * rate1;
  }
  else {
    p_ref(2) += 2 * step_height * (1-rate1);
  }
  EXPECT_TRUE(p.isApprox(p_ref));
  EXPECT_TRUE(preiodic_foot_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases + inactive_phases;
  preiodic_foot_ref->update_x3d_ref(grid_info, p);
  const int steps = std::floor((grid_info.contact_phase-start_phase)/(active_phases+inactive_phases));
  const double rate2 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  p_ref = p0;
  p_ref(0) += (steps + rate2) * step_length;
  if (rate2 < 0.5) {
    p_ref(2) += 2 * step_height * rate2;
  }
  else {
    p_ref(2) += 2 * step_height * (1-rate2);
  }
  EXPECT_TRUE(p.isApprox(p_ref));
  EXPECT_TRUE(preiodic_foot_ref->isActive(grid_info));

  grid_info.contact_phase = end_phase;
  EXPECT_FALSE(preiodic_foot_ref->isActive(grid_info));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
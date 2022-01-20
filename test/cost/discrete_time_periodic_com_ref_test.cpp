#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/discrete_time_periodic_com_ref.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

class DiscreteTimePeriodicCoMRefTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    com_ref0 = Eigen::Vector3d::Random();
    com_step_ref = Eigen::Vector3d::Random();

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

  Eigen::Vector3d com_ref0, com_step_ref;
  int start_phase, end_phase, active_phases, inactive_phases;
  GridInfo grid_info;
};


TEST_F(DiscreteTimePeriodicCoMRefTest, first_mode_half_true) {
  auto preiodic_com_ref = std::make_shared<DiscreteTimePeriodicCoMRef>(com_ref0, 
                                                                       com_step_ref, 
                                                                       start_phase, 
                                                                       end_phase, 
                                                                       active_phases, 
                                                                       inactive_phases, true);
  Eigen::VectorXd com(3), com_ref(3);
  grid_info.contact_phase = start_phase - 1;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
  grid_info.contact_phase = start_phase + 1;
  preiodic_com_ref->update_com_ref(grid_info, com);
  const double rate1 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  com_ref = com_ref0;
  com_ref += 0.5 * rate1 * com_step_ref;
  EXPECT_TRUE(com.isApprox(com_ref));
  EXPECT_TRUE(preiodic_com_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases + inactive_phases;
  preiodic_com_ref->update_com_ref(grid_info, com);
  const int steps = std::floor((grid_info.contact_phase-start_phase)/(active_phases+inactive_phases));
  const double rate2 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  com_ref = com_ref0;
  com_ref += ((steps-0.5) + rate2) * com_step_ref;
  EXPECT_TRUE(com.isApprox(com_ref));
  EXPECT_TRUE(preiodic_com_ref->isActive(grid_info));

  grid_info.contact_phase = end_phase;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
}


TEST_F(DiscreteTimePeriodicCoMRefTest, first_mode_half_false) {
  auto preiodic_com_ref = std::make_shared<DiscreteTimePeriodicCoMRef>(com_ref0, 
                                                                       com_step_ref, 
                                                                       start_phase, 
                                                                       end_phase, 
                                                                       active_phases, 
                                                                       inactive_phases, false);
  Eigen::VectorXd com(3), com_ref(3);
  grid_info.contact_phase = start_phase - 1;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
  grid_info.contact_phase = start_phase + 1;
  preiodic_com_ref->update_com_ref(grid_info, com);
  const double rate1 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  com_ref = com_ref0;
  com_ref += rate1 * com_step_ref;
  EXPECT_TRUE(com.isApprox(com_ref));
  EXPECT_TRUE(preiodic_com_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));

  grid_info.contact_phase = start_phase + 1 + active_phases + inactive_phases;
  preiodic_com_ref->update_com_ref(grid_info, com);
  const int steps = std::floor((grid_info.contact_phase-start_phase)/(active_phases+inactive_phases));
  const double rate2 = static_cast<double>(grid_info.grid_count_in_phase) 
                        / static_cast<double>(grid_info.N_phase);
  com_ref = com_ref0;
  com_ref += (steps + rate2) * com_step_ref;
  EXPECT_TRUE(com.isApprox(com_ref));
  EXPECT_TRUE(preiodic_com_ref->isActive(grid_info));

  grid_info.contact_phase = end_phase;
  EXPECT_FALSE(preiodic_com_ref->isActive(grid_info));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
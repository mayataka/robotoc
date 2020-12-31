#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_temporary_solution.hpp"


namespace idocp {

class SplitTemporarySolutionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ContactStatus& contact_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void SplitTemporarySolutionTest::test(
    const Robot& robot, const ContactStatus& contact_status) { 
  std::random_device rnd;
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  const SplitSolution s_next = SplitSolution::Random(robot, contact_status);
  const SplitDirection d = SplitDirection::Random(robot, contact_status);
  const SplitDirection d_next = SplitDirection::Random(robot, contact_status);
  SplitTemporarySolution s_tmp(robot);
  const double step_size = 0.3;
  s_tmp.setTemporarySolution(robot, contact_status, step_size, s, d, s_next, d_next);
  SplitSolution s_tmp_ref(robot);
  s_tmp_ref.setContactStatus(contact_status);
  s_tmp_ref.q = s.q;
  robot.integrateConfiguration(d.dq(), step_size, s_tmp_ref.q);
  s_tmp_ref.v = s.v + step_size * d.dv();
  s_tmp_ref.a = s.a + step_size * d.da();
  s_tmp_ref.f_stack() = s.f_stack() + step_size * d.df();
  s_tmp_ref.u = s.u + step_size * d.du();
  if (robot.hasFloatingBase()) {
    s_tmp_ref.u_passive = s.u_passive + step_size * d.du_passive;
  }
  Eigen::VectorXd q_next_ref = s_next.q;
  robot.integrateConfiguration(d_next.dq(), step_size, q_next_ref);
  Eigen::VectorXd v_next_ref = s_next.v + step_size * d_next.dv();
  EXPECT_TRUE(s_tmp.splitSolution().isApprox(s_tmp_ref));
  EXPECT_TRUE(s_tmp.q_next().isApprox(q_next_ref));
  EXPECT_TRUE(s_tmp.v_next().isApprox(v_next_ref));
}


TEST_F(SplitTemporarySolutionTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  test(robot, contact_status);
  contact_status.activateContact(0);
  test(robot, contact_status);
}


TEST_F(SplitTemporarySolutionTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false, false, false, false};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  test(robot, contact_status);
  std::random_device rnd;
  is_contact_active.clear();
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
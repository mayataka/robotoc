#include <string>

#include <gtest/gtest.h>

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class KKTCompositionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(KKTCompositionTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  KKTComposition composition(robot);
  composition.setContactStatus(robot);
  EXPECT_EQ(composition.Fq_begin(), 0);
  EXPECT_EQ(composition.Fq_size(), robot.dimv());
  EXPECT_EQ(composition.Fv_begin(), composition.Fq_begin()+composition.Fq_size());
  EXPECT_EQ(composition.Fv_size(), robot.dimv());
  EXPECT_EQ(composition.C_begin(), composition.Fv_begin()+composition.Fv_size());
  EXPECT_EQ(composition.C_size(), robot.dim_passive()+robot.dimf());
  EXPECT_EQ(composition.Qa_begin(), composition.C_begin()+composition.C_size());
  EXPECT_EQ(composition.Qa_size(), robot.dimv());
  EXPECT_EQ(composition.Qf_begin(), composition.Qa_begin()+composition.Qa_size());
  EXPECT_EQ(composition.Qf_size(), robot.dimf());
  EXPECT_EQ(composition.Qq_begin(), composition.Qf_begin()+composition.Qf_size());
  EXPECT_EQ(composition.Qq_size(), robot.dimv());
  EXPECT_EQ(composition.Qv_begin(), composition.Qq_begin()+composition.Qq_size());
  EXPECT_EQ(composition.Qv_size(), robot.dimv());
  EXPECT_EQ(composition.dimKKT(), composition.Qv_begin()+composition.Qv_size());
  EXPECT_EQ(composition.max_dimKKT(), 5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive());
  EXPECT_EQ(composition.Qx_begin(), composition.Qq_begin());
  EXPECT_EQ(composition.Qx_size(), 2*robot.dimv());
}


TEST_F(KKTCompositionTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  KKTComposition composition(robot);
  composition.setContactStatus(robot);
  EXPECT_EQ(composition.Fq_begin(), 0);
  EXPECT_EQ(composition.Fq_size(), robot.dimv());
  EXPECT_EQ(composition.Fv_begin(), composition.Fq_begin()+composition.Fq_size());
  EXPECT_EQ(composition.Fv_size(), robot.dimv());
  EXPECT_EQ(composition.C_begin(), composition.Fv_begin()+composition.Fv_size());
  EXPECT_EQ(composition.C_size(), robot.dim_passive()+robot.dimf());
  EXPECT_EQ(composition.Qa_begin(), composition.C_begin()+composition.C_size());
  EXPECT_EQ(composition.Qa_size(), robot.dimv());
  EXPECT_EQ(composition.Qf_begin(), composition.Qa_begin()+composition.Qa_size());
  EXPECT_EQ(composition.Qf_size(), robot.dimf());
  EXPECT_EQ(composition.Qq_begin(), composition.Qf_begin()+composition.Qf_size());
  EXPECT_EQ(composition.Qq_size(), robot.dimv());
  EXPECT_EQ(composition.Qv_begin(), composition.Qq_begin()+composition.Qq_size());
  EXPECT_EQ(composition.Qv_size(), robot.dimv());
  EXPECT_EQ(composition.dimKKT(), composition.Qv_begin()+composition.Qv_size());
  EXPECT_EQ(composition.max_dimKKT(), 5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive());
  EXPECT_EQ(composition.Qx_begin(), composition.Qq_begin());
  EXPECT_EQ(composition.Qx_size(), 2*robot.dimv());
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
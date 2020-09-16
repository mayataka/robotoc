#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class FloatingBaseContactCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    contact_frames_ = {14, 24, 34, 44};
    robot_ = Robot(urdf_, contact_frames_, 0, 0);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    data_ = CostFunctionData(robot_);
    s = SplitSolution(robot_);
    kkt_res = KKTResidual(robot_);
    kkt_mat = KKTMatrix(robot_);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  std::vector<int> contact_frames_;
  std::string urdf_;
  Robot robot_;
  CostFunctionData data_;
  SplitSolution s;
  KKTResidual kkt_res;
  KKTMatrix kkt_mat;
};


TEST_F(FloatingBaseContactCostTest, setWeights) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<contact_frames_.size(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
  }
  ContactCost cost(robot_);
  EXPECT_FALSE(cost.useKinematics());
  cost.set_f_weight(f_weight);
  cost.set_f_ref(f_ref);
  const Eigen::VectorXd q = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(dimv);
  std::vector<Eigen::Vector3d> f;
  for (int i=0; i<contact_frames_.size(); ++i) {
    f.push_back(Eigen::Vector3d::Random());
  }
  const Eigen::VectorXd u = Eigen::VectorXd::Random(dimv);
  ASSERT_EQ(robot_.dimf(), 0);
  kkt_res.setContactStatus(robot_);
  kkt_mat.setContactStatus(robot_);
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, s), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(robot_.max_dimf());
  cost.lf(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lf().isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(robot_.max_dimf(), robot_.max_dimf());
  cost.lff(robot_, data_, t_, dtau_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  std::random_device rnd;
  std::vector<bool> active_contacts;
  for (int i=0; i<contact_frames_.size(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setContactStatus(active_contacts);
  kkt_res.setContactStatus(robot_);
  kkt_mat.setContactStatus(robot_);
  double l_ref = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (robot_.is_contact_active(i)) {
      l_ref += (f_weight[i].array() * (s.f[i].array()-f_ref[i].array()) 
                                    * (s.f[i].array()-f_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, s), 0.5*dtau_*l_ref);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(robot_.dimf());
  int dimf_stack = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (robot_.is_contact_active(i)) {
      lf_ref.segment<3>(dimf_stack).array()
          = dtau_ * f_weight[i].array() * (s.f[i].array()-f_ref[i].array());
      dimf_stack += 3;
    }
  }
  cost.lf(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lf().isApprox(lf_ref));
  cost.lf(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lf().isApprox(2*lf_ref));
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(robot_.dimf(), robot_.dimf());
  dimf_stack = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (robot_.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        lff_ref.diagonal().segment<3>(dimf_stack) = dtau_ * f_weight[i];
      }
      dimf_stack += 3;
    }
  }
  cost.lff(robot_, data_, t_, dtau_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qff().isApprox(lff_ref));
  cost.lff(robot_, data_, t_, dtau_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qff().isApprox(2*lff_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
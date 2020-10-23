#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class FloatingBaseContactForceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    contact_frames_ = {14, 24, 34, 44};
    robot_ = Robot(urdf_, contact_frames_);
    contact_status_ = ContactStatus(contact_frames_.size());
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
  ContactStatus contact_status_;
  CostFunctionData data_;
  SplitSolution s;
  KKTResidual kkt_res;
  KKTMatrix kkt_mat;
};


TEST_F(FloatingBaseContactForceCostTest, setWeights) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<contact_frames_.size(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
  }
  ContactForceCost cost(robot_);
  EXPECT_FALSE(cost.useKinematics());
  cost.set_f_weight(f_weight);
  cost.set_f_ref(f_ref);
  for (int i=0; i<contact_status_.max_point_contacts(); ++i) {
    contact_status_.deactivateContact(i);
  }
  s = SplitSolution::Random(robot_, contact_status_);
  kkt_res.setContactStatus(contact_status_);
  kkt_mat.setContactStatus(contact_status_);
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, s), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(robot_.max_dimf());
  cost.lf(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lf().isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(robot_.max_dimf(), robot_.max_dimf());
  cost.lff(robot_, data_, t_, dtau_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  std::random_device rnd;
  for (int i=0; i<contact_status_.max_point_contacts(); ++i) {
    if (rnd()%2 == 0) 
      contact_status_.activateContact(i);
  }
  s = SplitSolution::Random(robot_, contact_status_);
  kkt_res.setContactStatus(contact_status_);
  kkt_mat.setContactStatus(contact_status_);
  double l_ref = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (contact_status_.isContactActive(i)) {
      l_ref += (f_weight[i].array() * (s.f[i].array()-f_ref[i].array()) 
                                    * (s.f[i].array()-f_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, s), 0.5*dtau_*l_ref);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(contact_status_.dimf());
  int dimf_stack = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (contact_status_.isContactActive(i)) {
      lf_ref.segment<3>(dimf_stack).array()
          = dtau_ * f_weight[i].array() * (s.f[i].array()-f_ref[i].array());
      dimf_stack += 3;
    }
  }
  cost.lf(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lf().isApprox(lf_ref));
  cost.lf(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lf().isApprox(2*lf_ref));
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(contact_status_.dimf(), contact_status_.dimf());
  dimf_stack = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (contact_status_.isContactActive(i)) {
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
#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/contact_cost.hpp"


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
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::vector<int> contact_frames_;
  std::string urdf_;
  Robot robot_;
};


TEST_F(FloatingBaseContactCostTest, zeroRefernceConstructor) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Zero(dimf);
  ContactCost cost(robot_, f_weight);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setContactStatus(active_contacts);
  double l_ref = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        l_ref += 0.5 * dtau_ * (f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j)) * (f(3*i+j)-f_ref(3*i+j)));
      }
    }
  }
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  int dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lf_ref(dimf_tmp+j) = dtau_ * f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j));
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lf_ref.isApprox(lf));
  cost.lff(robot_, dtau_, lff);
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(dimf, dimf);
  dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lff_ref(dimf_tmp+j, dimf_tmp+j) = dtau_ * f_weight(3*i+j);
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lff_ref.isApprox(lff));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lf_ref.isApprox(lf));
}


TEST_F(FloatingBaseContactCostTest, withRefernceConstructor) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(dimf);
  ContactCost cost(robot_, f_ref, f_weight);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setContactStatus(active_contacts);
  double l_ref = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        l_ref += 0.5 * dtau_ * (f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j)) * (f(3*i+j)-f_ref(3*i+j)));
      }
    }
  }
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  int dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lf_ref(dimf_tmp+j) = dtau_ * f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j));
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lf_ref.isApprox(lf));
  cost.lff(robot_, dtau_, lff);
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(dimf, dimf);
  dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lff_ref(dimf_tmp+j, dimf_tmp+j) = dtau_ * f_weight(3*i+j);
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lff_ref.isApprox(lff));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lf_ref.isApprox(lf));
}


TEST_F(FloatingBaseContactCostTest, setReference) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(dimf);
  ContactCost cost(robot_, Eigen::VectorXd::Zero(dimf), f_weight);
  cost.set_f_ref(f_ref);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setContactStatus(active_contacts);
  double l_ref = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        l_ref += 0.5 * dtau_ * (f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j)) * (f(3*i+j)-f_ref(3*i+j)));
      }
    }
  }
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  int dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lf_ref(dimf_tmp+j) = dtau_ * f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j));
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lf_ref.isApprox(lf));
  cost.lff(robot_, dtau_, lff);
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(dimf, dimf);
  dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lff_ref(dimf_tmp+j, dimf_tmp+j) = dtau_ * f_weight(3*i+j);
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lff_ref.isApprox(lff));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lf_ref.isApprox(lf));
}


TEST_F(FloatingBaseContactCostTest, setWeights) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(dimf);
  ContactCost cost(robot_, Eigen::VectorXd::Zero(dimf), Eigen::VectorXd::Zero(dimf));
  cost.set_f_weight(f_weight);
  cost.set_f_ref(f_ref);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setContactStatus(active_contacts);
  double l_ref = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        l_ref += 0.5 * dtau_ * (f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j)) * (f(3*i+j)-f_ref(3*i+j)));
      }
    }
  }
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  int dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lf_ref(dimf_tmp+j) = dtau_ * f_weight(3*i+j) * (f(3*i+j)-f_ref(3*i+j));
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lf_ref.isApprox(lf));
  cost.lff(robot_, dtau_, lff);
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(dimf, dimf);
  dimf_tmp = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (active_contacts[i]) {
      for (int j=0; j<3; ++j) {
        lff_ref(dimf_tmp+j, dimf_tmp+j) = dtau_ * f_weight(3*i+j);
      }
      dimf_tmp += 3;
    }
  }
  EXPECT_TRUE(lff_ref.isApprox(lff));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lf_ref.isApprox(lf));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
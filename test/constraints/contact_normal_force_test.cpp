#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/constraints/pdipm_func.hpp"
#include "idocp/constraints/contact_normal_force.hpp"

namespace idocp {

class ContactNormalForceTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    const std::string urdf = "../urdf/anymal/anymal.urdf";
    contact_frames = {14, 24, 34, 44};
    robot = Robot(urdf, contact_frames);
    for (int i=0; i<contact_frames.size(); ++i) {
      is_contact_active.push_back(rnd()%2==0);
    }
    for (int i=0; i<contact_frames.size(); ++i) {
      if (is_contact_active[i]) {
        active_contact_indices.push_back(i);
      }
    }
    contact_status = ContactStatus(robot.max_point_contacts());
    contact_status.setContactStatus(is_contact_active);
    barrier = 1.0e-04;
    fraction_to_boundary_rate = 0.995;
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    mu = 0.8;
    data = ConstraintComponentData(contact_frames.size());
    data_ref = ConstraintComponentData(contact_frames.size());
    s = SplitSolution::Random(robot, contact_status);
    d = SplitDirection::Random(robot, contact_status);
    kkt_residual = KKTResidual(robot);
    kkt_matrix = KKTMatrix(robot);
    kkt_residual.setContactStatus(contact_status);
    kkt_matrix.setContactStatus(contact_status);
  }

  virtual void TearDown() {
  }

  double barrier, fraction_to_boundary_rate, dtau, mu;
  Eigen::VectorXd slack, dual, dslack, ddual;
  Robot robot;
  ContactStatus contact_status;
  std::vector<int> contact_frames, active_contact_indices;
  std::vector<bool> is_contact_active;
  ConstraintComponentData data, data_ref;
  SplitSolution s;
  SplitDirection d;
  KKTResidual kkt_residual;
  KKTMatrix kkt_matrix;
};


TEST_F(ContactNormalForceTest, isFeasible) {
  ContactNormalForce contact_normal_force(robot); 
  is_contact_active[0] = true;
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(2) = -1;
  }
  EXPECT_FALSE(contact_normal_force.isFeasible(robot, data, s));
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(2) = 1;
  }
  EXPECT_TRUE(contact_normal_force.isFeasible(robot, data, s));
}


TEST_F(ContactNormalForceTest, setSlackAndDual) {
  ContactNormalForce contact_normal_force(robot); 
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(2) = 1;
  }
  EXPECT_TRUE(contact_normal_force.isFeasible(robot, data, s));
  contact_normal_force.setSlackAndDual(robot, data, dtau, s);
  for (int i=0; i<contact_frames.size(); ++i) {
    EXPECT_DOUBLE_EQ(data.slack.coeff(i), dtau*s.f[i].coeffRef(2));
    EXPECT_TRUE(data.dual.coeff(i) > 0);
  }
}


TEST_F(ContactNormalForceTest, augmentDualResidual) {
  ContactNormalForce contact_normal_force(robot); 
  contact_normal_force.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.lf().isZero());
  data.dual.fill(0.1);
  contact_normal_force.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  for (int i=0; i<active_contact_indices.size(); ++i) {
    EXPECT_DOUBLE_EQ(kkt_residual.lf().coeff(3*i  ), 0);
    EXPECT_DOUBLE_EQ(kkt_residual.lf().coeff(3*i+1), 0);
    EXPECT_DOUBLE_EQ(kkt_residual.lf().coeff(3*i+2), -0.1*dtau);
  }
}


TEST_F(ContactNormalForceTest, condenseSlackAndDual) {
  ContactNormalForce contact_normal_force(robot); 
  const double slack = 0.3;
  const double dual = 0.2;
  data.slack.fill(slack);
  data.dual.fill(dual);
  contact_normal_force.condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd Qff_ref = Eigen::MatrixXd::Zero(contact_status.dimf(), contact_status.dimf());
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(contact_status.dimf());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    Qff_ref.coeffRef(3*i+2, 3*i+2) = dtau * dtau * dual / slack;
    const double residual = - dtau * s.f[active_contact_indices[i]].coeff(2) + data.slack.coeff(active_contact_indices[i]);
    const double duality = slack * dual - barrier;
    EXPECT_DOUBLE_EQ(residual, data.residual.coeff(active_contact_indices[i]));
    EXPECT_DOUBLE_EQ(duality, data.duality.coeff(active_contact_indices[i]));
    lf_ref.coeffRef(3*i+2) = - dtau * (dual * residual - duality) / slack;
  }
  EXPECT_TRUE(Qff_ref.isApprox(kkt_matrix.Qff()));
  EXPECT_TRUE(lf_ref.isApprox(kkt_residual.lf()));
}


TEST_F(ContactNormalForceTest, computeSlackAndDualDirection) {
  ContactNormalForce contact_normal_force(robot); 
  const double slack = 0.3;
  const double dual = 0.2;
  const double residual = 0.1;
  const double duality = - 0.1;
  data.slack.fill(slack);
  data.dual.fill(dual);
  data.residual.fill(residual);
  data.duality.fill(duality);
  contact_normal_force.computeSlackAndDualDirection(robot, data, dtau, s, d);
  Eigen::VectorXd dslack_ref = Eigen::VectorXd::Ones(robot.max_point_contacts());
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Ones(robot.max_point_contacts());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    const int idx = active_contact_indices[i];
    dslack_ref.coeffRef(idx) = dtau * d.df().coeff(3*i+2) - residual;
    ddual_ref.coeffRef(idx) = pdipmfunc::ComputeDualDirection(slack, dual, 
                                                              data.dslack.coeff(idx), 
                                                              data.duality.coeff(idx));
  }
  EXPECT_TRUE(dslack_ref.isApprox(data.dslack));
  EXPECT_TRUE(ddual_ref.isApprox(data.ddual));
}


TEST_F(ContactNormalForceTest, dimc) {
  ContactNormalForce contact_normal_force(robot); 
  EXPECT_EQ(contact_normal_force.dimc(), contact_frames.size());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
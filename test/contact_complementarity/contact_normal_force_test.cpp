#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/constraints/pdipm_func.hpp"
#include "idocp/contact_complementarity/contact_normal_force.hpp"

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
    robot.setContactStatus(is_contact_active);
    for (int i=0; i<contact_frames.size(); ++i) {
      if (is_contact_active[i]) {
        active_contact_indices.push_back(i);
      }
    }
    barrier = 1.0e-04;
    fraction_to_boundary_rate = 0.995;
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    data = ConstraintComponentData(contact_frames.size());
    data_ref = ConstraintComponentData(contact_frames.size());
    s = SplitSolution::Random(robot);
    d = SplitDirection::Random(robot);
    kkt_residual = KKTResidual(robot);
    kkt_matrix = KKTMatrix(robot);
  }

  virtual void TearDown() {
  }

  double barrier, fraction_to_boundary_rate, dtau;
  Eigen::VectorXd slack, dual, dslack, ddual;
  Robot robot;
  std::vector<int> contact_frames, active_contact_indices;
  std::vector<bool> is_contact_active;
  ConstraintComponentData data, data_ref;
  SplitSolution s;
  SplitDirection d;
  KKTResidual kkt_residual;
  KKTMatrix kkt_matrix;
};


TEST_F(ContactNormalForceTest, isFeasible) {
  ContactNormalForce constraint(robot); 
  is_contact_active[0] = true;
  robot.setContactStatus(is_contact_active);
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(2) = -1;
  }
  EXPECT_FALSE(constraint.isFeasible(robot, data, s));
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(2) = 1;
  }
  EXPECT_TRUE(constraint.isFeasible(robot, data, s));
}


TEST_F(ContactNormalForceTest, setSlackAndDual) {
  ContactNormalForce constraint(robot); 
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(2) = 1;
  }
  EXPECT_TRUE(constraint.isFeasible(robot, data, s));
  constraint.setSlackAndDual(robot, data, dtau, s);
  for (int i=0; i<contact_frames.size(); ++i) {
    EXPECT_DOUBLE_EQ(data.slack.coeff(i), s.f[i].coeffRef(2));
    EXPECT_TRUE(data.dual.coeff(i) > 0);
  }
}


TEST_F(ContactNormalForceTest, augmentDualResidual) {
  ContactNormalForce constraint(robot); 
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.lf().isZero());
  data.dual.fill(0.1);
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  for (int i=0; i<active_contact_indices.size(); ++i) {
    EXPECT_DOUBLE_EQ(kkt_residual.lf().coeff(3*i  ), 0);
    EXPECT_DOUBLE_EQ(kkt_residual.lf().coeff(3*i+1), 0);
    EXPECT_DOUBLE_EQ(kkt_residual.lf().coeff(3*i+2), -0.1*dtau);
  }
}


TEST_F(ContactNormalForceTest, condenseSlackAndDual) {
  ContactNormalForce constraint(robot); 
  const double slack = 0.3;
  const double dual = 0.2;
  data.slack.fill(slack);
  data.dual.fill(dual);
  constraint.condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd Qff_ref = Eigen::MatrixXd::Zero(robot.dimf(), robot.dimf());
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(robot.dimf());
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
  ContactNormalForce constraint(robot); 
  const double slack = 0.3;
  const double dual = 0.2;
  const double residual = 0.1;
  const double duality = - 0.1;
  data.slack.fill(slack);
  data.dual.fill(dual);
  data.residual.fill(residual);
  data.duality.fill(duality);
  constraint.computeSlackAndDualDirection(robot, data, dtau, s, d);
  Eigen::VectorXd dslack_ref = Eigen::VectorXd::Zero(robot.max_point_contacts());
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(robot.max_point_contacts());
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


TEST_F(ContactNormalForceTest, residualL1Nrom) {
  ContactNormalForce constraint(robot); 
  const double slack = 0.3;
  const double dual = 0.2;
  data.slack.fill(slack);
  data.dual.fill(dual);
  double norm = 0;
  for (int i=0; i<active_contact_indices.size(); ++i) {
    norm += std::abs(slack - dtau * s.f[active_contact_indices[i]].coeff(2));
    norm += std::abs(slack * dual - barrier);
  }
  EXPECT_DOUBLE_EQ(norm, constraint.residualL1Nrom(robot, data, dtau, s));
}


TEST_F(ContactNormalForceTest, squaredKKTErrorNorm) {
  ContactNormalForce constraint(robot); 
  const double slack = 0.3;
  const double dual = 0.2;
  data.slack.fill(slack);
  data.dual.fill(dual);
  double norm = 0;
  for (int i=0; i<active_contact_indices.size(); ++i) {
    const double residual = slack - dtau * s.f[active_contact_indices[i]].coeff(2);
    const double duality = slack * dual - barrier;
    norm += residual * residual + duality * duality;
  }
  EXPECT_DOUBLE_EQ(norm, constraint.squaredKKTErrorNorm(robot, data, dtau, s));
}


TEST_F(ContactNormalForceTest, maxStepSize) {
  ContactNormalForce constraint(robot); 
  data.slack = Eigen::VectorXd::Random(robot.max_point_contacts()).array().abs();
  data.dual = Eigen::VectorXd::Random(robot.max_point_contacts()).array().abs();
  data.dslack = Eigen::VectorXd::Random(robot.max_point_contacts());
  data.ddual = Eigen::VectorXd::Random(robot.max_point_contacts());
  Eigen::VectorXd slack_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd dual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd dslack_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd ddual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    slack_active.coeffRef(i) = data.slack.coeff(active_contact_indices[i]);
    dual_active.coeffRef(i) = data.dual.coeff(active_contact_indices[i]);
    dslack_active.coeffRef(i) = data.dslack.coeff(active_contact_indices[i]);
    ddual_active.coeffRef(i) = data.ddual.coeff(active_contact_indices[i]);
  }
  const double max_slack_step_size = pdipmfunc::FractionToBoundary(
      active_contact_indices.size(), fraction_to_boundary_rate, slack_active, dslack_active);
  const double max_dual_step_size = pdipmfunc::FractionToBoundary(
      active_contact_indices.size(), fraction_to_boundary_rate, dual_active, ddual_active);
  EXPECT_DOUBLE_EQ(max_slack_step_size, constraint.maxSlackStepSize(data, is_contact_active));
  EXPECT_DOUBLE_EQ(max_dual_step_size, constraint.maxDualStepSize(data, is_contact_active));
}


TEST_F(ContactNormalForceTest, updateSlackAndDual) {
  ContactNormalForce constraint(robot); 
  data.slack = Eigen::VectorXd::Random(robot.max_point_contacts()).array().abs();
  data.dual = Eigen::VectorXd::Random(robot.max_point_contacts()).array().abs();
  data.dslack = Eigen::VectorXd::Random(robot.max_point_contacts());
  data.ddual = Eigen::VectorXd::Random(robot.max_point_contacts());
  data_ref.slack = data.slack;
  data_ref.dual = data.dual;
  data_ref.dslack = data.dslack;
  data_ref.ddual = data.ddual;
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  constraint.updateSlack(data, is_contact_active, step_size);
  constraint.updateDual(data, is_contact_active, step_size);
  for (const auto i : active_contact_indices) {
    data_ref.slack.coeffRef(i) += step_size * data_ref.dslack.coeff(i);
    data_ref.dual.coeffRef(i) += step_size * data_ref.ddual.coeff(i);
  }
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
}


TEST_F(ContactNormalForceTest, dimc) {
  ContactNormalForce constraint(robot); 
  EXPECT_EQ(constraint.dimc(), contact_frames.size());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
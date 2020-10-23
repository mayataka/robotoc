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
#include "idocp/contact_complementarity/distance_to_contact_surface.hpp"

namespace idocp {

class ContactDistanceTest : public ::testing::Test {
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
      if (!is_contact_active[i]) {
        inactive_contact_indices.push_back(i);
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
    robot.updateKinematics(s.q, s.v, s.a);
  }

  virtual void TearDown() {
  }

  double barrier, fraction_to_boundary_rate, dtau;
  Eigen::VectorXd slack, dual, dslack, ddual;
  Robot robot;
  std::vector<int> contact_frames, inactive_contact_indices;
  std::vector<bool> is_contact_active;
  ConstraintComponentData data, data_ref;
  SplitSolution s;
  SplitDirection d;
  KKTResidual kkt_residual;
  KKTMatrix kkt_matrix;
};


TEST_F(ContactDistanceTest, isFeasible) {
  ContactDistance constraint(robot); 
  bool feasible = true;
  for (int i=0; i<contact_frames.size(); ++i) {
    if (!robot.is_contact_active(i)) {
      Eigen::Vector3d frame_position(Eigen::Vector3d::Zero());
      robot.computeContactResidual(i, frame_position);
      const double distance = frame_position.coeff(2);
      if (distance < 0) {
        feasible = false;
      }
    }
  }
  EXPECT_EQ(constraint.isFeasible(robot, data, s), feasible);
}


TEST_F(ContactDistanceTest, setSlackAndDual) {
  ContactDistance constraint(robot); 
  constraint.setSlackAndDual(robot, data, dtau, s);
  bool feasible = true;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d frame_position(Eigen::Vector3d::Zero());
    robot.computeContactResidual(i, frame_position);
    const double distance = frame_position.coeff(2);
    if (distance > 0) {
      EXPECT_DOUBLE_EQ(data.slack.coeff(i), dtau*distance);
      EXPECT_TRUE(data.dual.coeff(i) > 0);
    }
    else {
      EXPECT_TRUE(data.slack.coeff(i) > 0);
      EXPECT_TRUE(data.dual.coeff(i) > 0);
    }
  }
}


TEST_F(ContactDistanceTest, augmentDualResidual) {
  ContactDistance constraint(robot); 
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.lq().isZero());
  data.slack = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  Eigen::MatrixXd distance_derivative = Eigen::MatrixXd::Zero(contact_frames.size(), robot.dimv());
  const int i_active = 0;
  for (const int i : inactive_contact_indices) {
    Eigen::MatrixXd der = Eigen::MatrixXd::Zero(3, robot.dimv());
    robot.computeContactDerivative(i, der);
    distance_derivative.row(i) = dtau * der.row(2);
  }
  Eigen::VectorXd lq_ref = - distance_derivative.transpose() * data.dual;
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
}


TEST_F(ContactDistanceTest, condenseSlackAndDual) {
  ContactDistance constraint(robot); 
  data.slack = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  Eigen::MatrixXd distance_derivative = Eigen::MatrixXd::Zero(contact_frames.size(), robot.dimv());
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  kkt_residual.setZero();
  constraint.condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, kkt_residual);
  for (const int i : inactive_contact_indices) {
    Eigen::MatrixXd der = Eigen::MatrixXd::Zero(3, robot.dimv());
    robot.computeContactDerivative(i, der);
    distance_derivative.row(i) = dtau * der.row(2);
    Eigen::Vector3d res = Eigen::Vector3d::Zero();
    robot.computeContactResidual(i, res);
    data_ref.residual.coeffRef(i) = - dtau * res.coeff(2) + data.slack.coeff(i);
  }
  Eigen::VectorXd dual_per_slack = Eigen::VectorXd::Zero(contact_frames.size());
  dual_per_slack.array() = data.dual.array() / data.slack.array();
  Eigen::MatrixXd Qqq_ref = distance_derivative.transpose() * dual_per_slack.asDiagonal() * distance_derivative;
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(robot.dimv());
}


TEST_F(ContactDistanceTest, computeSlackAndDualDirection) {
  ContactDistance constraint(robot); 
  data.slack = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.residual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.duality= Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  constraint.computeSlackAndDualDirection(robot, data, dtau, s, d);
  Eigen::VectorXd dslack_ref = Eigen::VectorXd::Zero(robot.max_point_contacts());
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(robot.max_point_contacts());
  Eigen::MatrixXd Jac = Eigen::MatrixXd::Zero(robot.max_point_contacts(), robot.dimv());
  for (const int i : inactive_contact_indices) {
    Eigen::MatrixXd der = Eigen::MatrixXd::Zero(3, robot.dimv());
    robot.computeContactDerivative(i, der);
    Jac.row(i) = der.row(2);
  }
  const Eigen::VectorXd dslack_ref_full = dtau * Jac * d.dq() - data.residual;
  for (const int i : inactive_contact_indices) {
    dslack_ref.coeffRef(i) = dslack_ref_full.coeff(i);
    ddual_ref.coeffRef(i) = pdipmfunc::ComputeDualDirection(data.slack.coeff(i),  
                                                            data.dual.coeff(i), 
                                                            data.dslack.coeff(i), 
                                                            data.duality.coeff(i));
  }
  EXPECT_TRUE(dslack_ref.isApprox(data.dslack));
  EXPECT_TRUE(ddual_ref.isApprox(data.ddual));
}


TEST_F(ContactDistanceTest, residualL1Nrom) {
  ContactDistance constraint(robot); 
  data.slack = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.residual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.duality= Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  double norm = 0;
  for (const int i : inactive_contact_indices) {
    Eigen::Vector3d res = Eigen::Vector3d::Zero();
    robot.computeContactResidual(i, res);
    const double distance = res.coeff(2);
    norm += std::abs(data.slack.coeff(i) - dtau * distance);
    norm += std::abs(data.slack.coeff(i) * data.dual.coeff(i) - barrier);
  }
  EXPECT_DOUBLE_EQ(norm, constraint.residualL1Nrom(robot, data, dtau, s));
}


TEST_F(ContactDistanceTest, squaredKKTErrorNorm) {
  ContactDistance constraint(robot); 
  data.slack = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.residual = Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  data.duality= Eigen::VectorXd::Random(contact_frames.size()).array().abs();
  double norm = 0;
  for (const int i : inactive_contact_indices) {
    Eigen::Vector3d res = Eigen::Vector3d::Zero();
    robot.computeContactResidual(i, res);
    const double distance = res.coeff(2);
    const double residual = data.slack.coeff(i) - dtau * distance;
    const double duality = data.slack.coeff(i) * data.dual.coeff(i) - barrier;
    norm += residual * residual + duality * duality;
  }
  EXPECT_DOUBLE_EQ(norm, constraint.squaredKKTErrorNorm(robot, data, dtau, s));
}


TEST_F(ContactDistanceTest, maxStepSize) {
  ContactDistance constraint(robot); 
  data.slack = Eigen::VectorXd::Random(robot.max_point_contacts()).array().abs();
  data.dual = Eigen::VectorXd::Random(robot.max_point_contacts()).array().abs();
  data.dslack = Eigen::VectorXd::Random(robot.max_point_contacts());
  data.ddual = Eigen::VectorXd::Random(robot.max_point_contacts());
  Eigen::VectorXd slack_active = Eigen::VectorXd::Zero(inactive_contact_indices.size());
  Eigen::VectorXd dual_active = Eigen::VectorXd::Zero(inactive_contact_indices.size());
  Eigen::VectorXd dslack_active = Eigen::VectorXd::Zero(inactive_contact_indices.size());
  Eigen::VectorXd ddual_active = Eigen::VectorXd::Zero(inactive_contact_indices.size());
  for (int i=0; i<inactive_contact_indices.size(); ++i) {
    slack_active.coeffRef(i) = data.slack.coeff(inactive_contact_indices[i]);
    dual_active.coeffRef(i) = data.dual.coeff(inactive_contact_indices[i]);
    dslack_active.coeffRef(i) = data.dslack.coeff(inactive_contact_indices[i]);
    ddual_active.coeffRef(i) = data.ddual.coeff(inactive_contact_indices[i]);
  }
  double max_slack_step_size, max_dual_step_size;
  if (robot.dimf() < robot.max_dimf()) {
    max_slack_step_size = pdipmfunc::FractionToBoundary(
        inactive_contact_indices.size(), fraction_to_boundary_rate, slack_active, dslack_active);
    max_dual_step_size = pdipmfunc::FractionToBoundary(
        inactive_contact_indices.size(), fraction_to_boundary_rate, dual_active, ddual_active);
  }
  else {
    max_slack_step_size = 1;
    max_dual_step_size = 1;
  }
  EXPECT_DOUBLE_EQ(max_slack_step_size, constraint.maxSlackStepSize(data, is_contact_active));
  EXPECT_DOUBLE_EQ(max_dual_step_size, constraint.maxDualStepSize(data, is_contact_active));
}


TEST_F(ContactDistanceTest, updateSlackAndDual) {
  ContactDistance constraint(robot); 
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
  for (const auto i : inactive_contact_indices) {
    data_ref.slack.coeffRef(i) += step_size * data_ref.dslack.coeff(i);
    data_ref.dual.coeffRef(i) += step_size * data_ref.ddual.coeff(i);
  }
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
}


TEST_F(ContactDistanceTest, dimc) {
  ContactDistance constraint(robot); 
  EXPECT_EQ(constraint.dimc(), contact_frames.size());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
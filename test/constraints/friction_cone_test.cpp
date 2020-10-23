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
#include "idocp/constraints/pdipm.hpp"
#include "idocp/constraints/friction_cone.hpp"

namespace idocp {

class FrictionConeTest : public ::testing::Test {
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


TEST_F(FrictionConeTest, isFeasible) {
  FrictionCone friction_cone(robot); 
  contact_status.activateContact(0);
  s.setContactStatus(contact_status);
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(0) = 0;
    s.f[i].coeffRef(1) = 0;
    s.f[i].coeffRef(2) = 1;
  }
  EXPECT_FALSE(friction_cone.isFeasible(robot, data, s));
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(0) = 1;
    s.f[i].coeffRef(1) = 1;
    s.f[i].coeffRef(2) = 0;
  }
  EXPECT_TRUE(friction_cone.isFeasible(robot, data, s));
}


TEST_F(FrictionConeTest, frictionConeResidual) {
  FrictionCone friction_cone(robot); 
  const double res = friction_cone.frictionConeResidual(mu, s.f[0]);
  const double fx = s.f[0].coeff(0);
  const double fy = s.f[0].coeff(1);
  const double fz = s.f[0].coeff(2);
  const double res_ref = fx * fx + fy * fy - mu * mu * fz * fz;
  EXPECT_DOUBLE_EQ(res, res_ref);
}


TEST_F(FrictionConeTest, setSlackAndDual) {
  FrictionCone friction_cone(robot); 
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].fill(1);
  }
  EXPECT_TRUE(friction_cone.isFeasible(robot, data, s));
  friction_cone.setSlackAndDual(robot, data, dtau, s);
  for (int i=0; i<contact_frames.size(); ++i) {
    EXPECT_DOUBLE_EQ(data.slack.coeff(i), dtau*(1+1-mu*mu));
    EXPECT_TRUE(data.dual.coeff(i) > 0);
  }
}


TEST_F(FrictionConeTest, augmentDualResidual) {
  FrictionCone friction_cone(robot); 
  friction_cone.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.lf().isZero());
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  friction_cone.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  Eigen::MatrixXd gf = Eigen::MatrixXd::Zero(active_contact_indices.size(), contact_status.dimf());
  Eigen::VectorXd dual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    const double fx = s.f[active_contact_indices[i]].coeff(0);
    const double fy = s.f[active_contact_indices[i]].coeff(1);
    const double fz = s.f[active_contact_indices[i]].coeff(2);
    const double mu = 0.8;
    gf.row(i).coeffRef(3*i+0) = 2 * dtau * fx;
    gf.row(i).coeffRef(3*i+1) = 2 * dtau * fy;
    gf.row(i).coeffRef(3*i+2) = - 2 * mu * mu * dtau * fz;
    dual_active.coeffRef(i) = data.dual(active_contact_indices[i]);
  }
  std::cout << gf << std::endl;
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(contact_status.dimf());
  lf_ref = gf.transpose() * dual_active;
  EXPECT_TRUE(kkt_residual.lf().isApprox(lf_ref));
}


TEST_F(FrictionConeTest, condenseSlackAndDual) {
  FrictionCone friction_cone(robot); 
  data.slack = Eigen::VectorXd::Random(data.slack.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  friction_cone.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  kkt_residual.setZero();
  friction_cone.condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd gf = Eigen::MatrixXd::Zero(active_contact_indices.size(), contact_status.dimf());
  Eigen::VectorXd slack_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd dual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd residual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd duality_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    const double fx = s.f[active_contact_indices[i]].coeff(0);
    const double fy = s.f[active_contact_indices[i]].coeff(1);
    const double fz = s.f[active_contact_indices[i]].coeff(2);
    const double mu = 0.8;
    gf.row(i).coeffRef(3*i+0) = 2 * dtau * fx;
    gf.row(i).coeffRef(3*i+1) = 2 * dtau * fy;
    gf.row(i).coeffRef(3*i+2) = - 2 * mu * mu * dtau * fz;
    slack_active.coeffRef(i) = data.slack(active_contact_indices[i]);
    dual_active.coeffRef(i) = data.dual(active_contact_indices[i]);
    residual_active.coeffRef(i) = - dtau * (fx*fx+fy*fy-mu*mu*fz*fz) + slack_active.coeff(i);
    duality_active.coeffRef(i) = slack_active.coeff(i) * dual_active.coeff(i) - barrier;
  }
  Eigen::VectorXd dual_per_slack = Eigen::VectorXd::Zero(active_contact_indices.size());
  dual_per_slack.array() = dual_active.array() / slack_active.array();
  Eigen::MatrixXd Qff_ref = gf.transpose() * dual_per_slack.asDiagonal() * gf;
  EXPECT_TRUE(Qff_ref.isApprox(kkt_matrix.Qff()));
  Eigen::VectorXd condensed_residual = Eigen::VectorXd::Zero(active_contact_indices.size());
  condensed_residual.array() 
      = (dual_active.array() * residual_active.array() - duality_active.array())
        / slack_active.array();
  Eigen::VectorXd lf_ref = gf.transpose() * condensed_residual;
  EXPECT_TRUE(lf_ref.isApprox(kkt_residual.lf()));
  std::cout << lf_ref.transpose() << std::endl;
  std::cout << kkt_residual.lf().transpose() << std::endl;
}


TEST_F(FrictionConeTest, computeSlackAndDualDirection) {
  FrictionCone friction_cone(robot); 
  data.slack = Eigen::VectorXd::Random(data.slack.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  data.residual = Eigen::VectorXd::Random(data.dual.size());
  data.duality = Eigen::VectorXd::Random(data.dual.size());
  friction_cone.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  friction_cone.condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, kkt_residual);
  friction_cone.computeSlackAndDualDirection(robot, data, dtau, s, d);
  Eigen::MatrixXd gf = Eigen::MatrixXd::Zero(active_contact_indices.size(), contact_status.dimf());
  Eigen::VectorXd slack_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd dual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd residual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd duality_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    const double fx = s.f[active_contact_indices[i]].coeff(0);
    const double fy = s.f[active_contact_indices[i]].coeff(1);
    const double fz = s.f[active_contact_indices[i]].coeff(2);
    const double mu = 0.8;
    gf.row(i).coeffRef(3*i+0) = 2 * dtau * fx;
    gf.row(i).coeffRef(3*i+1) = 2 * dtau * fy;
    gf.row(i).coeffRef(3*i+2) = - 2 * mu * mu * dtau * fz;
    slack_active.coeffRef(i) = data.slack(active_contact_indices[i]);
    dual_active.coeffRef(i) = data.dual(active_contact_indices[i]);
    residual_active.coeffRef(i) = - dtau * (fx*fx+fy*fy-mu*mu*fz*fz) + slack_active.coeff(i);
    duality_active.coeffRef(i) = slack_active.coeff(i) * dual_active.coeff(i) - barrier;
  }
  Eigen::VectorXd dslack_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd ddual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  dslack_active = - gf * d.df() - residual_active;

  pdipmfunc::ComputeDualDirection(slack_active, dual_active, dslack_active, 
                                  duality_active, ddual_active);
  Eigen::VectorXd dslack_ref = Eigen::VectorXd::Ones(robot.max_point_contacts());
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Ones(robot.max_point_contacts());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    dslack_ref.coeffRef(active_contact_indices[i]) = dslack_active.coeff(i);
    ddual_ref.coeffRef(active_contact_indices[i]) = ddual_active.coeff(i);
  }
  EXPECT_TRUE(dslack_ref.isApprox(data.dslack));
  EXPECT_TRUE(ddual_ref.isApprox(data.ddual));
}


TEST_F(FrictionConeTest, dimc) {
  FrictionCone friction_cone(robot); 
  EXPECT_EQ(friction_cone.dimc(), contact_frames.size());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
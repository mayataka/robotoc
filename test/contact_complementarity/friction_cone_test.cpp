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
#include "idocp/contact_complementarity/friction_cone.hpp"

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


TEST_F(FrictionConeTest, isFeasible) {
  FrictionCone constraint(robot); 
  is_contact_active[0] = true;
  robot.setContactStatus(is_contact_active);
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(0) = 0;
    s.f[i].coeffRef(1) = 0;
    s.f[i].coeffRef(2) = 1;
  }
  EXPECT_FALSE(constraint.isFeasible(robot, data, s));
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].coeffRef(0) = 1;
    s.f[i].coeffRef(1) = 1;
    s.f[i].coeffRef(2) = 0;
  }
  EXPECT_TRUE(constraint.isFeasible(robot, data, s));
}


TEST_F(FrictionConeTest, frictionConeResidual) {
  FrictionCone constraint(robot); 
  const double mu = 0.8;
  const double res = constraint.frictionConeResidual(mu, s.f[0]);
  const double fx = s.f[0].coeff(0);
  const double fy = s.f[0].coeff(1);
  const double fz = s.f[0].coeff(2);
  const double res_ref = fx * fx + fy * fy - mu * mu * fz * fz;
  EXPECT_DOUBLE_EQ(res, res_ref);
}


TEST_F(FrictionConeTest, setSlackAndDual) {
  FrictionCone constraint(robot); 
  for (int i=0; i<contact_frames.size(); ++i) {
    s.f[i].fill(1);
  }
  EXPECT_TRUE(constraint.isFeasible(robot, data, s));
  constraint.setSlackAndDual(robot, data, dtau, s);
  for (int i=0; i<contact_frames.size(); ++i) {
    EXPECT_DOUBLE_EQ(data.slack.coeff(i), dtau*(1+1-0.8*0.8));
    EXPECT_TRUE(data.dual.coeff(i) > 0);
  }
}


TEST_F(FrictionConeTest, augmentDualResidual) {
  FrictionCone constraint(robot); 
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  EXPECT_TRUE(kkt_residual.lf().isZero());
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  Eigen::MatrixXd gf = Eigen::MatrixXd::Zero(active_contact_indices.size(), robot.dimf());
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
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(robot.dimf());
  lf_ref = gf.transpose() * dual_active;
  EXPECT_TRUE(kkt_residual.lf().isApprox(lf_ref));
}


TEST_F(FrictionConeTest, condenseSlackAndDual) {
  FrictionCone constraint(robot); 
  data.slack = Eigen::VectorXd::Random(data.slack.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  kkt_residual.setZero();
  constraint.condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd gf = Eigen::MatrixXd::Zero(active_contact_indices.size(), robot.dimf());
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
}


TEST_F(FrictionConeTest, computeSlackAndDualDirection) {
  FrictionCone constraint(robot); 
  data.slack = Eigen::VectorXd::Random(data.slack.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  data.residual = Eigen::VectorXd::Random(data.dual.size());
  data.duality = Eigen::VectorXd::Random(data.dual.size());
  constraint.augmentDualResidual(robot, data, dtau, s, kkt_residual);
  constraint.condenseSlackAndDual(robot, data, dtau, s, kkt_matrix, kkt_residual);
  constraint.computeSlackAndDualDirection(robot, data, dtau, s, d);
  Eigen::MatrixXd gf = Eigen::MatrixXd::Zero(active_contact_indices.size(), robot.dimf());
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
  Eigen::VectorXd dslack_ref = Eigen::VectorXd::Zero(robot.max_point_contacts());
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(robot.max_point_contacts());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    dslack_ref.coeffRef(active_contact_indices[i]) = dslack_active.coeff(i);
    ddual_ref.coeffRef(active_contact_indices[i]) = ddual_active.coeff(i);
  }
  EXPECT_TRUE(dslack_ref.isApprox(data.dslack));
  EXPECT_TRUE(ddual_ref.isApprox(data.ddual));
}


TEST_F(FrictionConeTest, residualL1Nrom) {
  FrictionCone constraint(robot); 
  data.slack = Eigen::VectorXd::Random(data.slack.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  data.residual = Eigen::VectorXd::Random(data.dual.size());
  data.duality = Eigen::VectorXd::Random(data.dual.size());
  const double norm = constraint.residualL1Nrom(robot, data, dtau, s);
  Eigen::VectorXd slack_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd dual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd residual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd duality_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    const double fx = s.f[active_contact_indices[i]].coeff(0);
    const double fy = s.f[active_contact_indices[i]].coeff(1);
    const double fz = s.f[active_contact_indices[i]].coeff(2);
    const double mu = 0.8;
    slack_active.coeffRef(i) = data.slack(active_contact_indices[i]);
    dual_active.coeffRef(i) = data.dual(active_contact_indices[i]);
    residual_active.coeffRef(i) = dtau * (fx*fx+fy*fy-mu*mu*fz*fz) + slack_active.coeff(i);
    duality_active.coeffRef(i) = slack_active.coeff(i) * dual_active.coeff(i) - barrier;
  }
  double norm_ref = residual_active.template lpNorm<1>();
  norm_ref += duality_active.template lpNorm<1>();
  EXPECT_DOUBLE_EQ(norm, norm_ref);
}


TEST_F(FrictionConeTest, squaredKKTErrorNorm) {
  FrictionCone constraint(robot); 
  data.slack = Eigen::VectorXd::Random(data.slack.size()).array().abs();
  data.dual = Eigen::VectorXd::Random(data.dual.size()).array().abs();
  data.residual = Eigen::VectorXd::Random(data.dual.size());
  data.duality = Eigen::VectorXd::Random(data.dual.size());
  const double norm = constraint.squaredKKTErrorNorm(robot, data, dtau, s);
  Eigen::VectorXd slack_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd dual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd residual_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  Eigen::VectorXd duality_active = Eigen::VectorXd::Zero(active_contact_indices.size());
  for (int i=0; i<active_contact_indices.size(); ++i) {
    const double fx = s.f[active_contact_indices[i]].coeff(0);
    const double fy = s.f[active_contact_indices[i]].coeff(1);
    const double fz = s.f[active_contact_indices[i]].coeff(2);
    const double mu = 0.8;
    slack_active.coeffRef(i) = data.slack(active_contact_indices[i]);
    dual_active.coeffRef(i) = data.dual(active_contact_indices[i]);
    residual_active.coeffRef(i) = dtau * (fx*fx+fy*fy-mu*mu*fz*fz) + slack_active.coeff(i);
    duality_active.coeffRef(i) = slack_active.coeff(i) * dual_active.coeff(i) - barrier;
  }
  double norm_ref = residual_active.squaredNorm();
  norm_ref += duality_active.squaredNorm();
  EXPECT_DOUBLE_EQ(norm, norm_ref);
}


TEST_F(FrictionConeTest, maxStepSize) {
  FrictionCone constraint(robot); 
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
  double max_slack_step_size, max_dual_step_size;
  if (robot.dimf() > 0) {
    max_slack_step_size = pdipmfunc::FractionToBoundary(
        active_contact_indices.size(), fraction_to_boundary_rate, slack_active, dslack_active);
    max_dual_step_size = pdipmfunc::FractionToBoundary(
        active_contact_indices.size(), fraction_to_boundary_rate, dual_active, ddual_active);
  }
  else {
    max_slack_step_size = 1;
    max_dual_step_size = 1;
  }
  EXPECT_DOUBLE_EQ(max_slack_step_size, constraint.maxSlackStepSize(data, is_contact_active));
  EXPECT_DOUBLE_EQ(max_dual_step_size, constraint.maxDualStepSize(data, is_contact_active));
}


TEST_F(FrictionConeTest, updateSlackAndDual) {
  FrictionCone constraint(robot); 
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


TEST_F(FrictionConeTest, dimc) {
  FrictionCone constraint(robot); 
  EXPECT_EQ(constraint.dimc(), contact_frames.size());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
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
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/pdipm.hpp"
#include "idocp/constraints/contact_distance.hpp"

namespace idocp {

class ContactDistanceTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    barrier = 1.0e-04;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    fraction_to_boundary_rate = 0.995;
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot, const ContactStatus& contact_status) const;
  void testIsFeasible(Robot& robot, const ContactStatus& contact_status) const;
  void testSetSlackAndDual(Robot& robot, const ContactStatus& contact_status) const;
  void testAugmentDualResidual(Robot& robot, const ContactStatus& contact_status) const;
  void testComputePrimalAndDualResidual(Robot& robot, const ContactStatus& contact_status) const;
  void testCondenseSlackAndDual(Robot& robot, const ContactStatus& contact_status) const;
  void testComputeSlackAndDualDirection(Robot& robot, const ContactStatus& contact_status) const;

  double barrier, dt, fraction_to_boundary_rate;
  std::string fixed_base_urdf, floating_base_urdf;
};


void ContactDistanceTest::testKinematics(Robot& robot, const ContactStatus& contact_status) const {
  ContactDistance limit(robot); 
  EXPECT_TRUE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::PositionLevel);
}


void ContactDistanceTest::testIsFeasible(Robot& robot, const ContactStatus& contact_status) const {
  ContactDistance limit(robot); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), contact_status.maxPointContacts());
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  robot.updateFrameKinematics(s.q);
  bool is_feasible_ref = true;
  for (int i=0; i<limit.dimc(); ++i) {
    if (!contact_status.isContactActive(i)) {
      if (robot.framePosition(robot.contactFramesIndices()[i]).coeff(2) <= 0) {
        is_feasible_ref = false;
      }
    }
  }
  EXPECT_EQ(limit.isFeasible(robot, data, s), is_feasible_ref);
}


void ContactDistanceTest::testSetSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  ContactDistance limit(robot); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<dimc; ++i) {
    data_ref.slack.coeffRef(i) = robot.framePosition(robot.contactFramesIndices()[i]).coeff(2);
  }
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void ContactDistanceTest::testAugmentDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  ContactDistance limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.lq().setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  robot.updateKinematics(s.q, s.v, s.a);
  limit.allocateExtraData(data);
  limit.augmentDualResidual(robot, data, dt, s, kkt_res);
  std::vector<Eigen::MatrixXd> J(dimc, Eigen::MatrixXd::Zero(6, robot.dimv()));
  for (int i=0; i<dimc; ++i) {
    robot.getFrameJacobian(robot.contactFramesIndices()[i], J[i]);
  }
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (!contact_status.isContactActive(i)) {
      const Eigen::VectorXd Jrow2 = J[i].row(2);
      kkt_res_ref.lq() -= dt * data_ref.dual.coeff(i) * Jrow2;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void ContactDistanceTest::testComputePrimalAndDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  ContactDistance limit(robot); 
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ConstraintComponentData data(limit.dimc());
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, s);
  robot.updateFrameKinematics(s.q);
  data_ref.residual.setZero();
  data_ref.duality.setZero();
  for (int i=0; i<dimc; ++i) {
    if (!contact_status.isContactActive(i)) {
      data_ref.residual.coeffRef(i) = - robot.framePosition(robot.contactFramesIndices()[i]).coeff(2) + data_ref.slack.coeff(i);
      data_ref.duality.coeffRef(i) = data_ref.slack.coeff(i) * data_ref.dual.coeff(i) - barrier;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void ContactDistanceTest::testCondenseSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  ContactDistance limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.setContactStatus(contact_status);
  kkt_mat.Qqq().setRandom();
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.lq().setRandom();
  robot.updateKinematics(s.q, s.v, s.a);
  limit.allocateExtraData(data);
  limit.augmentDualResidual(robot, data, dt, s, kkt_res);
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  std::vector<Eigen::MatrixXd> J(dimc, Eigen::MatrixXd::Zero(6, robot.dimv()));
  for (int i=0; i<dimc; ++i) {
    robot.getFrameJacobian(robot.contactFramesIndices()[i], J[i]);
  }
  for (int i=0; i<dimc; ++i) {
    if (!contact_status.isContactActive(i)) {
      const Eigen::VectorXd Jrow2 = J[i].row(2);
      kkt_mat_ref.Qqq() += (dt * data_ref.dual.coeff(i) / data_ref.slack.coeff(i)) 
                            * Jrow2 * Jrow2.transpose();
      data_ref.residual.coeffRef(i) 
          = - robot.framePosition(robot.contactFramesIndices()[i]).coeff(2) 
              + data_ref.slack.coeff(i);
      data_ref.duality.coeffRef(i) = data_ref.slack.coeff(i) * data_ref.dual.coeff(i) - barrier;
      kkt_res_ref.lq().noalias()
          -= (dt * (data_ref.dual.coeff(i)*data_ref.residual.coeff(i)-data_ref.duality.coeff(i)) 
                  / data_ref.slack.coeff(i))
              * Jrow2;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void ContactDistanceTest::testComputeSlackAndDualDirection(Robot& robot, const ContactStatus& contact_status) const {
  ContactDistance limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  limit.allocateExtraData(data);
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  limit.augmentDualResidual(robot, data, dt, s, kkt_res);
  data.residual.setRandom();
  data.duality.setRandom();
  ConstraintComponentData data_ref = data;
  const SplitDirection d = SplitDirection::Random(robot, contact_status);
  limit.computeSlackAndDualDirection(robot, data, s, d);
  std::vector<Eigen::MatrixXd> J(dimc, Eigen::MatrixXd::Zero(6, robot.dimv()));
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  for (int i=0; i<dimc; ++i) {
    robot.getFrameJacobian(robot.contactFramesIndices()[i], J[i]);
  }
  for (int i=0; i<dimc; ++i) {
    if (!contact_status.isContactActive(i)) {
      const Eigen::VectorXd Jrow2 = J[i].row(2);
      data_ref.dslack.coeffRef(i) = Jrow2.dot(d.dq()) - data_ref.residual.coeff(i);
      data_ref.ddual.coeffRef(i) = - (data_ref.dual.coeff(i)*data_ref.dslack.coeff(i)+data_ref.duality.coeff(i))
                                      / data_ref.slack.coeff(i);
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(ContactDistanceTest, fixedBase) {
  const std::vector<int> frames = {18};
  Robot robot(fixed_base_urdf, frames);
  auto contact_status = robot.createContactStatus();
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
  contact_status.activateContact(0);
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
}


TEST_F(ContactDistanceTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, frames);
  auto contact_status = robot.createContactStatus();
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
  contact_status.setRandom();
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
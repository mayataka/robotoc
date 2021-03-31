#include "kkt_factory.hpp"


namespace idocp {
namespace testhelper {

SplitKKTMatrix CreateSplitKKTMatrix(const Robot& robot, const double dt) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const Eigen::MatrixXd seed = Eigen::MatrixXd::Random(2*dimv+dimu, 2*dimv+dimu);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qss() = seed * seed.transpose();
  kkt_matrix.Qvq().setZero();
  kkt_matrix.Quq().setZero();
  kkt_matrix.Quv().setZero();
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fqv() = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
    kkt_matrix.Fqv().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  return kkt_matrix;
}


SplitKKTResidual CreateSplitKKTResidual(const Robot& robot) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lx().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  return kkt_residual;
}


SplitKKTResidual CreateSplitKKTResidual(const Robot& robot, 
                                        const ImpulseStatus& impulse_status) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.lx().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  kkt_residual.P().setRandom();
  return kkt_residual;
}


KKTMatrix CreateKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence, 
                          const int N, const int max_num_impulse, const bool is_parnmpc) {
  KKTMatrix kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  for (int i=0; i<=N; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.hasFloatingBase()) {
      kkt_matrix[i].Fqq().setIdentity();
      kkt_matrix[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix[i].Fvq().setRandom();
    kkt_matrix[i].Fvv().setRandom();
    kkt_matrix[i].Fvu().setRandom();
  }
  const int num_impulse = contact_sequence.numImpulseEvents();
  for (int i=0; i<num_impulse; ++i) {
    kkt_matrix.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx, dimx);
    kkt_matrix.impulse[i].Qxx() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    if (robot.hasFloatingBase()) {
      kkt_matrix.impulse[i].Fqq().setIdentity();
      kkt_matrix.impulse[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.impulse[i].Fvq().setRandom();
    kkt_matrix.impulse[i].Fvv().setRandom();
  }
  for (int i=0; i<num_impulse; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix.aux[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix.aux[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix.aux[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.hasFloatingBase()) {
      kkt_matrix.aux[i].Fqq().setIdentity();
      kkt_matrix.aux[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.aux[i].Fvq().setRandom();
    kkt_matrix.aux[i].Fvv().setRandom();
    kkt_matrix.aux[i].Fvu().setRandom();
  }
  const int num_lift = contact_sequence.numLiftEvents();
  for (int i=0; i<num_lift; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix.lift[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix.lift[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix.lift[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.hasFloatingBase()) {
      kkt_matrix.lift[i].Fqq().setIdentity();
      kkt_matrix.lift[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.lift[i].Fvq().setRandom();
    kkt_matrix.lift[i].Fvv().setRandom();
    kkt_matrix.lift[i].Fvu().setRandom();
  }
  if (is_parnmpc) {
    for (int i=0; i<num_impulse; ++i) {
      kkt_matrix.aux[i].setImpulseStatus(contact_sequence.impulseStatus(i));
      kkt_matrix.aux[i].Pq().setRandom();
    }
  }
  else {
  }
  return kkt_matrix;
}


KKTResidual CreateKKTResidual(const Robot& robot, const ContactSequence& contact_sequence, 
                              const int N, const int max_num_impulse, const bool is_parnmpc) {
  KKTResidual kkt_residual = KKTResidual(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    kkt_residual[i].lx().setRandom();
    kkt_residual[i].lu().setRandom();
    kkt_residual[i].Fx().setRandom();
  }
  const int num_impulse = contact_sequence.numImpulseEvents();
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    kkt_residual.impulse[i].lx().setRandom();
    kkt_residual.impulse[i].Fx().setRandom();
  }
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.aux[i].lx().setRandom();
    kkt_residual.aux[i].lu().setRandom();
    kkt_residual.aux[i].Fx().setRandom();
  }
  const int num_lift = contact_sequence.numLiftEvents();
  for (int i=0; i<num_lift; ++i) {
    kkt_residual.lift[i].lx().setRandom();
    kkt_residual.lift[i].lu().setRandom();
    kkt_residual.lift[i].Fx().setRandom();
  }
  if (is_parnmpc) {
    for (int i=0; i<num_impulse; ++i) {
      kkt_residual.aux[i].setImpulseStatus(contact_sequence.impulseStatus(i));
      kkt_residual.aux[i].P().setRandom();
    }
  }
  else {

  }
  return kkt_residual;
}

} // namespace testhelper
} // namespace idocp
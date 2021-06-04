#include "kkt_factory.hpp"


namespace idocp {
namespace testhelper {

SplitKKTMatrix CreateSplitKKTMatrix(const Robot& robot, const double dt) {
  const int dimv = robot.dimv();
  const int dimx = 2*dimv;
  const int dimu = robot.dimu();
  SplitKKTMatrix kkt_matrix(robot);
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
    kkt_matrix.Fqq().topLeftCorner(6, 6).setRandom();
    kkt_matrix.Fqv() = dt * Eigen::MatrixXd::Identity(dimv, dimv);
    kkt_matrix.Fqv().topLeftCorner(6, 6).setRandom();
  }
  else {
    kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
    kkt_matrix.Fqv() = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  }
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu.setRandom();
  const Eigen::MatrixXd H_seed = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
  const Eigen::MatrixXd H = H_seed * H_seed.transpose();
  kkt_matrix.Qxx = H.topLeftCorner(dimx, dimx);
  kkt_matrix.Qxu = H.topRightCorner(dimx, dimu);
  kkt_matrix.Quu = H.bottomRightCorner(dimu, dimu);
  return kkt_matrix;
}


ImpulseSplitKKTMatrix CreateImpulseSplitKKTMatrix(const Robot& robot) {
  const int dimv = robot.dimv();
  const int dimx = 2*dimv;
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  else {
    kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  }
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  const Eigen::MatrixXd Qxx_seed = Eigen::MatrixXd::Random(dimx, dimx);
  kkt_matrix.Qxx = Qxx_seed * Qxx_seed.transpose();
  return kkt_matrix;
}


SplitKKTResidual CreateSplitKKTResidual(const Robot& robot) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.Fx.setRandom();
  kkt_residual.lx.setRandom();
  kkt_residual.lu.setRandom();
  return kkt_residual;
}


ImpulseSplitKKTResidual CreateImpulseSplitKKTResidual(const Robot& robot) {
  ImpulseSplitKKTResidual  kkt_residual(robot);
  kkt_residual.Fx.setRandom();
  kkt_residual.lx.setRandom();
  return kkt_residual;
}


KKTMatrix CreateKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence, 
                          const int N, const int max_num_impulse) {
  KKTMatrix kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const double dt = 0.01;
  for (int i=0; i<=N; ++i) {
    kkt_matrix[i] = CreateSplitKKTMatrix(robot, dt);
  }
  const int num_impulse = contact_sequence.numImpulseEvents();
  for (int i=0; i<num_impulse; ++i) {
    kkt_matrix.impulse[i] = CreateImpulseSplitKKTMatrix(robot);
  }
  for (int i=0; i<num_impulse; ++i) {
    kkt_matrix.aux[i] = CreateSplitKKTMatrix(robot, dt);
  }
  const int num_lift = contact_sequence.numLiftEvents();
  for (int i=0; i<num_lift; ++i) {
    kkt_matrix.lift[i] = CreateSplitKKTMatrix(robot, dt);
  }
  for (int i=0; i<num_impulse; ++i) {
    kkt_matrix.switching[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    kkt_matrix.switching[i].Pq().setRandom();
    kkt_matrix.switching[i].Phix().setRandom();
    kkt_matrix.switching[i].Phia().setRandom();
  }
  return kkt_matrix;
}


KKTResidual CreateKKTResidual(const Robot& robot, const ContactSequence& contact_sequence, 
                              const int N, const int max_num_impulse) {
  KKTResidual kkt_residual = KKTResidual(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    kkt_residual[i] = CreateSplitKKTResidual(robot);
  }
  const int num_impulse = contact_sequence.numImpulseEvents();
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.impulse[i] = CreateImpulseSplitKKTResidual(robot);
  }
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.aux[i] = CreateSplitKKTResidual(robot);
  }
  const int num_lift = contact_sequence.numLiftEvents();
  for (int i=0; i<num_lift; ++i) {
    kkt_residual.lift[i] = CreateSplitKKTResidual(robot);
  }
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.switching[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    kkt_residual.switching[i].P().setRandom();
  }
  return kkt_residual;
}

} // namespace testhelper
} // namespace idocp
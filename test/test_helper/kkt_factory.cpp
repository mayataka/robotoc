#include "kkt_factory.hpp"


namespace robotoc {
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
  kkt_matrix.hx.setRandom();
  kkt_matrix.hu.setRandom();
  kkt_matrix.ha.setRandom();
  kkt_matrix.hf().setRandom();
  kkt_matrix.fx.setRandom();
  return kkt_matrix;
}


SplitKKTMatrix CreateSplitKKTMatrix(const Robot& robot) {
  const int dimv = robot.dimv();
  const int dimx = 2*dimv;
  SplitKKTMatrix kkt_matrix(robot);
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
  kkt_residual.h = Eigen::VectorXd::Random(1)[0];
  return kkt_residual;
}


KKTMatrix CreateKKTMatrix(const Robot& robot, const int N) {
  KKTMatrix kkt_matrix(N+1, SplitKKTMatrix(robot));
  for (int i=0; i<=N; ++i) {
    kkt_matrix[i] = CreateSplitKKTMatrix(robot);
  }
  return kkt_matrix;
}


KKTMatrix CreateKKTMatrix(const Robot& robot,  
                          const std::shared_ptr<ContactSequence>& contact_sequence, 
                          const TimeDiscretization& time_discretization) {
  KKTMatrix kkt_matrix(time_discretization.N_grids()+1, SplitKKTMatrix(robot));
  for (int i=0; i<=time_discretization.N_grids(); ++i) {
    kkt_matrix[i] = CreateSplitKKTMatrix(robot);
    if (time_discretization.grid(i).switching_constraint) {
      const int impulse_index = time_discretization.grid(i).impulse_index + 1;
      kkt_matrix[i].setSwitchingConstraintDimension(contact_sequence->impulseStatus(impulse_index).dimf());
      kkt_matrix[i].Phix().setRandom();
      kkt_matrix[i].Phia().setRandom();
    }
  }
  return kkt_matrix;
}


KKTResidual CreateKKTResidual(const Robot& robot, const int N) {
  KKTResidual kkt_residual(N+1, SplitKKTResidual(robot));
  for (int i=0; i<=N; ++i) {
    kkt_residual[i] = CreateSplitKKTResidual(robot);
  }
  return kkt_residual;
}


KKTResidual CreateKKTResidual(const Robot& robot,  
                              const std::shared_ptr<ContactSequence>& contact_sequence, 
                              const TimeDiscretization& time_discretization) {
  KKTResidual kkt_residual(time_discretization.N_grids()+1, SplitKKTResidual(robot));
  for (int i=0; i<=time_discretization.N_grids(); ++i) {
    kkt_residual[i] = CreateSplitKKTResidual(robot);
    if (time_discretization.grid(i).switching_constraint) {
      const int impulse_index = time_discretization.grid(i).impulse_index + 1;
      kkt_residual[i].setSwitchingConstraintDimension(contact_sequence->impulseStatus(impulse_index).dimf());
      kkt_residual[i].P().setRandom();
    }
  }
  return kkt_residual;
}

} // namespace testhelper
} // namespace robotoc
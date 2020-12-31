#ifndef IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_DATA_HPP_ 
#define IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_DATA_HPP_

#include "Eigen/Core"
#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

class ImpulseDynamicsBackwardEulerData {
public:
  ImpulseDynamicsBackwardEulerData(const Robot& robot);

  ImpulseDynamicsBackwardEulerData();

  ~ImpulseDynamicsBackwardEulerData();

  ImpulseDynamicsBackwardEulerData(
      const ImpulseDynamicsBackwardEulerData&) = default;

  ImpulseDynamicsBackwardEulerData& operator=(
      const ImpulseDynamicsBackwardEulerData&) = default;
 
  ImpulseDynamicsBackwardEulerData(
      ImpulseDynamicsBackwardEulerData&&) noexcept = default;

  ImpulseDynamicsBackwardEulerData& operator=(
      ImpulseDynamicsBackwardEulerData&&) noexcept = default;

  void setImpulseStatus(const ImpulseStatus& Impulse_status);

  Eigen::Block<Eigen::MatrixXd> Qdvf();

  const Eigen::Block<const Eigen::MatrixXd> Qdvf() const;

  bool checkDimensions() const;

  Eigen::MatrixXd dImDdq;

  Eigen::MatrixXd dImDddv;

  Eigen::MatrixXd Minv;

  Eigen::MatrixXd Qdvq;

  Eigen::VectorXd ImD;

  Eigen::VectorXd Minv_ImD;

  Eigen::VectorXd ldv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd Qdvf_full_;
  int dimv_, dimf_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_backward_euler_data.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_DATA_HPP_ 
#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_DATA_HPP_ 
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_DATA_HPP_

#include "Eigen/Core"
#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

class ImpulseDynamicsForwardEulerData {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ImpulseDynamicsForwardEulerData(const Robot& robot);

  ImpulseDynamicsForwardEulerData();

  ~ImpulseDynamicsForwardEulerData();

  ImpulseDynamicsForwardEulerData(
      const ImpulseDynamicsForwardEulerData&) = default;

  ImpulseDynamicsForwardEulerData& operator=(
      const ImpulseDynamicsForwardEulerData&) = default;
 
  ImpulseDynamicsForwardEulerData(
      ImpulseDynamicsForwardEulerData&&) noexcept = default;

  ImpulseDynamicsForwardEulerData& operator=(
      ImpulseDynamicsForwardEulerData&&) noexcept = default;

  void setImpulseStatus(const ImpulseStatus& Impulse_status);

  Eigen::MatrixXd dImDdq;

  Eigen::MatrixXd dImDddv;

  Eigen::Block<Eigen::MatrixXd> dCdqv();

  const Eigen::Block<const Eigen::MatrixXd> dCdqv() const;

  Eigen::Block<Eigen::MatrixXd> dCdq();

  const Eigen::Block<const Eigen::MatrixXd> dCdq() const;

  Eigen::Block<Eigen::MatrixXd> dCdv();

  const Eigen::Block<const Eigen::MatrixXd> dCdv() const;

  Eigen::Block<Eigen::MatrixXd> dCddv();

  const Eigen::Block<const Eigen::MatrixXd> dCddv() const;

  Eigen::Block<Eigen::MatrixXd> MJtJinv();

  const Eigen::Block<const Eigen::MatrixXd> MJtJinv() const;

  Eigen::Block<Eigen::MatrixXd> MJtJinv_dImDCdqv();

  const Eigen::Block<const Eigen::MatrixXd> MJtJinv_dImDCdqv() const;

  Eigen::Block<Eigen::MatrixXd> Qdvfqv();

  const Eigen::Block<const Eigen::MatrixXd> Qdvfqv() const;

  Eigen::VectorBlock<Eigen::VectorXd> ImDC();

  const Eigen::VectorBlock<const Eigen::VectorXd> ImDC() const;

  Eigen::VectorBlock<Eigen::VectorXd> ImD();

  const Eigen::VectorBlock<const Eigen::VectorXd> ImD() const;

  Eigen::VectorBlock<Eigen::VectorXd> C();

  const Eigen::VectorBlock<const Eigen::VectorXd> C() const;

  Eigen::VectorBlock<Eigen::VectorXd> MJtJinv_ImDC();

  const Eigen::VectorBlock<const Eigen::VectorXd> MJtJinv_ImDC() const;

  Eigen::VectorBlock<Eigen::VectorXd> ldvf();

  const Eigen::VectorBlock<const Eigen::VectorXd> ldvf() const;

  Eigen::VectorBlock<Eigen::VectorXd> ldv();

  const Eigen::VectorBlock<const Eigen::VectorXd> ldv() const;

  Eigen::VectorBlock<Eigen::VectorXd> lf();

  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd dCdqv_full_, dCddv_full_, MJtJinv_full_, 
                  MJtJinv_dImDCdqv_full_, Qdvfqv_full_;
  Eigen::VectorXd ImDC_full_, MJtJinv_ImDC_full_, ldvf_full_;
  int dimv_, dimf_, dimvf_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_forward_euler_data.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_DATA_HPP_ 
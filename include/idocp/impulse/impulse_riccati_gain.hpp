#ifndef IDOCP_IMPULSE_RICCATI_GAIN_HPP_
#define IDOCP_IMPULSE_RICCATI_GAIN_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


namespace idocp {

class ImpulseRiccatiGain {
public:
  ImpulseRiccatiGain(const Robot& robot);

  ImpulseRiccatiGain();

  ~ImpulseRiccatiGain();

  ImpulseRiccatiGain(const ImpulseRiccatiGain&) = default;

  ImpulseRiccatiGain& operator=(const ImpulseRiccatiGain&) = default;
 
  ImpulseRiccatiGain(ImpulseRiccatiGain&&) noexcept = default;

  ImpulseRiccatiGain& operator=(ImpulseRiccatiGain&&) noexcept = default;

  void setContactStatus(const ContactStatus& contact_status);

  const Eigen::Block<const Eigen::MatrixXd> Kfq() const;

  const Eigen::Block<const Eigen::MatrixXd> Kfv() const;

  const Eigen::Block<const Eigen::MatrixXd> Kmuq() const;

  const Eigen::Block<const Eigen::MatrixXd> Kmuv() const;

  const Eigen::VectorBlock<const Eigen::VectorXd> kf() const;

  const Eigen::VectorBlock<const Eigen::VectorXd> kmu() const;

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Ginv, 
                           const Eigen::MatrixBase<MatrixType2>& Qfqv, 
                           const Eigen::MatrixBase<MatrixType3>& Cqv);

  template <typename MatrixType, typename VectorType1, typename VectorType2>
  void computeFeedforward(const Eigen::MatrixBase<MatrixType>& Ginv, 
                          const Eigen::MatrixBase<VectorType1>& lf, 
                          const Eigen::MatrixBase<VectorType2>& C);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_, dim_passive_, dimf_, dimc_;
  Eigen::MatrixXd K_;
  Eigen::VectorXd k_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_riccati_gain.hxx"

#endif // IDOCP_IMPULSE_RICCATI_GAIN_HPP_ 
#ifndef IDOCP_RICCATI_GAIN_HPP_
#define IDOCP_RICCATI_GAIN_HPP_

#include <assert.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"


namespace idocp {

class RiccatiGain {
public:
  RiccatiGain(const Robot& robot);

  RiccatiGain();

  ~RiccatiGain();

  RiccatiGain(const RiccatiGain&) = default;

  RiccatiGain& operator=(const RiccatiGain&) = default;
 
  RiccatiGain(RiccatiGain&&) noexcept = default;

  RiccatiGain& operator=(RiccatiGain&&) noexcept = default;

  void computeFeedbackGainAndFeedforward(const KKTMatrix& kkt_matrix, 
                                         const KKTResidual& kkt_residual);

  const Eigen::Block<const Eigen::MatrixXd> Kq() const;

  const Eigen::Block<const Eigen::MatrixXd> Kv() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::MatrixXd K, Ginv;
  Eigen::VectorXd k;

private:
  Eigen::LLT<Eigen::MatrixXd> llt_;
  int dimv_, dimu_;

};

} // namespace idocp 

#include "idocp/ocp/riccati_gain.hxx"

#endif // IDOCP_RICCATI_GAIN_HPP_
#ifndef IDOCP_RICCATI_GAIN_HPP_
#define IDOCP_RICCATI_GAIN_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class RiccatiGain {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RiccatiGain(const Robot& robot) 
    : Kaq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Kav(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Kfq_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
      Kfv_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
      Kmuq_(Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), 
                                  robot.dimv())),
      Kmuv_(Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), 
                                  robot.dimv())),
      ka(Eigen::VectorXd::Zero(robot.dimv())),
      kf_(Eigen::VectorXd::Zero(robot.max_dimf())),
      kmu_(Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf())),
      dimf_(robot.dimf()),
      dimc_(robot.dim_passive()+robot.dimf()) {
  }

  RiccatiGain() 
    : Kaq(),
      Kav(),
      Kfq_(),
      Kfv_(),
      Kmuq_(),
      Kmuv_(),
      ka(),
      kf_(),
      kmu_(),
      dimf_(0),
      dimc_(0) {
  }

  ~RiccatiGain() {
  }

  RiccatiGain(const RiccatiGain&) = default;

  RiccatiGain& operator=(const RiccatiGain&) = default;
 
  RiccatiGain(RiccatiGain&&) noexcept = default;

  RiccatiGain& operator=(RiccatiGain&&) noexcept = default;

  inline Eigen::Block<Eigen::MatrixXd> Kfq() {
    return Kfq_.topRows(dimf_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Kfv() {
    return Kfv_.topRows(dimf_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Kmuq() {
    return Kmuq_.topRows(dimc_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Kmuv() {
    return Kmuv_.topRows(dimc_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> kf() {
    return kf_.head(dimf_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> kmu() {
    return kmu_.head(dimc_);
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Kfq() const {
    return Kfq_.topRows(dimf_);
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Kfv() const {
    return Kfv_.topRows(dimf_);
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Kmuq() const {
    return Kmuq_.topRows(dimc_);
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Kmuv() const {
    return Kmuv_.topRows(dimc_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> kf() const {
    return kf_.head(dimf_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> kmu() const {
    return kmu_.head(dimc_);
  }

  Eigen::MatrixXd Kaq, Kav;
  Eigen::VectorXd ka;

private:
  int dimf_, dimc_;
  Eigen::MatrixXd Kfq_, Kfv_, Kmuq_, Kmuv_;
  Eigen::VectorXd kf_, kmu_;

};

} // namespace idocp 


#endif // IDOCP_RICCATI_GAIN_HPP_
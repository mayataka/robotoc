#ifndef IDOCP_KKT_RESIDUAL_HPP_
#define IDOCP_KKT_RESIDUAL_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class KKTResidual {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTResidual(const Robot& robot) 
    : lu(Eigen::VectorXd::Zero(robot.dimv())),
      u_res(Eigen::VectorXd::Zero(robot.dimv())),
      kkt_residual_(Eigen::VectorXd::Zero(
          5*robot.dimv()+robot.dim_passive()+2*robot.max_dimf())),
      dimv_(robot.dimv()), 
      dimx_(2*robot.dimv()), 
      dim_passive_(robot.dim_passive()), 
      max_dimf_(robot.max_dimf()), 
      dimf_(robot.dimf()), 
      max_dimc_(robot.dim_passive()+robot.max_dimf()), 
      dimc_(robot.dim_passive()+robot.dimf()),
      max_dimKKT_(5*robot.dimv()+robot.dim_passive()+2*robot.max_dimf()),
      dimKKT_(5*robot.dimv()+robot.dim_passive()+2*robot.dimf()) {
  }

  KKTResidual() 
    : lu(),
      u_res(),
      kkt_residual_(), 
      dimv_(0), 
      dimx_(0), 
      dim_passive_(0), 
      max_dimf_(0), 
      dimf_(0), 
      max_dimc_(0), 
      dimc_(0),
      max_dimKKT_(0),
      dimKKT_(0) {
  }

  ~KKTResidual() {
  }

  KKTResidual(const KKTResidual&) = default;

  KKTResidual& operator=(const KKTResidual&) = default;
 
  KKTResidual(KKTResidual&&) noexcept = default;

  KKTResidual& operator=(KKTResidual&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    dimf_ = robot.dimf();
    dimc_ = dim_passive_ + robot.dimf();
    dimKKT_ = 5*dimv_ + dim_passive_ + 2*robot.dimf();
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> KKT_residual() {
    return kkt_residual_.head(dimKKT_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> Fq() {
    return kkt_residual_.head(dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> Fv() {
    return kkt_residual_.segment(dimv_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> Fx() {
    return kkt_residual_.segment(0, dimx_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> C() {
    return kkt_residual_.segment(dimx_, dimc_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> la() {
    return kkt_residual_.segment(dimx_+dimc_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> lf() {
    return kkt_residual_.segment(dimx_+dimc_+dimv_, dimf_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> lq() {
    return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> lv() {
    return kkt_residual_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> lx() {
    return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> laf() {
    return kkt_residual_.segment(dimx_+dimc_, dimv_+dimf_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> KKT_residual() const {
    return kkt_residual_.head(dimKKT_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const {
    return kkt_residual_.head(dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const {
    return kkt_residual_.segment(dimv_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> Fx() const {
    return kkt_residual_.segment(0, dimx_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> C() const {
    return kkt_residual_.segment(dimx_, dimc_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> la() const {
    return kkt_residual_.segment(dimx_+dimc_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> lf() const {
    return kkt_residual_.segment(dimx_+dimc_+dimv_, dimf_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> lq() const {
    return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> lv() const {
    return kkt_residual_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> lx() const {
    return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
  }

  inline double squaredKKTErrorNorm(const double dtau) const {
    double error = kkt_residual_.head(dimKKT_).squaredNorm();
    error += lu.squaredNorm();
    error += dtau * dtau * u_res.squaredNorm();
    return error;
  }

  inline void setZeroMinimum() {
    lu.setZero();
    kkt_residual_.segment(dimx_+dimc_, 3*dimv_+dimf_).setZero();
  }

  inline void setZero() {
    lu.setZero();
    u_res.setZero();
    kkt_residual_.setZero();
  }

  inline int dimKKT() const {
    return dimKKT_;
  }

  inline int max_dimKKT() const {
    return max_dimKKT_;
  }

  inline int dimc() const {
    return dimc_;
  }

  inline int dimf() const {
    return dimf_;
  }

  Eigen::VectorXd lu, u_res;

private:
  Eigen::VectorXd kkt_residual_;
  int dimv_, dimx_, dim_passive_, max_dimf_, dimf_, max_dimc_, dimc_, dimKKT_, 
      max_dimKKT_;

};

} // namespace idocp 


#endif // IDOCP_KKT_RESIDUAL_HPP_
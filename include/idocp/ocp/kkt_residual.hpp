#ifndef IDOCP_KKT_RESIDUAL_HPP_
#define IDOCP_KKT_RESIDUAL_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints_interface.hpp"


namespace idocp {

class KKTResidual {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTResidual(const Robot& robot) 
    : kkt_matrix_(Eigen::MatrixXd::Zero(
          5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive(),
          5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive())),
      identity_matrix_(Eigen::MatrixXd::Identity(
          5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive(),
          5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive()))
      dimv_(robot.dimv()),
      dim_passive_(robot.dim_passive()),
      dimf_(robot.dimf()),
      dimc_(robot.dim_passive()+robot.dimf()) {
  }

  KKTResidual() 
    : kkt_matrix_(),
      identity_matrix_(),
      dimv_(0),
      dim_passive_(0),
      dimf_(0),
      dimc_(0) {
  }

  ~KKTResidual() {
  }

  KKTResidual(const KKTResidual&) = default;

  KKTResidual& operator=(const KKTResidual&) = default;
 
  KKTResidual(KKTResidual&&) noexcept = default;

  KKTResidual& operator=(KKTResidual&&) noexcept = default;

  void setContactStatus(const Robot& robot) {
    assert(robot.dim_passive() == dim_passive_);
    dimf_ = robot.dimf();
    dimc_ = robot.dim_passive() + robot.dimf();
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fqq() {
    return kkt_matrix_.block(Fq_begin_, Qq_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fqv() {
    return kkt_matrix_.block(Fq_begin_, Qv_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fvq() {
    return kkt_matrix_.block(Fv_begin_, Qq_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fvv() {
    return kkt_matrix_.block(Fv_begin_, Qv_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Cq() {
    return kkt_matrix_.block(C_begin_, Qq_begin_, dimc_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Cv {
    return kkt_matrix_.block(C_begin_, Qv_begin_, dimc_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Ca {
    return kkt_matrix_.block(C_begin_, Qa_begin_, dimc_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Cf {
    return kkt_matrix_.block(C_begin_+dimf_, Qf_begin_, dim_passive_, dimf_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qaa {
    return kkt_matrix_.block(Qa_begin_, Qa_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qaf {
    return kkt_matrix_.block(Qa_begin_, Qf_begin_, dimv_, dimf_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qaq {
    return kkt_matrix_.block(Qa_begin_, Qq_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qav {
    return kkt_matrix_.block(Qa_begin_, Qv_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qfa {
    return Qaf().transpose();
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qff {
    return kkt_matrix_.block(Qf_begin_, Qf_begin_, dimf_, dimf_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qfq {
    return kkt_matrix_.block(Qf_begin_, Qq_begin_, dimf_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qfq {
    return kkt_matrix_.block(Qf_begin_, Qq_begin_, dimf_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qfv {
    return kkt_matrix_.block(Qf_begin_, Qv_begin_, dimf_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqa {
    return Qaq().transpose();
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqf {
    return Qfq().transpose();
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqq {
    return kkt_matrix_.block(Qq_begin_, Qq_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqv {
    return kkt_matrix_.block(Qq_begin_, Qv_begin_, dimv_, dimv_);
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qva {
    return Qav().transpose();
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qvf {
    return Qfv().transpose();
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qvq {
    return Qqv().transpose();
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qvv {
    return kkt_matrix_.block(Qv_begin_, Qv_begin_, dimv_, dimv_);
  }

  void setZero() {
    kkt_residual_.setZero();
  }

  int size() const {
    return kkt_residual_.size();
  }

private:
  Eigen::VectorXd kkt_residual_;
  int dimv_, dim_passive_, dimf_, dimc_, size_, max_size_,
      Fq_begin_, Fv_begin_, C_begin_, Qa_begin_, Qf_begin_, Qq_begin_, Qv_begin_;

};

} // namespace idocp 


#endif // IDOCP_KKT_RESIDUAL_HPP_
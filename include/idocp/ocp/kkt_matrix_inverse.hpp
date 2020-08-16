#ifndef IDOCP_KKT_MATRIX_INVERSE_HPP_
#define IDOCP_KKT_MATRIX_INVERSE_HPP_

#include <assert.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class KKTMatrixInverse {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTMatrixInverse(const Robot& robot) 
    : kkt_composition_(robot),
      kkt_matrix_inverse_(Eigen::MatrixXd::Zero(kkt_composition_.max_dimKKT(), 
                                                kkt_composition_.max_dimKKT())) {
  }

  KKTMatrixInverse() 
    : kkt_composition_(),
      kkt_matrix_inverse_() {
  }

  ~KKTMatrixInverse() {
  }

  KKTMatrixInverse(const KKTMatrixInverse&) = default;

  KKTMatrixInverse& operator=(const KKTMatrixInverse&) = default;
 
  KKTMatrixInverse(KKTMatrixInverse&&) noexcept = default;

  KKTMatrixInverse& operator=(KKTMatrixInverse&&) noexcept = default;

  void setContactStatus(const Robot& robot) {
    kkt_composition_.setContactStatus(robot);
  }

  inline Eigen::Ref<Eigen::MatrixXd> KKT_matrix_inverse() {
    return kkt_matrix_inverse_.topLeftCorner(kkt_composition_.dimKKT(), 
                                             kkt_composition_.dimKKT());
  }

  inline Eigen::Ref<Eigen::MatrixXd> auxiliaryMatrix() {
    return kkt_matrix_inverse_.topLeftCorner(kkt_composition_.Fx_size(), 
                                             kkt_composition_.Fx_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> backwardCorrectionSerialCoeff() {
    return kkt_matrix_inverse_.block(kkt_composition_.Fx_begin(), 
                                     kkt_composition_.Qx_begin(), 
                                     kkt_composition_.Fx_size(), 
                                     kkt_composition_.Qx_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> backwardCorrectionParallelCoeff() {
    return kkt_matrix_inverse_.block(kkt_composition_.C_begin(), 
                                     kkt_composition_.Qx_begin(), 
                                     kkt_composition_.dimKKT()-kkt_composition_.Fx_size(), 
                                     kkt_composition_.Qx_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> forwardCorrectionSerialCoeff() {
    return kkt_matrix_inverse_.block(kkt_composition_.Qx_begin(), 
                                     kkt_composition_.Fx_begin(), 
                                     kkt_composition_.Qx_size(), 
                                     kkt_composition_.Fx_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> forwardCorrectionParallelCoeff() {
    return kkt_matrix_inverse_.topLeftCorner(kkt_composition_.dimKKT()-kkt_composition_.Fx_size(), 
                                             kkt_composition_.Qx_size());
  }

  inline Eigen::Ref<const Eigen::MatrixXd> KKT_matrix_inverse() const {
    return kkt_matrix_inverse_.topLeftCorner(kkt_composition_.dimKKT(), 
                                             kkt_composition_.dimKKT());
  }

  inline Eigen::Ref<const Eigen::MatrixXd> auxiliaryMatrix() const {
    return kkt_matrix_inverse_.topLeftCorner(kkt_composition_.Fx_size(), 
                                             kkt_composition_.Fx_size());
  }

  inline Eigen::Ref<const Eigen::MatrixXd> backwardCorrectionSerialCoeff() const {
    return kkt_matrix_inverse_.block(kkt_composition_.Fx_begin(), 
                                     kkt_composition_.Qx_begin(), 
                                     kkt_composition_.Fx_size(), 
                                     kkt_composition_.Qx_size());
  }

  inline Eigen::Ref<const Eigen::MatrixXd> backwardCorrectionParallelCoeff() const {
    return kkt_matrix_inverse_.block(kkt_composition_.C_begin(), 
                                     kkt_composition_.Qx_begin(), 
                                     kkt_composition_.dimKKT()-kkt_composition_.Fx_size(), 
                                     kkt_composition_.Qx_size());
  }

  inline Eigen::Ref<const Eigen::MatrixXd> forwardCorrectionSerialCoeff() const {
    return kkt_matrix_inverse_.block(kkt_composition_.Qx_begin(), 
                                     kkt_composition_.Fx_begin(), 
                                     kkt_composition_.Qx_size(), 
                                     kkt_composition_.Fx_size());
  }

  inline Eigen::Ref<const Eigen::MatrixXd> forwardCorrectionParallelCoeff() const {
    return kkt_matrix_inverse_.topLeftCorner(kkt_composition_.dimKKT()-kkt_composition_.Fx_size(), 
                                             kkt_composition_.Qx_size());
  }

  inline void setZero() {
    kkt_matrix_inverse_.setZero();
  }

  inline int dimKKT() const {
    return kkt_composition_.dimKKT();
  }

  inline int max_dimKKT() const {
    return kkt_composition_.max_dimKKT();
  }

private:
  KKTComposition kkt_composition_;
  Eigen::MatrixXd kkt_matrix_inverse_;

};

} // namespace idocp 


#endif // IDOCP_KKT_MATRIX_INVERSE_HPP_
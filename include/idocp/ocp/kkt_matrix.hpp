#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class KKTMatrix {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTMatrix(const Robot& robot) 
    : kkt_composition_(robot),
      kkt_matrix_(Eigen::MatrixXd::Zero(kkt_composition_.max_dimKKT(), 
                                        kkt_composition_.max_dimKKT())) {
  }

  KKTMatrix() 
    : kkt_composition_(),
      kkt_matrix_() {
  }

  ~KKTMatrix() {
  }

  KKTMatrix(const KKTMatrix&) = default;

  KKTMatrix& operator=(const KKTMatrix&) = default;
 
  KKTMatrix(KKTMatrix&&) noexcept = default;

  KKTMatrix& operator=(KKTMatrix&&) noexcept = default;

  void setContactStatus(const Robot& robot) {
    kkt_composition_.set(robot);
  }

  inline Eigen::Ref<Eigen::MatrixXd> KKT_matrix() {
    return kkt_matrix_.topLeftCorner(kkt_composition_.dimKKT(), 
                                     kkt_composition_.dimKKT());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fqq() {
    return kkt_matrix_.block(kkt_composition_.Fq_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Fq_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fqv() {
    return kkt_matrix_.block(kkt_composition_.Fq_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Fq_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fvq() {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Fvv() {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Cq() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Cv() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Ca() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Cf() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qaa() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qaf() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qaq() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qav() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qfa() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qff() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qfq() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qfv() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqa() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqf() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqq() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qqv() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qva() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qvf() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qvq() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qvv() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::MatrixXd> Qxx() {
    return kkt_matrix_.block(kkt_composition_.Qx_begin(), 
                             kkt_composition_.Qx_begin(), 
                             kkt_composition_.Qx_size(), 
                             kkt_composition_.Qx_size());
  }

  inline void symmetrize() {
    const int size = kkt_composition_.dimKKT();
    kkt_matrix_.topLeftCorner(size, size)
               .triangularView<Eigen::StrictlyLower>() 
        = kkt_matrix_.topLeftCorner(size, size).transpose()
                     .triangularView<Eigen::StrictlyLower>();
  }

  inline void invert(Eigen::Ref<Eigen::MatrixXd> kkt_matrix_inverse) {
    const int size = kkt_composition_.dimKKT();
    kkt_matrix_inverse 
        = kkt_matrix_.topLeftCorner(size, size)
                     .llt().solve(Eigen::MatrixXd::Identity(size, size));
  }

  inline void setZero() {
    kkt_matrix_.setZero();
  }

  inline int dimKKT() const {
    return kkt_composition_.dimKKT();
  }

  inline int max_dimKKT() const {
    return kkt_composition_.max_dimKKT();
  }

private:
  KKTComposition kkt_composition_;
  Eigen::MatrixXd kkt_matrix_;

};

} // namespace idocp 


#endif // IDOCP_KKT_MATRIX_HPP_
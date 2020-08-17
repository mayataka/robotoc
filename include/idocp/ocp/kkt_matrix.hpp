#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include <assert.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class KKTMatrix {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTMatrix(const Robot& robot) 
    : kkt_composition_(robot),
      kkt_matrix_(Eigen::MatrixXd::Zero(kkt_composition_.max_dimKKT(), 
                                        kkt_composition_.max_dimKKT())),
      luu(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      du_dq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      du_dv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      du_da(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
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
    kkt_composition_.setContactStatus(robot);
  }

  inline Eigen::Block<Eigen::MatrixXd> KKT_matrix() {
    return kkt_matrix_.topLeftCorner(kkt_composition_.dimKKT(), 
                                     kkt_composition_.dimKKT());
  }

  inline Eigen::Block<Eigen::MatrixXd> Fqq() {
    return kkt_matrix_.block(kkt_composition_.Fq_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Fq_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Fqv() {
    return kkt_matrix_.block(kkt_composition_.Fq_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Fq_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Fvq() {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Fvv() {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Fva() {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> CC() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.C_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.C_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Cq() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Cv() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Ca() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Cf() {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qaa() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qaf() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qaq() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qav() {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qfa() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qff() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qfq() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qfv() {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqa() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqf() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqq() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqv() {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qva() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qa_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qvf() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qf_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qvq() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qq_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qvv() {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qv_size());
  }

  inline Eigen::Block<Eigen::MatrixXd> Qxx() {
    return kkt_matrix_.block(kkt_composition_.Qx_begin(), 
                             kkt_composition_.Qx_begin(), 
                             kkt_composition_.Qx_size(), 
                             kkt_composition_.Qx_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> KKT_matrix() const {
    return kkt_matrix_.topLeftCorner(kkt_composition_.dimKKT(), 
                                     kkt_composition_.dimKKT());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Fqq() const {
    return kkt_matrix_.block(kkt_composition_.Fq_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Fq_size(), 
                             kkt_composition_.Qq_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Fqv() const {
    return kkt_matrix_.block(kkt_composition_.Fq_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Fq_size(), 
                             kkt_composition_.Qv_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Fvq() const {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qq_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Fvv() const {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qv_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Fva() const {
    return kkt_matrix_.block(kkt_composition_.Fv_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Fv_size(), 
                             kkt_composition_.Qa_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> CC() const {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.C_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.C_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Cq() const {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qq_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Cv() const {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qv_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Ca() const {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qa_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Cf() const {
    return kkt_matrix_.block(kkt_composition_.C_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.C_size(), 
                             kkt_composition_.Qf_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qaa() const {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qa_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qaf() const {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qf_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qaq() const {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qq_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qav() const {
    return kkt_matrix_.block(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qa_size(), 
                             kkt_composition_.Qv_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qfa() const {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qa_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qff() const {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qf_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qfq() const {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qq_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qfv() const {
    return kkt_matrix_.block(kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qf_size(), 
                             kkt_composition_.Qv_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qqa() const {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qa_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qqf() const {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qf_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qqq() const {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qq_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qqv() const {
    return kkt_matrix_.block(kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qq_size(), 
                             kkt_composition_.Qv_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qva() const {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qa_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qvf() const {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qf_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qf_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qvq() const {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qq_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qq_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qvv() const {
    return kkt_matrix_.block(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qv_size(), 
                             kkt_composition_.Qv_size());
  }

  inline const Eigen::Block<const Eigen::MatrixXd> Qxx() const {
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

  template <typename MatrixType>
  inline void invert(Eigen::MatrixBase<MatrixType>& kkt_matrix_inverse) {
    const int size = kkt_composition_.dimKKT();
    assert(kkt_matrix_inverse.rows() == size);
    assert(kkt_matrix_inverse.cols() == size);
    kkt_matrix_inverse 
        = kkt_matrix_.topLeftCorner(size, size)
                     .ldlt().solve(Eigen::MatrixXd::Identity(size, size));
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

  Eigen::MatrixXd luu, du_dq, du_dv, du_da, du_df;

private:
  KKTComposition kkt_composition_;
  Eigen::MatrixXd kkt_matrix_;

};

} // namespace idocp 


#endif // IDOCP_KKT_MATRIX_HPP_
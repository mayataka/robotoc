#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class KKTMatrix
/// @brief The KKT matrix of a time stage.
///
class KKTMatrix {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTMatrix(const Robot& robot);

  KKTMatrix();

  ~KKTMatrix();

  KKTMatrix(const KKTMatrix&) = default;

  KKTMatrix& operator=(const KKTMatrix&) = default;
 
  KKTMatrix(KKTMatrix&&) noexcept = default;

  KKTMatrix& operator=(KKTMatrix&&) noexcept = default;

  void setContactStatus(const Robot& robot);

  Eigen::Block<Eigen::MatrixXd> Ca();

  Eigen::Block<Eigen::MatrixXd> Cf();

  Eigen::Block<Eigen::MatrixXd> Cq();

  Eigen::Block<Eigen::MatrixXd> Cv();

  Eigen::Block<Eigen::MatrixXd> Caf();

  Eigen::Block<Eigen::MatrixXd> Cqv();

  Eigen::Block<Eigen::MatrixXd> Qaa();

  Eigen::Block<Eigen::MatrixXd> Qaf();

  Eigen::Block<Eigen::MatrixXd> Qaq();

  Eigen::Block<Eigen::MatrixXd> Qav();

  Eigen::Block<Eigen::MatrixXd> Qfa();

  Eigen::Block<Eigen::MatrixXd> Qff();

  Eigen::Block<Eigen::MatrixXd> Qfq();

  Eigen::Block<Eigen::MatrixXd> Qfv();

  Eigen::Block<Eigen::MatrixXd> Qqa();

  Eigen::Block<Eigen::MatrixXd> Qqf();

  Eigen::Block<Eigen::MatrixXd> Qqq();

  Eigen::Block<Eigen::MatrixXd> Qqv();

  Eigen::Block<Eigen::MatrixXd> Qva();

  Eigen::Block<Eigen::MatrixXd> Qvf();

  Eigen::Block<Eigen::MatrixXd> Qvq();

  Eigen::Block<Eigen::MatrixXd> Qvv();

  Eigen::Block<Eigen::MatrixXd> Qxx();

  Eigen::Block<Eigen::MatrixXd> Qafaf();

  Eigen::Block<Eigen::MatrixXd> Qafqv();

  Eigen::Block<Eigen::MatrixXd> costHessian();

  Eigen::Block<Eigen::MatrixXd> constraintsJacobian();

  void symmetrize();

  template <typename MatrixType>
  void invert(const double dtau, 
              const Eigen::MatrixBase<MatrixType>& kkt_matrix_inverse);

  void setZeroMinimum();

  void setZero();

  Eigen::MatrixXd Quu, Fqq, Fqq_prev;

private:
  Eigen::MatrixXd C_, Q_, Sc_, Sx_, FMinv_, C_H_inv_;
  bool has_floating_base_;
  int dimv_, dimx_, dimf_, dimc_, a_begin_, f_begin_, q_begin_, v_begin_, dimQ_;
  static constexpr int kDimFloatingBase = 6;

  template <typename MatrixType>
  void invertConstrainedHessian(
      const Eigen::MatrixBase<MatrixType>& hessian_inverse);

};

} // namespace idocp 

#include "idocp/ocp/kkt_matrix.hxx"

#endif // IDOCP_KKT_MATRIX_HPP_
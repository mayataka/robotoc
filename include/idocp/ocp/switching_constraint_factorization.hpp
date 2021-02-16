#ifndef IDOCP_SWITCHING_CONSTRAINT_FACTORIZATION_HPP_ 
#define IDOCP_SWITCHING_CONSTRAINT_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

class SwitchingConstraintFactorization {
public:
  SwitchingConstraintFactorization(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SwitchingConstraintFactorization();

  ///
  /// @brief Destructor. 
  ///
  ~SwitchingConstraintFactorization();

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingConstraintFactorization(
      const SwitchingConstraintFactorization&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SwitchingConstraintFactorization& operator=(
      const SwitchingConstraintFactorization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingConstraintFactorization(
      SwitchingConstraintFactorization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingConstraintFactorization& operator=(
      SwitchingConstraintFactorization&&) noexcept = default;

  ///
  /// @brief Set the dimension of the impulse. 
  /// @param[in] impulse_status Impulse status. 
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set the dimension of the impulse 0. 
  ///
  void setImpulseStatus();

  Eigen::Block<Eigen::MatrixXd> Phia();

  const Eigen::Block<const Eigen::MatrixXd> Phia() const;

  Eigen::Block<Eigen::MatrixXd> Phix();

  const Eigen::Block<const Eigen::MatrixXd> Phix() const;

  Eigen::Block<Eigen::MatrixXd> Phiq();

  const Eigen::Block<const Eigen::MatrixXd> Phiq() const;

  Eigen::Block<Eigen::MatrixXd> Phiv();

  const Eigen::Block<const Eigen::MatrixXd> Phiv() const;

  Eigen::Block<Eigen::MatrixXd> Phiu();

  const Eigen::Block<const Eigen::MatrixXd> Phiu() const;

  Eigen::Block<Eigen::MatrixXd> PhiuGinv();

  const Eigen::Block<const Eigen::MatrixXd> PhiuGinv() const;

  Eigen::Block<Eigen::MatrixXd> S();

  const Eigen::Block<const Eigen::MatrixXd> S() const;

  Eigen::Block<Eigen::MatrixXd> Sinv();

  const Eigen::Block<const Eigen::MatrixXd> Sinv() const;

  Eigen::Block<Eigen::MatrixXd> DGinv();

  const Eigen::Block<const Eigen::MatrixXd> DGinv() const; 

  Eigen::Block<Eigen::MatrixXd> SinvDGinv();

  const Eigen::Block<const Eigen::MatrixXd> SinvDGinv() const; 

  Eigen::Block<Eigen::MatrixXd> M();

  const Eigen::Block<const Eigen::MatrixXd> M() const;

  Eigen::VectorBlock<Eigen::VectorXd> m();

  const Eigen::VectorBlock<const Eigen::VectorXd> m() const;

  Eigen::VectorXd q;

  Eigen::VectorXd dq;

  Eigen::MatrixXd dintegrate_dq;

  Eigen::MatrixXd dintegrate_dv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd Phia_full_, Phix_full_, Phiu_full_, DGinv_full_, 
                  S_full_, Sinv_full_, SinvDGinv_full_, M_full_;
  Eigen::VectorXd m_full_;
  int dimv_, dimx_, dimu_, dimi_;

};

} // namespace idocp

#include "idocp/ocp/switching_constraint_factorization.hxx"

#endif // IDOCP_SWITCHING_CONSTRAINT_FACTORIZATION_HPP_
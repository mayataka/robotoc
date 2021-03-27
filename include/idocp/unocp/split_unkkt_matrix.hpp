#ifndef IDOCP_SPLIT_UNKKT_MATRIX_HPP_
#define IDOCP_SPLIT_UNKKT_MATRIX_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitUnKKTMatrix
/// @brief The KKT matrix split into a time stage.
///
class SplitUnKKTMatrix {
public:
  ///
  /// @brief Construct a KKT matrix.
  /// @param[in] robot Robot model. 
  ///
  SplitUnKKTMatrix(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitUnKKTMatrix();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnKKTMatrix();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitUnKKTMatrix(const SplitUnKKTMatrix&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitUnKKTMatrix& operator=(const SplitUnKKTMatrix&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  SplitUnKKTMatrix(SplitUnKKTMatrix&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnKKTMatrix& operator=(SplitUnKKTMatrix&&) noexcept = default;

  ///
  /// @brief Hessian of the Lagrangian with respect to the acceleration 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaa();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qaa().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaa() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the acceleration and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaq();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qaq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the acceleration and 
  /// generalized velocity. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qav();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qav().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qav() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and the 
  /// acceleratioin. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqa();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qqa().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqa() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqq();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// the acceleration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qva();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qva().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qva() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the acceleration 
  /// and state. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qax();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qax().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qax() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the state and 
  /// acceleration. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxa();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qxa().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxa() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxx();

  ///
  /// @brief const version of SplitUnKKTMatrix::Qxx().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxx() const;

  ///
  /// @brief Symmetrize the Hessian for matrix inversion. 
  ///
  void symmetrize();

  ///
  /// @brief Set the all components zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the condensed KKT condition.
  /// @return Dimension of the condensed KKT condition.
  ///
  int dimKKT() const;

  ///
  /// @brief Checks the equivalence of two SplitUnKKTMatrix.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitUnKKTMatrix& other) const;

  ///
  /// @brief Checks this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  Eigen::MatrixXd Q;

private:
  int dimv_, dimx_, dimKKT_;

};

} // namespace idocp 

#include "idocp/unocp/split_unkkt_matrix.hxx"

#endif // IDOCP_SPLIT_UNKKT_MATRIX_HPP_ 
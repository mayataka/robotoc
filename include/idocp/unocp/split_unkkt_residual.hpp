#ifndef IDOCP_SPLIT_UNKKT_RESIDUAL_HPP_ 
#define IDOCP_SPLIT_UNKKT_RESIDUAL_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class SplitUnKKTResidual;
using UnKKTResidual = std::vector<SplitUnKKTResidual>;

///
/// @brief KKT residual split into each time stage. 
///
class SplitUnKKTResidual {
public:
  ///
  /// @brief Construct a KKT residual.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  SplitUnKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  SplitUnKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnKKTResidual();

  ///
  /// @brief Use default copy constructor. 
  ///
  SplitUnKKTResidual(const SplitUnKKTResidual&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  SplitUnKKTResidual& operator=(const SplitUnKKTResidual&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  SplitUnKKTResidual(SplitUnKKTResidual&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  SplitUnKKTResidual& operator=(SplitUnKKTResidual&&) noexcept = default;

  ///
  /// @brief Residual with respect of the state equation.
  /// @return Reference to the residual of the state equation. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fq();

  ///
  /// @brief const version of SplitUnKKTResidual::Fq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const;

  ///
  /// @brief Residual with respect of the state equation.
  /// @return Reference to the residual of the state equation. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fv();

  ///
  /// @brief const version of SplitUnKKTResidual::Fv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const;

  ///
  /// @brief Residual with respect of the state equation.
  /// @return Reference to the residual of the state equation. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fx();

  ///
  /// @brief const version of SplitUnKKTResidual::Fx().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fx() const;

  ///
  /// @brief Residual with respect to the acceleration.
  /// @return Reference to the residual with respect to a. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> la();

  ///
  /// @brief const version of SplitUnKKTResidual::la().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> la() const;

  ///
  /// @brief Residual with respect to configuration q.
  /// @return Reference to the residual with respect to configuration q. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lq();

  ///
  /// @brief const version of SplitUnKKTResidual::lq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lq() const;

  ///
  /// @brief Residual with respect to generalized velocity v.
  /// @return Reference to the residual with respect to generalized velocity v.
  /// Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lv();

  ///
  /// @brief const version of SplitUnKKTResidual::lv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lv() const;

  ///
  /// @brief Residual with respect to state x.
  /// @return Reference to the residual with respect to state x. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lx();

  ///
  /// @brief const version of SplitUnKKTResidual::lv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lx() const;

  ///
  /// @brief Set the KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the KKT at the current contact status.
  /// @return Dimension of the KKT at the current contact status.
  ///
  int dimKKT() const;

  ///
  /// @brief Chech the equivalence of two SplitUnKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitUnKKTResidual& other) const;

  ///
  /// @brief Chech this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  /// @brief KKT residual.
  Eigen::VectorXd KKT_residual;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_, dimx_, dimKKT_;

};

} // namespace idocp 

#include "idocp/unocp/split_unkkt_residual.hxx"

#endif // IDOCP_SPLIT_UNKKT_RESIDUAL_HPP_ 
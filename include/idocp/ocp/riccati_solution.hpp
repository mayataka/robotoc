#ifndef IDOCP_RICCATI_SOLUTION_HPP_
#define IDOCP_RICCATI_SOLUTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class RiccatiSolution
/// @brief Riccati solution matrix and vector of a time stage.
///
class RiccatiSolution {
public:
  ///
  /// @brief Allocate Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  RiccatiSolution(const Robot& robot) 
    : Pqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pqv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      sq(Eigen::VectorXd::Zero(robot.dimv())),
      sv(Eigen::VectorXd::Zero(robot.dimv())) {
  }

  ///
  /// @brief Default constructor. 
  ///
  RiccatiSolution() 
    : Pqq(),
      Pqv(),
      Pvq(),
      Pvv(),
      sq(),
      sv() {
  }

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiSolution() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiSolution(const RiccatiSolution&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  RiccatiSolution& operator=(const RiccatiSolution&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiSolution(RiccatiSolution&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiSolution& operator=(RiccatiSolution&&) noexcept = default;

  ///
  /// @brief Top left corner of the Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pqq;

  ///
  /// @brief Top right corner of the Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pqv;

  ///
  /// @brief Bottom left corner of the Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pvq;

  ///
  /// @brief Bottom right corner of the Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pvv;

  ///
  /// @brief Head segment of the Riccati factorization vector. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorXd sq;

  ///
  /// @brief Tail segment of the Riccati factorization vector. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorXd sv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

} // namespace idocp 


#endif // IDOCP_RICCATI_SOLUTION_HPP_ 
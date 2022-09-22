#ifndef ROBOTOC_IMPULSE_DYNAMICS_DATA_HPP_ 
#define ROBOTOC_IMPULSE_DYNAMICS_DATA_HPP_

#include "Eigen/Core"
#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"


namespace robotoc {

///
/// @class ImpulseDynamicsData
/// @brief Data used in ImpulseDynamics.
///
class ImpulseDynamicsData {
public:
  ///
  /// @brief Constructs a data.
  /// @param[in] robot Robot model. 
  ///
  ImpulseDynamicsData(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseDynamicsData();

  ///
  /// @brief Default destructor. 
  ///
  ~ImpulseDynamicsData() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseDynamicsData(const ImpulseDynamicsData&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseDynamicsData& operator=(const ImpulseDynamicsData&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  ImpulseDynamicsData(ImpulseDynamicsData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseDynamicsData& operator=(ImpulseDynamicsData&&) noexcept = default;

  ///
  /// @brief Set the impulse status, i.e., set dimension of the impulses.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  Eigen::MatrixXd dImDddv;

  Eigen::Block<Eigen::MatrixXd> dImDCdqv();

  const Eigen::Block<const Eigen::MatrixXd> dImDCdqv() const;

  Eigen::Block<Eigen::MatrixXd> dImDCdq();

  const Eigen::Block<const Eigen::MatrixXd> dImDCdq() const;

  Eigen::Block<Eigen::MatrixXd> dImDdq();

  const Eigen::Block<const Eigen::MatrixXd> dImDdq() const;

  Eigen::Block<Eigen::MatrixXd> dCdq();

  const Eigen::Block<const Eigen::MatrixXd> dCdq() const;

  Eigen::Block<Eigen::MatrixXd> dCdv();

  const Eigen::Block<const Eigen::MatrixXd> dCdv() const;

  Eigen::Block<Eigen::MatrixXd> MJtJinv();

  const Eigen::Block<const Eigen::MatrixXd> MJtJinv() const;

  Eigen::Block<Eigen::MatrixXd> MJtJinv_dImDCdqv();

  const Eigen::Block<const Eigen::MatrixXd> MJtJinv_dImDCdqv() const;

  Eigen::Block<Eigen::MatrixXd> Qdvfqv();

  const Eigen::Block<const Eigen::MatrixXd> Qdvfqv() const;

  Eigen::VectorBlock<Eigen::VectorXd> ImDC();

  const Eigen::VectorBlock<const Eigen::VectorXd> ImDC() const;

  Eigen::VectorBlock<Eigen::VectorXd> ImD();

  const Eigen::VectorBlock<const Eigen::VectorXd> ImD() const;

  Eigen::VectorBlock<Eigen::VectorXd> C();

  const Eigen::VectorBlock<const Eigen::VectorXd> C() const;

  Eigen::VectorBlock<Eigen::VectorXd> MJtJinv_ImDC();

  const Eigen::VectorBlock<const Eigen::VectorXd> MJtJinv_ImDC() const;

  Eigen::VectorBlock<Eigen::VectorXd> ldvf();

  const Eigen::VectorBlock<const Eigen::VectorXd> ldvf() const;

  Eigen::VectorBlock<Eigen::VectorXd> ldv();

  const Eigen::VectorBlock<const Eigen::VectorXd> ldv() const;

  Eigen::VectorBlock<Eigen::VectorXd> lf();

  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

private:
  Eigen::MatrixXd dImDCdqv_full_, dCddv_full_, MJtJinv_full_, 
                  MJtJinv_dImDCdqv_full_, Qdvfqv_full_;
  Eigen::VectorXd ImDC_full_, MJtJinv_ImDC_full_, ldvf_full_;
  int dimv_, dimf_, dimvf_;

};

} // namespace robotoc 

#include "robotoc/dynamics/impulse_dynamics_data.hxx"

#endif // ROBOTOC_IMPULSE_DYNAMICS_DATA_HPP_ 
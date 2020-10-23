#ifndef IDOCP_IMPULSE_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_IMPULSE_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {

///
/// @class ImpulseCostFunctionComponentBase
/// @brief Base class of components of cost function.
///
class ImpulseCostFunctionComponentBase {
public:

  ///
  /// @brief Default constructor. 
  ///
  ImpulseCostFunctionComponentBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~ImpulseCostFunctionComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseCostFunctionComponentBase(const ImpulseCostFunctionComponentBase&) 
      = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseCostFunctionComponentBase& operator=(const ImpulseCostFunctionComponentBase&) 
      = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseCostFunctionComponentBase(ImpulseCostFunctionComponentBase&&) noexcept 
      = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseCostFunctionComponentBase& operator=(ImpulseCostFunctionComponentBase&&) noexcept 
      = default;

  ///
  /// @brief Computes and returns stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  virtual double l(Robot& robot, CostFunctionData& data, const double t, 
                   const ImpulseSplitSolution& s) const = 0;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  virtual void lq(Robot& robot, CostFunctionData& data, const double t, 
                  const ImpulseSplitSolution& s, 
                  ImpulseKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  virtual void lv(Robot& robot, CostFunctionData& data, const double t, 
                  const ImpulseSplitSolution& s, 
                  ImpulseKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the impulse change in the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  virtual void ldv(Robot& robot, CostFunctionData& data, const double t, 
                   const ImpulseSplitSolution& s, 
                   ImpulseKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] f Impulse force. Size must be ConstactStatus::dimf().
  /// @param[out] lf The KKT residual with respect to f. Size must be 
  /// ConstactStatus::dimf().
  ///
  virtual void lf(Robot& robot, CostFunctionData& data, const double t, 
                  const ImpulseSplitSolution& s, 
                  ImpulseKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void lqq(Robot& robot, CostFunctionData& data, const double t, 
                   const ImpulseSplitSolution& s, 
                   ImpulseKKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void lvv(Robot& robot, CostFunctionData& data, const double t, 
                   const ImpulseSplitSolution& s, 
                   ImpulseKKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the impulse change in the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void ldvdv(Robot& robot, CostFunctionData& data, const double t, 
                     const ImpulseSplitSolution& s, 
                     ImpulseKKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the Hessian of the stage cost with respect
  /// to the impulse velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] f Impulse force. Size must be ConstactStatus::dimf().
  /// @param[out] Qff The Hessian of the KKT residual with respect to f.  
  /// Size must be ConstactStatus::dimf() x ConstactStatus::dimf().
  ///
  virtual void lff(Robot& robot, CostFunctionData& data, const double t, 
                   const ImpulseSplitSolution& s, 
                   ImpulseKKTMatrix& kkt_matrix) const = 0;

};

} // namespace idocp

#endif // IDOCP_IMPULSE_COST_FUNCTION_COMPONENT_BASE_HPP_ 
#ifndef ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_
#define ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

#define DEFINE_DEFAULT_CLONE_COST_FUNCTION_COMPONENT(Derived) \
  std::shared_ptr<CostFunctionComponentBase> clone() const override { return std::make_shared<Derived>(*this); } 

///
/// @class CostFunctionComponentBase
/// @brief Base class of components of cost function.
///
class CostFunctionComponentBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  CostFunctionComponentBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~CostFunctionComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  CostFunctionComponentBase(const CostFunctionComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CostFunctionComponentBase& operator=(const CostFunctionComponentBase&) 
      = default;

  ///
  /// @brief Default move constructor. 
  ///
  CostFunctionComponentBase(CostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CostFunctionComponentBase& operator=(CostFunctionComponentBase&&) noexcept 
      = default;

  ///
  /// @brief Clones this to a shared ptr. 
  ///
  virtual std::shared_ptr<CostFunctionComponentBase> clone() const = 0;

  ///
  /// @brief Check if the cost function component requres kinematics
  /// (forward kinematics and its Jacobians) of robot model.
  /// @return true if the cost function component requres kinematics of 
  /// Robot model. false if not.
  ///
  virtual bool useKinematics() const = 0;

  ///
  /// @brief Computes the stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  virtual double evalStageCost(Robot& robot, const ContactStatus& contact_status,
                               CostFunctionData& data, const GridInfo& grid_info,
                               const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the first-order partial derivatives of the stage cost. 
  /// This function is always called just after evalStageCost().
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  ///
  virtual void evalStageCostDerivatives(Robot& robot, 
                                        const ContactStatus& contact_status,
                                        CostFunctionData& data, 
                                        const GridInfo& grid_info, 
                                        const SplitSolution& s, 
                                        SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessian, i.e., the second-order partial derivatives of 
  /// the stage cost.  This function is always called just after 
  /// evalStageCostDerivatives().
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_matrix Split KKT matrix. The Hessians are added to 
  /// this object.
  ///
  virtual void evalStageCostHessian(Robot& robot, 
                                    const ContactStatus& contact_status,
                                    CostFunctionData& data, 
                                    const GridInfo& grid_info, 
                                    const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the terminal cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @return Terminal cost.
  ///
  virtual double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                  const GridInfo& grid_info, 
                                  const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the first-order partial derivatives of the terminal cost. 
  /// This function is always called just after evalTerminalCost().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  ///
  virtual void evalTerminalCostDerivatives(Robot& robot, 
                                           CostFunctionData& data,  
                                           const GridInfo& grid_info, 
                                           const SplitSolution& s, 
                                           SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessian, i.e., the second-order partial derivatives of 
  /// the teminal cost. This function is always called just after 
  /// evalTerminalCostDerivatives().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_matrix Split KKT matrix. The Hessians are added to 
  /// this object.
  ///
  virtual void evalTerminalCostHessian(Robot& robot, 
                                       CostFunctionData& data, 
                                       const GridInfo& grid_info, 
                                       const SplitSolution& s, 
                                       SplitKKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the impulse cost. 
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @return Impulse cost.
  ///
  virtual double evalImpulseCost(Robot& robot, 
                                 const ImpulseStatus& impulse_status,
                                 CostFunctionData& data, 
                                 const GridInfo& grid_info, 
                                 const ImpulseSplitSolution& s) const = 0;

  ///
  /// @brief Computes the first-order partial derivatives of the impulse cost. 
  /// This function is always called just after evalImpulseCost().
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  ///
  virtual void evalImpulseCostDerivatives(Robot& robot, 
                                          const ImpulseStatus& impulse_status, 
                                          CostFunctionData& data,  
                                          const GridInfo& grid_info, 
                                          const ImpulseSplitSolution& s, 
                                          ImpulseSplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessian, i.e., the second-order partial derivatives of 
  /// the impulse cost. This function is always called just after 
  /// evalImpulseCostDerivatives().
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_matrix Impulse split KKT matrix. The Hessians are  
  /// added to this object.
  ///
  virtual void evalImpulseCostHessian(Robot& robot, 
                                      const ImpulseStatus& impulse_status, 
                                      CostFunctionData& data,  
                                      const GridInfo& grid_info, 
                                      const ImpulseSplitSolution& s, 
                                      ImpulseSplitKKTMatrix& kkt_matrix) const = 0; 

};

} // namespace robotoc

#endif // ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_
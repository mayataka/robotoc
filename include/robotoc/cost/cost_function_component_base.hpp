#ifndef ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_
#define ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_

#include <memory>
#include <stdexcept>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/ocp/grid_info.hpp"
#include "robotoc/cost/cost_function_data.hpp"


namespace robotoc {

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
  /// @brief Computes the impact cost. 
  /// @param[in] robot Robot model.
  /// @param[in] impact_status Impact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @return Impact cost.
  ///
  virtual double evalImpactCost(Robot& robot, 
                                 const ImpactStatus& impact_status,
                                 CostFunctionData& data, 
                                 const GridInfo& grid_info, 
                                 const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the first-order partial derivatives of the impact cost. 
  /// This function is always called just after evalImpactCost().
  /// @param[in] robot Robot model.
  /// @param[in] impact_status Impact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  ///
  virtual void evalImpactCostDerivatives(Robot& robot, 
                                          const ImpactStatus& impact_status, 
                                          CostFunctionData& data,  
                                          const GridInfo& grid_info, 
                                          const SplitSolution& s, 
                                          SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessian, i.e., the second-order partial derivatives of 
  /// the impact cost. This function is always called just after 
  /// evalImpactCostDerivatives().
  /// @param[in] robot Robot model.
  /// @param[in] impact_status Impact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_matrix Impact split KKT matrix. The Hessians are  
  /// added to this object.
  ///
  virtual void evalImpactCostHessian(Robot& robot, 
                                      const ImpactStatus& impact_status, 
                                      CostFunctionData& data,  
                                      const GridInfo& grid_info, 
                                      const SplitSolution& s, 
                                      SplitKKTMatrix& kkt_matrix) const = 0; 

  ///
  /// @brief Gets the shared ptr of this object as the specified type. If this 
  /// fails in dynamic casting, throws an exception.
  /// @tparam Derived The derived type.
  /// @return shared ptr of this object as the specified type. 
  ///
  template <typename Derived>
  std::shared_ptr<Derived> as() const {
    Derived* derived_ptr = dynamic_cast<Derived*>(this);
    if (derived_ptr == nullptr) {
        throw std::runtime_error("[CostFunctionComponentBase] runtime error: failed in down-casting!");
    }
    return std::shared_ptr<Derived>(derived_ptr);
  }

};

} // namespace robotoc

#endif // ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_
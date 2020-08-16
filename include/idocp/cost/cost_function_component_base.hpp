#ifndef IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class CostFunctionComponentBase {
public:
  CostFunctionComponentBase() {}

  virtual ~CostFunctionComponentBase() {}

  // Use default copy constructor.
  CostFunctionComponentBase(const CostFunctionComponentBase&) = default;

  // Use default copy coperator.
  CostFunctionComponentBase& operator=(const CostFunctionComponentBase&) 
      = default;

  // Use default move constructor.
  CostFunctionComponentBase(CostFunctionComponentBase&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionComponentBase& operator=(CostFunctionComponentBase&&) noexcept 
      = default;

  virtual double l(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s) const = 0;

  virtual double phi(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s) const = 0;

  template <typename VectorType>
  virtual void lq(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  const Eigen::MatrixBase<VectorType>& lq) const = 0;

  template <typename VectorType>
  virtual void lv(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  const Eigen::MatrixBase<VectorType>& lv) const = 0;

  template <typename VectorType>
  virtual void la(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  const Eigen::MatrixBase<VectorType>& la) const = 0;

  template <typename VectorType>
  virtual void lf(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  const Eigen::MatrixBase<VectorType>& lf) const = 0;

  template <typename VectorType>
  virtual void lu(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  const Eigen::MatrixBase<VectorType>& lu) const = 0;

  template <typename MatrixType>
  virtual void lqq(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& lqq) const = 0;

  template <typename MatrixType>
  virtual void lvv(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& lvv) const = 0;

  template <typename MatrixType>
  virtual void laa(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& laa) const = 0;

  template <typename MatrixType>
  virtual void lff(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& lff) const = 0;

  template <typename MatrixType>
  virtual void luu(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& luu) const = 0;

  template <typename MatrixType>
  virtual void augment_lqq(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const SplitSolution& s,
                           const Eigen::MatrixBase<MatrixType>& lqq) const = 0;

  template <typename MatrixType>
  virtual void augment_lvv(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const SplitSolution& s,
                           const Eigen::MatrixBase<MatrixType>& lvv) const = 0;

  template <typename MatrixType>
  virtual void augment_laa(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const SplitSolution& s,
                           const Eigen::MatrixBase<MatrixType>& laa) const = 0;

  template <typename MatrixType>
  virtual void augment_lff(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const SplitSolution& s,
                           const Eigen::MatrixBase<MatrixType>& lff) const = 0;

  template <typename MatrixType>
  virtual void augment_luu(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const SplitSolution& s,
                           const Eigen::MatrixBase<MatrixType>& luu) const = 0;

  template <typename VectorType>
  virtual void phiq(const Robot& robot, CostFunctionData& data, const double t, 
                    const SplitSolution& s,
                    const Eigen::MatrixBase<VectorType>& phiq) const = 0;

  template <typename VectorType>
  virtual void phiv(const Robot& robot, CostFunctionData& data, const double t, 
                    const SplitSolution& s,
                    const Eigen::MatrixBase<VectorType>& phiv) const = 0;

  template <typename MatrixType>
  virtual void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s,
                     const Eigen::MatrixBase<MatrixType>& phiqq) const = 0;

  template <typename MatrixType>
  virtual void phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s,
                     const Eigen::MatrixBase<MatrixType>& phivv) const = 0;

  template <typename MatrixType>
  virtual void augment_phiqq(
      const Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, 
      const Eigen::MatrixBase<MatrixType>& phiqq) const = 0; 

  template <typename MatrixType>
  virtual void augment_phivv(
      const Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, 
      const Eigen::MatrixBase<MatrixType>& phivv) const = 0; 

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_
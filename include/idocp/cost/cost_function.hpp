#ifndef IDOCP_COST_FUNCTION_HPP_
#define IDOCP_COST_FUNCTION_HPP_

#include <vector>
#include <memory>
#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class CostFunction {
public:
  CostFunction() 
    : costs_() {
  }

  ~CostFunction() {}

  // Use default copy constructor.
  CostFunction(const CostFunction&) = default;

  // Use default copy coperator.
  CostFunction& operator=(const CostFunction&) = default;

  // Use default move constructor.
  CostFunction(CostFunction&&) noexcept = default;

  // Use default move assign coperator.
  CostFunction& operator=(CostFunction&&) noexcept = default;

  void push_back(const std::shared_ptr<CostFunctionComponentBase>& cost) {
    costs_.push_back(cost);
  }

  void clear() {
    costs_.clear();
  }

  bool isEmpty() {
    return costs_.empty();
  }

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const {
    assert(dtau > 0);
    double l = 0;
    for (const auto cost : costs_) {
      l += cost->l(robot, data, t, dtau, s);
    }
    return l;
  }

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const {
    double phi = 0;
    for (const auto cost : costs_) {
      phi += cost->phi(robot, data, t, s);
    }
    return phi;
  }

  template <typename VectorType>
  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          const Eigen::MatrixBase<VectorType>& lq) const {
    assert(dtau > 0);
    assert(lq.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lq(robot, data, t, dtau, s, 
               const_cast<Eigen::MatrixBase<VectorType>&>(lq));
    }
  }

  template <typename VectorType>
  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          const Eigen::MatrixBase<VectorType>& lv) const {
    assert(dtau > 0);
    assert(lv.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lv(robot, data, t, dtau, s, 
               const_cast<Eigen::MatrixBase<VectorType>&>(lv));
    }
  }

  template <typename VectorType>
  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          const Eigen::MatrixBase<VectorType>& la) const {
    assert(dtau > 0);
    assert(la.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->la(robot, data, t, dtau, s, 
               const_cast<Eigen::MatrixBase<VectorType>&>(la));
    }
  }

  template <typename VectorType>
  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          const Eigen::MatrixBase<VectorType>& lf) const {
    assert(dtau > 0);
    assert(lf.size() <= robot.max_dimf());
    assert(lf.size() >= robot.dimf());
    for (const auto cost : costs_) {
      cost->lf(robot, data, t, dtau, s, 
               const_cast<Eigen::MatrixBase<VectorType>&>(lf));
    }
  }

  template <typename VectorType>
  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          const Eigen::MatrixBase<VectorType>& lu) const {
    assert(dtau > 0);
    assert(lu.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lu(robot, data, t, dtau, s, 
               const_cast<Eigen::MatrixBase<VectorType>&>(lu));
    }
  }

  template <typename MatrixType>
  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lqq) const {
    assert(dtau > 0);
    assert(lqq.rows() == robot.dimv());
    assert(lqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lqq(robot, data, t, dtau, s, 
                const_cast<Eigen::MatrixBase<MatrixType>&>(lqq));
    }
  }

  template <typename MatrixType>
  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lvv) const {
    assert(dtau > 0);
    assert(lvv.rows() == robot.dimv());
    assert(lvv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lvv(robot, data, t, dtau, s,
                const_cast<Eigen::MatrixBase<MatrixType>&>(lvv));
    }
  }

  template <typename MatrixType>
  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& laa) const {
    assert(dtau > 0);
    assert(laa.rows() == robot.dimv());
    assert(laa.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->laa(robot, data, t, dtau, s, 
                const_cast<Eigen::MatrixBase<MatrixType>&>(laa));
    }
  }

  template <typename MatrixType>
  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lff) const {
    assert(dtau > 0);
    assert(lff.rows() >= robot.dimf());
    assert(lff.rows() <= robot.max_dimf());
    assert(lff.cols() >= robot.dimf());
    assert(lff.cols() <= robot.max_dimf());
    for (const auto cost : costs_) {
      cost->lff(robot, data, t, dtau, s,
                const_cast<Eigen::MatrixBase<MatrixType>&>(lff));
    }
  }

  template <typename MatrixType>
  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& luu) const {
    assert(dtau > 0);
    assert(luu.rows() == robot.dimv());
    assert(luu.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->luu(robot, data, t, dtau, s,
                const_cast<Eigen::MatrixBase<MatrixType>&>(luu));
    }
  }

  template <typename MatrixType>
  void augment_lqq(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& lqq) const {
    assert(dtau > 0);
    assert(lqq.rows() == robot.dimv());
    assert(lqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_lqq(robot, data, t, dtau, s,
                        const_cast<Eigen::MatrixBase<MatrixType>&>(lqq));
    }
  }

  template <typename MatrixType>
  void augment_lvv(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& lvv) const {
    assert(dtau > 0);
    assert(lvv.rows() == robot.dimv());
    assert(lvv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_lvv(robot, data, t, dtau, s, 
                        const_cast<Eigen::MatrixBase<MatrixType>&>(lvv));
    }
  }

  template <typename MatrixType>
  void augment_laa(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& laa) const {
    assert(dtau > 0);
    assert(laa.rows() == robot.dimv());
    assert(laa.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_laa(robot, data, t, dtau, s, 
                        const_cast<Eigen::MatrixBase<MatrixType>&>(laa));
    }
  }

  template <typename MatrixType>
  void augment_lff(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& lff) const {
    assert(dtau > 0);
    assert(lff.rows() >= robot.dimf());
    assert(lff.rows() <= robot.max_dimf());
    assert(lff.cols() >= robot.dimf());
    assert(lff.cols() <= robot.max_dimf());
    for (const auto cost : costs_) {
      cost->augment_lff(robot, data, t, dtau, s, 
                        const_cast<Eigen::MatrixBase<MatrixType>&>(lff));
    }
  }

  template <typename MatrixType>
  void augment_luu(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& luu) const {
    assert(dtau > 0);
    assert(luu.rows() == robot.dimv());
    assert(luu.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_luu(robot, data, t, dtau, s, 
                        const_cast<Eigen::MatrixBase<MatrixType>&>(luu));
    }
  }

  template <typename VectorType>
  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, 
            const Eigen::MatrixBase<VectorType>& phiq) const {
    assert(phiq.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phiq(robot, data, t, s, 
                 const_cast<Eigen::MatrixBase<VectorType>&>(phiq));
    }
  }

  template <typename VectorType>
  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, 
            const Eigen::MatrixBase<VectorType>& phiv) const {
    assert(phiv.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phiv(robot, data, t, s, 
                 const_cast<Eigen::MatrixBase<VectorType>&>(phiv));
    }
  }

  template <typename MatrixType>
  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             const Eigen::MatrixBase<MatrixType>& phiqq) const {
    assert(phiqq.rows() == robot.dimv());
    assert(phiqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phiqq(robot, data, t, s, 
                  const_cast<Eigen::MatrixBase<MatrixType>&>(phiqq));
    }
  }

  template <typename MatrixType>
  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             const Eigen::MatrixBase<MatrixType>& phivv) const {
    assert(phivv.rows() == robot.dimv());
    assert(phivv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phivv(robot, data, t, s, 
                  const_cast<Eigen::MatrixBase<MatrixType>&>(phivv));
    }
  }

  template <typename MatrixType>
  void augment_phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const double dtau, const SplitSolution& s, 
                     const Eigen::MatrixBase<MatrixType>& phiqq) const {
    assert(phiqq.rows() == robot.dimv());
    assert(phiqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_phiqq(robot, data, t, s, 
                          const_cast<Eigen::MatrixBase<MatrixType>&>(phiqq));
    }
  }

  template <typename MatrixType>
  void augment_phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const double dtau, const SplitSolution& s, 
                     const Eigen::MatrixBase<MatrixType>& phivv) const {
    assert(phivv.rows() == robot.dimv());
    assert(phivv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_phivv(robot, data, t, s, 
                          const_cast<Eigen::MatrixBase<MatrixType>&>(phivv));
    }
  }

private:
  std::vector<std::shared_ptr<CostFunctionComponentBase>> costs_;

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_HPP_
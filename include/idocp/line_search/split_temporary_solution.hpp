#ifndef IDOCP_SPLIT_TEMPORARY_SOLUTION_HPP_
#define IDOCP_SPLIT_TEMPORARY_SOLUTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

///
/// @class SplitTemporarySolution
/// @brief Split temporary solution used in line search. 
///
class SplitTemporarySolution {
public:
  ///
  /// @brief Construct SplitTemporarySolution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  SplitTemporarySolution(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitTemporarySolution();

  ///
  /// @brief Destructor. 
  ///
  ~SplitTemporarySolution();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitTemporarySolution(const SplitTemporarySolution&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitTemporarySolution& operator=(const SplitTemporarySolution&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitTemporarySolution(SplitTemporarySolution&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitTemporarySolution& operator=(SplitTemporarySolution&&) noexcept = default;

  void setTemporarySolution(const Robot& robot, 
                            const ContactStatus& contact_status,
                            const double step_size, const SplitSolution& s, 
                            const SplitDirection& d, 
                            const SplitSolution& s_next,
                            const SplitDirection& d_next);

  void setTemporarySolution(const Robot& robot, const double step_size, 
                            const SplitSolution& s, const SplitDirection& d);

  const SplitSolution& splitSolution() const;

  const Eigen::VectorXd& q_next() const;

  const Eigen::VectorXd& v_next() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  SplitSolution s_tmp_;
  Eigen::VectorXd q_next_tmp_, v_next_tmp_;
};

} // namespace idocp

#include "idocp/ocp/split_temporary_solution.hxx"

#endif // IDOCP_SPLIT_TEMPORARY_SOLUTION_HPP_ 
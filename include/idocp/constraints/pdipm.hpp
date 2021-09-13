#ifndef IDOCP_CONSTRAINTS_PDIPM_HPP_
#define IDOCP_CONSTRAINTS_PDIPM_HPP_

#include "Eigen/Core"

#include "idocp/constraints/constraint_component_data.hpp"


namespace idocp {
namespace pdipm {

///
/// @brief Sets the slack and dual variables positive.
/// @param[in] barrier Barrier parameter. Must be positive. 
/// @param[in, out] data Constraint component data.
///
void setSlackAndDualPositive(const double barrier, 
                             ConstraintComponentData& data);

///
/// @brief Computes the residual in the complementarity slackness between  
/// the slack and dual variables.
/// @param[in] barrier Barrier parameter. Must be positive. 
/// @param[in, out] data Constraint component data.
///
void computeComplementarySlackness(const double barrier, 
                                   ConstraintComponentData& data);

///
/// @brief Computes the residual in the complementarity slackness between  
/// the slack and dual variables.
/// @param[in] barrier Barrier parameter. Must be positive. 
/// @param[in, out] data Constraint data.
/// @param[in] start Start position of the segment.
/// @param[in] size Size of the segment.
///
void computeComplementarySlackness(const double barrier, 
                                   ConstraintComponentData& data,
                                   const int start, const int size);

///
/// @brief Computes the residual in the complementarity slackness between  
/// the slack and dual variables.
/// @param[in] barrier Barrier parameter. Must be positive. 
/// @param[in, out] data Constraint data.
/// @param[in] start Start position of the segment.
/// @tparam Size Size of the segment.
///
template <int Size>
void computeComplementarySlackness(const double barrier, 
                                   ConstraintComponentData& data,
                                   const int start);

///
/// @brief Computes the residual in the complementarity slackness between  
/// the slack and dual variables.
/// @param[in] barrier Barrier parameter. Must be positive. 
/// @param[in] slack An element of the slack variable.
/// @param[in] dual An element of the dual variable.
/// @return The complementarity slackness between the slack and dual variables.
///
double computeComplementarySlackness(const double barrier, const double slack, 
                                     const double dual);

///
/// @brief Computes the coefficient of the condensing.
/// @param[in, out] data Constraint component data.
///
void computeCondensingCoeffcient(ConstraintComponentData& data);

///
/// @brief Computes the coefficient of the condensing.
/// @param[in, out] data Constraint data.
/// @param[in] start Start position of the segment.
/// @param[in] size Size of the segment.
///
void computeCondensingCoeffcient(ConstraintComponentData& data,
                                 const int start, const int size);

///
/// @brief Computes the coefficient of the condensing.
/// @param[in, out] data Constraint data.
/// @param[in] start Start position of the segment.
/// @tparam Size Size of the segment.
///
template <int Size>
void computeCondensingCoeffcient(ConstraintComponentData& data,
                                 const int start);

///
/// @brief Computes the residual in the complementarity slackness between  
/// the slack and dual variables.
/// @param[in] slack An element of the slack variable.
/// @param[in] dual An element of the dual variable.
/// @param[in] residual An element of the primal residual.
/// @param[in] cmpl An element of the complementarity slackness.
/// @return Coefficient of the condensing. 
///
double computeCondensingCoeffcient(const double slack, const double dual,
                                   const double residual, const double cmpl);

///
/// @brief Applies the fraction-to-boundary-rule to the directions of the slack 
/// variables.
/// @param[in] fraction_rate Must be larger than 0 and smaller than 1. Should be 
/// between 0.9 and 0.995.
/// @param[in] data Constraint component data.
/// @return Fraction-to-boundary of the direction of the slack variables. 
///
double fractionToBoundarySlack(const double fraction_rate, 
                               const ConstraintComponentData& data);

///
/// @brief Applies the fraction-to-boundary-rule to the directions of the dual
/// variables.
/// @param[in] fraction_rate Must be larger than 0 and smaller than 1. Should be 
/// between 0.9 and 0.995.
/// @param[in] data Constraint component data.
/// @return Fraction-to-boundary of the direction of the dual variables. 
///
double fractionToBoundaryDual(const double fraction_rate, 
                              const ConstraintComponentData& data);

///
/// @brief Applies the fraction-to-boundary-rule.
/// @param[in] dim Dimension of vec and dvec. 
/// @param[in] fraction_rate Must be larger than 0 and smaller than 1. Should be 
/// between 0.9 and 0.995.
/// @param[in] vec A vector. 
/// @param[in] dvec A direction vector of vec. 
/// @return Fraction-to-boundary of dvec.
///
double fractionToBoundary(const int dim, const double fraction_rate, 
                          const Eigen::VectorXd& vec,
                          const Eigen::VectorXd& dvec);

///
/// @brief Applies the fraction-to-boundary-rule.
/// @param[in] fraction_rate Must be larger than 0 and smaller than 1. Should be 
/// between 0.9 and 0.995.
/// @param[in] var A variable. 
/// @param[in] dvar A direction the variable. 
/// @return Fraction-to-boundary of dvec.
///
double fractionToBoundary(const double fraction_rate, 
                          const double var, const double dvar);

///
/// @brief Computes the direction of the dual variable from slack, primal 
/// residual, complementary slackness, and the direction of the slack.
/// @param[in, out] data Constraint component data.
///
void computeDualDirection(ConstraintComponentData& data);

///
/// @brief Computes the direction of the dual variable from slack, primal 
/// residual, complementary slackness, and the direction of the slack.
/// @param[in, out] data Constraint data.
/// @param[in] start Start position of the segment.
/// @param[in] size Size of the segment.
///
void computeDualDirection(ConstraintComponentData& data, 
                          const int start, const int size);

///
/// @brief Computes the direction of the dual variable from slack, primal 
/// residual, complementary slackness, and the direction of the slack.
/// @param[in, out] data Constraint data.
/// @param[in] start Start position of the segment.
///
template <int Size>
void computeDualDirection(ConstraintComponentData& data, const int start);

///
/// @brief Computes the direction of the dual variable from slack, primal 
/// residual, complementary slackness, and the direction of the slack.
/// @param[in] slack The slack variable.
/// @param[in] dual The dual variable.
/// @param[in] dslack The direction of the slack variable.
/// @param[in] cmpl The complementary slackness.
/// @return The direction of the dual variable.
///
double computeDualDirection(const double slack, const double dual,
                            const double dslack, const double cmpl);

///
/// @brief Computes the log barrier function.
/// @param[in] barrier Barrier parameter. Must be positive. 
/// @param[in] vec Argument of the log function. All the components must be 
/// positive.
/// @return log barrier function.
///
template <typename VectorType>
double logBarrier(const double barrier, 
                   const Eigen::MatrixBase<VectorType>& vec);

} // namespace pdipm
} // namespace idocp

#include "idocp/constraints/pdipm.hxx"

#endif // IDOCP_CONSTRAINTS_PDIPM_HPP_
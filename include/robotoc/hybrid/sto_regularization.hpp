#ifndef ROBOTOC_STO_REGULARIZATION_HPP_
#define ROBOTOC_STO_REGULARIZATION_HPP_

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"


namespace robotoc {

/// 
/// @enum STORegularizationType
/// @brief Type of the regularization for switching time optimization (STO) 
/// problem.
///
enum class STORegularizationType {
  Const,
  Abs,
  Square,
  Quad,
  None
};

///
/// @class STORegularization
/// @brief Regularization for switching time optimization (STO) problem.
///
class STORegularization {
public:
  ///
  /// @brief Constructs regularization.
  /// @param[in] reg_type regularization type. 
  /// @param[in] w weight parameter (a.k.a. scaling parameter) of the 
  /// regularization. Must be non-negative.
  ///
  STORegularization(const STORegularizationType& reg_type, const double w);

  ///
  /// @brief Default constructor. 
  ///
  STORegularization();

  ///
  /// @brief Destructor. 
  ///
  ~STORegularization();

  ///
  /// @brief Default copy constructor. 
  ///
  STORegularization(const STORegularization&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  STORegularization& operator=(const STORegularization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  STORegularization(STORegularization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  STORegularization& operator=(STORegularization&&) noexcept = default;

  ///
  /// @brief Sets regularization method.
  /// @param[in] reg_type regularization type. 
  /// @param[in] w weight parameter (a.k.a. scaling parameter) of the 
  /// regularization. Must be non-negative.
  ///
  void setRegularization(const STORegularizationType& reg_type,
                         const double w);

  ///
  /// @return Applies the regularization on the KKT matrix.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_error Squared norm of the KKT residual. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @note KKT error must be the squared norm of the KKT residual.
  ///
  void applyRegularization(const OCP& ocp, const double kkt_error, 
                           KKTMatrix& kkt_matrix) const; 

  ///
  /// @return true if the regularization is valid. faise if not.
  ///
  bool isRegularizationValid() const;

  ///
  /// @return Get the regularization from the KKt error.
  /// @param[in] kkt_error KKT error. 
  ///
  double getRegularization(const double kkt_error) const;

  ///
  /// @return Default STO regularization. STORegularizationType is 
  /// STORegularizationType::Quad and weight parameter is 0.1.
  ///
  static STORegularization defaultSTORegularization();

private:
  STORegularizationType reg_type_;
  double w_;
  bool is_reg_valid_;

};

} // namespace robotoc

#include "robotoc/hybrid/sto_regularization.hxx"

#endif // ROBOTOC_STO_REGULARIZATION_HPP_ 
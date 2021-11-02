#ifndef robotoc_sto_regularization_hpp_
#define robotoc_sto_regularization_hpp_

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"


namespace robotoc {

/// 
/// @enum STORegularizationType
/// @brief type of the regularization for switching time optimization (sto).
///
enum class STORegularizationType {
  Const,
  Abs,
  Quad,
  Exp,
  Exp2,
  None
};

///
/// @class STORegularization
/// @brief regularization for switching time optimization (sto) problems.
///
class STORegularization {
public:
  ///
  /// @brief constructs regularization.
  /// @param[in] reg_type regularization type. 
  /// @param[in] w weight parameter (a.k.a. scaling parameter) of the 
  /// regularization.
  ///
  STORegularization(const STORegularizationType& reg_type, const double w);

  ///
  /// @brief default constructor. 
  ///
  STORegularization();

  ///
  /// @brief destructor. 
  ///
  ~STORegularization();

  ///
  /// @brief default copy constructor. 
  ///
  STORegularization(const STORegularization&) = default;

  ///
  /// @brief default copy assign operator. 
  ///
  STORegularization& operator=(const STORegularization&) = default;

  ///
  /// @brief default move constructor. 
  ///
  STORegularization(STORegularization&&) noexcept = default;

  ///
  /// @brief default move assign operator. 
  ///
  STORegularization& operator=(STORegularization&&) noexcept = default;

  ///
  /// @brief sets regularization method.
  /// @param[in] reg_type regularization type. 
  /// @param[in] w weight parameter (a.k.a. scaling parameter) of the 
  /// regularization.
  ///
  void setRegularization(const STORegularizationType& reg_type,
                         const double w);

  ///
  /// @return Applies the regularization on the KKT matrix.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_error KKT error. 
  /// @param[in, out] kkt_matrix KKT matrix. 
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

private:
  STORegularizationType reg_type_;
  double w_;
  bool is_reg_valid_;

};

} // namespace robotoc

#include "robotoc/hybrid/sto_regularization.hxx"

#endif // ROBOTOC_STO_REGULARIZATION_HPP_ 
#ifndef ROBOTOC_STO_REGULARIZATION_HXX_
#define ROBOTOC_STO_REGULARIZATION_HXX_

#include "robotoc/hybrid/sto_regularization.hpp"

#include <stdexcept>
#include <iostream>
#include <cmath>


namespace robotoc {

inline STORegularization::STORegularization(
    const STORegularizationType& reg_type, const double w)
  : reg_type_(reg_type),
    w_(w),
    is_reg_valid_(false) {
  try {
    if (w < 0) {
      throw std::out_of_range(
          "invalid value: w must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline STORegularization::STORegularization()
  : reg_type_(STORegularizationType::None),
    w_(0),
    is_reg_valid_(false) {
}


inline STORegularization::~STORegularization() {
}


inline void STORegularization::setRegularization(
    const STORegularizationType& reg_type, const double w) {
  try {
    if (w < 0) {
      throw std::out_of_range(
          "invalid value: w must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  reg_type_ = reg_type;
  w_ = w;
  if (reg_type == STORegularizationType::None) {
    is_reg_valid_ = false;
  }
  else {
    is_reg_valid_ = true;
  }
}


inline void STORegularization::applyRegularization(const OCP& ocp, 
                                                   const double kkt_error, 
                                                   KKTMatrix& kkt_matrix) const {
  const double reg = getRegularization(kkt_error);
  for (int i=0; i<ocp.discrete().N_impulse(); ++i) {
    kkt_matrix.aux[i].Qtt += reg;
  }
  for (int i=0; i<ocp.discrete().N_lift(); ++i) {
    kkt_matrix.lift[i].Qtt += reg;
  }
}


inline bool STORegularization::isRegularizationValid() const {
  return is_reg_valid_;
}


inline double STORegularization::getRegularization(const double kkt_error) const {
  switch (reg_type_) {
    case STORegularizationType::Const:
      return w_;
    case STORegularizationType::Abs:
      return (w_ * std::abs(kkt_error));
    case STORegularizationType::Quad:
      return (w_ * kkt_error*kkt_error);
    case STORegularizationType::Exp:
      return (w_ * std::exp(std::abs(kkt_error)));
    case STORegularizationType::Exp2:
      return (w_ * std::exp(kkt_error*kkt_error));
    default:
      return 0.0;
  }
}

} // namespace robotoc 

#endif // ROBOTOC_STO_REGULARIZATION_HXX_ 
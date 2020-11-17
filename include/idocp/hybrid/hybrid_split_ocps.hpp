#ifndef IDOCP_HYBRID_SPLIT_OCPS_HPP_
#define IDOCP_HYBRID_SPLIT_OCPS_HPP_

#include <vector>

namespace idocp {

///
/// @class HybridSplitOCPs
/// @brief A container that is useful to formulate the hybrid optimal control 
/// problem. This container has the standard data (with Type), data for lift 
/// stages (with Type), data for aux stages (with Type), and data for impulse 
/// stages (with ImpulseType).
/// @tparam Type The type name of the standard data type.
/// @tparam ImpulseType The type name of the impulse data type.
/// 
///
struct HybridSplitOCPs {
  ///
  /// @brief Construct the standard data, impulse data, and lift data. 
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  ///
  HybridSplitOCPs(const int N, const int N_impulse) 
    : data(N, Type()), 
      aux(N_impulse, Type()), 
      lift(N_impulse, Type()),
      impulse(N_impulse, ImpulseType()) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  /// @param[in] obj An object of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  /// @param[in] impulse_obj An object of the impulse data.
  ///
  HybridSplitOCPs(const int N, const Type& obj, const int N_impulse, 
                   const ImpulseType& impulse_obj) 
    : data(N, obj), 
      aux(N_impulse, obj), 
      lift(N_impulse, obj),
      impulse(N_impulse, impulse_obj) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  ///
  HybridSplitOCPs(const int N) 
    : data(N, Type()), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  /// @param[in] obj An object of the standard data.
  ///
  HybridSplitOCPs(const int N, const Type& obj) 
    : data(N, obj), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Default Constructor.
  ///
  HybridSplitOCPs() 
    : data(), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  HybridSplitOCPs(const HybridSplitOCPs&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  HybridSplitOCPs& operator=(const HybridSplitOCPs&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  HybridSplitOCPs(HybridSplitOCPs&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  HybridSplitOCPs& operator=(HybridSplitOCPs&&) noexcept = default;

  ///
  /// @brief Overload operator[] to access the standard data, i.e., 
  /// HybridSplitOCPs::data as std::vector. 
  ///
  Type& operator[] (const int i) {
    return data[i];
  }

  ///
  /// @brief const version of HybridSplitOCPs::operator[]. 
  ///
  const Type& operator[] (const int i) const {
    return data[i];
  }

  std::vector<SplitOCP> ocps, aux_ocps, lift_ocps;
  std::vector<SplitImpulseOCP> impulse_ocps;
};
  
} // namespace idocp

#endif // IDOCP_HYBRID_SPLIT_OCP_HPP_ 
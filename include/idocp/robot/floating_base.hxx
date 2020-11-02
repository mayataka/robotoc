#ifndef IDOCP_FLOATING_BASE_HXX_
#define IDOCP_FLOATING_BASE_HXX_

#include "idocp/robot/floating_base.hpp"

namespace idocp {

inline int FloatingBase::dim_passive() const {
  return dim_passive_;
}


inline std::vector<int> FloatingBase::passive_joint_indices() const {
  return passive_joint_indices_;
}


inline bool FloatingBase::has_floating_base() const {
  return has_floating_base_;
}

} // namespace idocp

#endif // IDOCP_FLOATING_BASE_HXX_
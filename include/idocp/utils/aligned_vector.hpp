#ifndef IDOCP_ALIGNED_VECTOR_HPP_
#define IDOCP_ALIGNED_VECTOR_HPP_

#include <vector>
#include "Eigen/Core"

namespace idocp {

template <typename T>
using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

} // namespace idocp

#endif // IDOCP_ALIGNED_VECTOR_HPP_
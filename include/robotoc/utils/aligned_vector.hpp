#ifndef ROBOTOC_ALIGNED_VECTOR_HPP_
#define ROBOTOC_ALIGNED_VECTOR_HPP_

#include <vector>
#include "Eigen/StdVector"

namespace robotoc {

///
/// @typedef aligned_vector
/// @brief std vector with Eigen::aligned_allocator. 
///
template <typename T>
using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

} // namespace robotoc

#endif // ROBOTOC_ALIGNED_VECTOR_HPP_
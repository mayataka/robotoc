#ifndef ROBOTOC_ALIGNED_DEQUE_HPP_
#define ROBOTOC_ALIGNED_DEQUE_HPP_

#include <deque>
#include "Eigen/StdDeque"

namespace robotoc {

///
/// @typedef aligned_deque
/// @brief std deque with Eigen::aligned_allocator. 
///
template <typename T>
using aligned_deque = std::deque<T, Eigen::aligned_allocator<T>>;

} // namespace robotoc

#endif // ROBOTOC_ALIGNED_DEQUE_HPP_
#ifndef ROBOTOC_ALIGNED_UNORDERED_MAP_HPP_
#define ROBOTOC_ALIGNED_UNORDERED_MAP_HPP_

#include <unordered_map>
#include "Eigen/StdVector"

namespace robotoc {

///
/// @typedef aligned_unordered_map
/// @brief std unordered_map with Eigen::aligned_allocator. 
///
template <typename Key, typename T>
using aligned_unordered_map = std::unordered_map<Key, T, std::hash<Key>, std::equal_to<Key>,
                                                 Eigen::aligned_allocator<std::pair<Key, T>>>;

} // namespace robotoc

#endif // ROBOTOC_ALIGNED_UNORDERED_MAP_HPP_
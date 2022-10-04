#ifndef ROBOTOC_KKT_MATRIX_HPP_
#define ROBOTOC_KKT_MATRIX_HPP_

#include <iostream>

#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @typedef KKTMatrix 
/// @brief The KKT matrix of the optimal control problem. 
///
using KKTMatrix = aligned_vector<SplitKKTMatrix>;

std::ostream& operator<<(std::ostream& os, const KKTMatrix& kkt_matrix);

} // namespace robotoc

#endif // ROBOTOC_KKT_MATRIX_HPP_ 
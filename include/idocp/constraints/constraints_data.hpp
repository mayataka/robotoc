#ifndef IDOCP_CONSTRAINTS_DATA_HPP_
#define IDOCP_CONSTRAINTS_DATA_HPP_

#include <vector>

#include "idocp/constraints/constraint_component_data.hpp"


namespace idocp {

class ConstraintsData {
public:
  ConstraintsData() 
    : data() {
  }

  ~ConstraintsData() {
  }

  // Use default copy constructor.
  ConstraintsData(const ConstraintsData&) = default;

  // Use default copy coperator.
  ConstraintsData& operator=(const ConstraintsData&) = default;

  // Use default move constructor.
  ConstraintsData(ConstraintsData&&) noexcept = default;

  // Use default move assign coperator.
  ConstraintsData& operator=(ConstraintsData&&) noexcept = default;

  std::vector<ConstraintComponentData> data;

};

} // namespace idocp


#endif // IDOCP_CONSTRAINTS_DATA_HPP_
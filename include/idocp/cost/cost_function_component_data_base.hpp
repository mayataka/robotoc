#ifndef IDOCP_COST_FUNCTION_COMPONENT_DATA_BASE_HPP_
#define IDOCP_COST_FUNCTION_COMPONENT_DATA_BASE_HPP_


namespace idocp {

class CostFunctionComponentDataBase {
public:
  CostFunctionComponentDataBase() {}

  virtual ~CostFunctionComponentDataBase() {}

  // Use default copy constructor.
  CostFunctionComponentDataBase(const CostFunctionComponentDataBase&) = default;

  // Use default copy coperator.
  CostFunctionComponentDataBase& operator=(const CostFunctionComponentDataBase&) 
      = default;

  // Use default move constructor.
  CostFunctionComponentDataBase(CostFunctionComponentDataBase&&) noexcept 
      = default;

  // Use default move assign coperator.
  CostFunctionComponentDataBase& operator=(CostFunctionComponentDataBase&&) noexcept 
      = default;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_COMPONENT_DATA_BASE_HPP_ 
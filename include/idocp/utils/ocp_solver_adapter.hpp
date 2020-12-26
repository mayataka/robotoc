#ifndef IDOCP_OCP_SOLVER_ADAPTER_HPP_
#define IDOCP_OCP_SOLVER_ADAPTER_HPP_

#include <string>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class OCPSolverAdapter {
public:
  OCPSolverAdapter();
  virtual ~OCPSolverAdapter() = 0;
  
private:

};

} // namespace idocp 

#endif // IDOCP_OCP_SOLVER_ADAPTER_HPP_ 
#include "robotoc/solver/solver_statistics.hpp"


namespace robotoc {

SolverStatistics::SolverStatistics()
  : convergence(false),
    iter(0),
    kkt_error(),
    primal_step_size(),
    dual_step_size(),
    ts(),
    mesh_refinement_iter() {
}


SolverStatistics::~SolverStatistics() {
}


void SolverStatistics::clear() {
  convergence = false;
  iter = 0;
  kkt_error.clear();
  primal_step_size.clear();
  dual_step_size.clear();
  ts.clear();
  mesh_refinement_iter.clear();
}


void SolverStatistics::disp(std::ostream& os) const {
  os << "Solver statistics:" << std::endl;
  os << "  convergence: " << std::boolalpha << convergence << std::endl;
  os << "  total no. of iteration: " << iter << std::endl;
  for (int i=0; i<iter; ++i) {
    os << "  iter: " << i+1 << std::endl;
    os << "    kkt_error:        " << kkt_error[i] << std::endl;
    os << "    primal_step_size: " << primal_step_size[i] << std::endl;
    os << "    dual_step_size:   " << dual_step_size[i] << std::endl;
    if (ts.size() > 0) {
    os << "    ts:               [";
      for (int j=0; j<ts[i].size()-1; ++j) {
        os << ts[i][j] << ", ";
      }
      os << ts[i][ts[i].size()-1] << "]" << std::endl;
    }
  }
  if (mesh_refinement_iter.size() > 0) {
    os << "  mesh_refinement_iter: [";
    for (int i=0; i<mesh_refinement_iter.size()-1; ++i) {
      os << mesh_refinement_iter[i] << ", ";
    }
    os << mesh_refinement_iter[mesh_refinement_iter.size()-1] << "]" << std::flush;
  }
}


std::ostream& operator<<(std::ostream& os, 
                         const SolverStatistics& solver_statistics) {
  solver_statistics.disp(os);
  return os;
}

} // namespace robotoc
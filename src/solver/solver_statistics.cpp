#include "robotoc/solver/solver_statistics.hpp"

#include <iomanip>
#include <algorithm>


namespace robotoc {

void SolverStatistics::clear() {
  convergence = false;
  iter = 0;
  kkt_error.clear();
  primal_step_size.clear();
  dual_step_size.clear();
  ts.clear();
  mesh_refinement_iter.clear();
  cpu_time = 0.0;
}


void SolverStatistics::disp(std::ostream& os) const {
  os << "Solver statistics:" << std::endl;
  os << "  convergence: " << std::boolalpha << convergence << std::endl;
  os << "  total no. of iteration: " << iter << std::endl;
  os << "  CPU time (positive if benchmark is enabled): " << cpu_time << std::endl;
  os << "  ------------------------------------------------------------------------------------ " << std::endl;
  os << "   iteration |      KKT error | primal step size |   dual step size |        ts        " << std::endl;
  os << "  ------------------------------------------------------------------------------------ " << std::endl;
  for (int i=0; i<iter; ++i) {
    if (std::find(mesh_refinement_iter.begin(), mesh_refinement_iter.end(), i) != mesh_refinement_iter.end()) {
      os << "        ======================== Mesh-refinement is carried out! ========================" << std::endl;
    }
    os << "         " << std::setw(3) << i+1 << " | ";
    os << std::scientific << std::setprecision(6);
    os << "  " << kkt_error[i];
    os << " |     " << primal_step_size[i];
    os << " |     " << dual_step_size[i] << " |";
    if (ts.size() > 0) {
      os << "  [";
      os << std::fixed << std::setprecision(4);
      for (int j=0; j<ts[i].size()-1; ++j) {
        os << std::setw(6) << ts[i][j] << ", ";
      }
      os << std::setw(6) << ts[i][ts[i].size()-1] << "]";
    }
    os << std::endl;
  }
  os << std::defaultfloat << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SolverStatistics& solver_statistics) {
  solver_statistics.disp(os);
  return os;
}

} // namespace robotoc
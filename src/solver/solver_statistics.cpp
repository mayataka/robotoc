#include "robotoc/solver/solver_statistics.hpp"

#include <cassert>
#include <iomanip>
#include <algorithm>
#include <cmath>


namespace robotoc {

void SolverStatistics::reserve(const int size) {
  assert(size >= 0);
  performance_index.reserve(size);
  primal_step_size.reserve(size);
  dual_step_size.reserve(size);
  ts.reserve(size);
  mesh_refinement_iter.reserve(size);
}


void SolverStatistics::clear() {
  convergence = false;
  iter = 0;
  performance_index.clear();
  primal_step_size.clear();
  dual_step_size.clear();
  ts.clear();
  mesh_refinement_iter.clear();
  cpu_time = 0.0;
}


void SolverStatistics::disp(std::ostream& os) const {
  os << "Solver statistics:" << "\n";
  os << "  convergence: " << std::boolalpha << convergence << "\n";
  os << "  total No. of iterations: " << iter << "\n";
  os << "  CPU time: " << std::setprecision(3) << cpu_time << " ms (non-zero if benchmark is enabled) \n";
  os << "  ------------------------------------------------------------------------------------------------------------------ " << "\n";
  os << "   iter |   KKT error  |      cost    |  primal feas |   dual feas  | primal alpha |   dual alpha |        ts        " << "\n";
  os << "  ------------------------------------------------------------------------------------------------------------------ " << "\n";
  for (int i=0; i<iter; ++i) {
    if (std::find(mesh_refinement_iter.begin(), mesh_refinement_iter.end(), i) != mesh_refinement_iter.end()) {
      os << "  ========================================= Mesh-refinement is carried out! ========================================= " << "\n";
    }
    os << "    " << std::setw(3) << i+1;
    os << std::scientific << std::setprecision(3);
    os << " |    " << std::sqrt(performance_index[i].kkt_error);
    const double cost = performance_index[i].cost + performance_index[i].cost_barrier;
    if (cost >= 0) os << " |    " << cost;
    else os << " |   " << cost;
    os << " |    " << performance_index[i].primal_feasibility;
    os << " |    " << performance_index[i].dual_feasibility;
    os << " |    " << primal_step_size[i];
    os << " |    " << dual_step_size[i] << " |";
    if (ts.size() > 0) {
      os << "  [";
      os << std::fixed << std::setprecision(4);
      for (int j=0; j<ts[i].size()-1; ++j) {
        os << std::setw(6) << ts[i][j] << ", ";
      }
      os << std::setw(6) << ts[i][ts[i].size()-1] << "]";
    }
    os << "\n";
  }
  os << std::defaultfloat << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SolverStatistics& solver_statistics) {
  solver_statistics.disp(os);
  return os;
}

} // namespace robotoc
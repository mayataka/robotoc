#ifndef ROBOTOC_LOGGER_HPP_
#define ROBOTOC_LOGGER_HPP_

#include <string>
#include <vector>

#include "Eigen/Core"

#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/solver/unconstr_ocp_solver.hpp"
#include "robotoc/solver/unconstr_parnmpc_solver.hpp"


namespace robotoc {

///
/// @class Logger
/// @brief Logger for optimal control solvers.
///
class Logger {
public:
  ///
  /// @brief Constructs the logger.
  /// @param[in] vars Names of the variables to be logged. Choose from "q", "v",
  /// "a", "f", "u", "ts", and "KKT".
  /// @param[in] save_dir Path to the directory where the log files are 
  /// saved. Default is "./".
  ///
  Logger(const std::vector<std::string>& vars, 
         const std::string& save_dir="./");

  ///
  /// @brief Destructor.
  ///
  ~Logger();

  ///
  /// @brief Takes the log.
  /// @param[in] solver The OCP solver.
  /// @param[in] precision Precision of float values. Default is 
  /// Eigen::FullPrecision.
  /// @param[in] delimiter Delimiter of values. Default is ", ".
  ///
  void takeLog(OCPSolver& solver, const int precision=Eigen::FullPrecision,
               const std::string& delimiter=", ");

  ///
  /// @brief Takes the log.
  /// @param[in] solver The constrained OCP solver.
  /// @param[in] precision Precision of float values. Default is 
  /// Eigen::FullPrecision.
  /// @param[in] delimiter Delimiter of values. Default is ", ".
  ///
  void takeLog(UnconstrOCPSolver& solver, 
               const int precision=Eigen::FullPrecision,
               const std::string& delimiter=", ");

  ///
  /// @brief Takes the log.
  /// @param[in] solver The constrained ParNMPC solver.
  /// @param[in] precision Precision of float values. Default is 
  /// Eigen::FullPrecision.
  /// @param[in] delimiter Delimiter of values. Default is ", ".
  ///
  void takeLog(UnconstrParNMPCSolver& solver, 
               const int precision=Eigen::FullPrecision,
               const std::string& delimiter=", ");

private:
  std::vector<std::string> vars_;
  std::string save_dir_;
  std::vector<std::ofstream> logs_;

  template <typename UnconstrOCPSolverType>
  void takeLog_unconstr_impl(UnconstrOCPSolverType& solver, 
                             const int precision,
                             const std::string& delimiter) {
    Eigen::IOFormat format(precision, 0, delimiter);
    for (int i=0; i<vars_.size(); ++i) {
      if (vars_[i] == "KKT") {
        const double kkt = solver.KKTError();
        Eigen::VectorXd vec(1); vec << kkt;
        logs_[i] << vec.format(format) << "\n";
      }
      else {
        const auto sols = solver.getSolution(vars_[i]);
        for (const auto& e : sols) {
          logs_[i] << e.transpose().format(format) << "\n";
        }
        logs_[i] << "\n";
      }
    }
  }

};

} // namespace robotoc

#endif // ROBOTOC_LOGGER_HPP_
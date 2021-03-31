#ifndef IDOCP_LOGGER_HPP_
#define IDOCP_LOGGER_HPP_

#include <string>
#include <vector>

#include "Eigen/Core"


namespace idocp {

///
/// @class Logger
/// @brief Logger of the solution of the optimal control problems.
///
class Logger {
public:
  ///
  /// @brief Constructs the logger.
  /// @param[in] save_dir Path to the directory where the log files are 
  /// saved. Default is "./".
  ///
  Logger(const std::string& save_dir="./");

  ///
  /// @brief Saves the log.
  /// @param[in] file_name File name.
  /// @param[in] var Std vector of variables.
  /// @param[in] precision Precision of float values. Default is 
  /// Eigen::FullPrecision.
  /// @param[in] delimiter Delimiter of values. Default is ", ".
  ///
  void save(const std::string& file_name, 
            const std::vector<Eigen::VectorXd>& var,
            const int precision=Eigen::FullPrecision,
            const std::string& delimiter=", ") const;

private:
  std::string save_dir_;

};

} // namespace idocp


#endif // IDOCP_LOGGER_HPP_
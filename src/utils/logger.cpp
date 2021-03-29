#include "idocp/utils/logger.hpp"

#include <fstream>


namespace idocp {

Logger::Logger(const std::string& save_dir) 
  : save_dir_(save_dir) {
}

void Logger::save(const std::string& file_name, 
                  const std::vector<Eigen::VectorXd>& var, const int precision, 
                  const std::string& delimiter) const {
 	Eigen::IOFormat format(precision, 0, delimiter);
  std::ofstream log(save_dir_+file_name);
  for (const auto& v : var) {
    log << v.transpose().format(format) << "\n";
  }
  log.close();
}

} // namespace idocp
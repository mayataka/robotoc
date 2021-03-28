#include "idocp/utils/logger.hpp"

#include <fstream>


namespace idocp {

Logger::Logger(const std::string& path_to_file) 
  : path_to_file_(path_to_file) {
}

void Logger::save(const std::string& file_name, 
                  const std::vector<Eigen::VectorXd>& var) const {
 	Eigen::IOFormat format(Eigen::FullPrecision, 0, ", ");
  std::ofstream log(path_to_file_+file_name);
  for (const auto& v : var) {
    log << v.transpose().format(format) << "\n";
  }
  log.close();
}

} // namespace idocp
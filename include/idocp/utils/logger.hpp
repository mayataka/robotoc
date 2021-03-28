#ifndef IDOCP_LOGGER_HPP_
#define IDOCP_LOGGER_HPP_

#include <string>
#include <vector>

#include "Eigen/Core"


namespace idocp {

class Logger {
public:
  Logger(const std::string& path_to_file="./");

  void save(const std::string& file_name, 
            const std::vector<Eigen::VectorXd>& var) const;

private:
  std::string path_to_file_;

};

} // namespace idocp


#endif // IDOCP_LOGGER_HPP_
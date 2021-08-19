#include "idocp/utils/logger.hpp"

#include <fstream>
#include <boost/filesystem.hpp>


namespace idocp {

Logger::Logger(const std::vector<std::string>& vars,
               const std::string& save_dir) 
  : vars_(vars),
    save_dir_(save_dir+"_log"),
    logs_() {
  namespace fs = boost::filesystem;
  const fs::path path = save_dir + "_log";
  if (!fs::exists(path)) {
    fs::create_directory(path);
  }
  for (const auto& var : vars) {
    const std::string log_file_path = save_dir + "_log" + "/" + var + ".log";
    if (fs::exists(fs::path(log_file_path))) {
      fs::remove(fs::path(log_file_path));
    }
    logs_.push_back(std::ofstream(log_file_path));
  }
}


Logger::~Logger() {
  for (auto& log : logs_) {
    log.close();
  }
}


void Logger::takeLog(OCPSolver& solver, const int precision,
                     const std::string& delimiter) {
  Eigen::IOFormat format(precision, 0, delimiter);
  for (int i=0; i<vars_.size(); ++i) {
    if (vars_[i] == "KKT") {
      const double kkt = solver.KKTError();
      Eigen::VectorXd vec(1); vec << kkt;
      logs_[i] << vec.format(format) << "\n";
    }
    else if (vars_[i] == "f") {
      const auto fs = solver.getSolution(vars_[i], "WORLD");
      for (const auto& e : fs) {
        logs_[i] << e.transpose().format(format) << "\n";
      }
      logs_[i] << "\n";
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


void Logger::takeLog(UnconstrOCPSolver& solver, const int precision,
                     const std::string& delimiter) {
  takeLog_unconstr_impl(solver, precision, delimiter);
}


void Logger::takeLog(UnconstrParNMPCSolver& solver, const int precision,
                     const std::string& delimiter) {
  takeLog_unconstr_impl(solver, precision, delimiter);
}

} // namespace idocp
#include "idocp/utils/csv_to_eigen.hpp"


namespace idocp {

CSVToEigen::CSVToEigen(const std::string& file_name, const int dim) 
  : ifs_(file_name),
    row_(),
    dim_(dim) {
}


CSVToEigen::~CSVToEigen() {
  ifs_.close();
}


bool CSVToEigen::readCSVLine() {
  std::string line;  
  if (!std::getline(ifs_, line)) {
    return false;
  }
  std::string value;  
  vec_.clear();
  for (std::stringstream ss(line); std::getline(ss, value, ' '); ) {
    vec_.push_back(std::stod(value));
  }
  if (vec_.size() != dim_) {
    return false;
  }
  return true;
}


Eigen::VectorXd CSVToEigen::get() const {
  return Eigen::VectorXd::Map(vec_.data(), vec_.size());
}

} // namespace idocp
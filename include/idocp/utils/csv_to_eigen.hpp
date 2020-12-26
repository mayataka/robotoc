#ifndef IDOCP_CSV_TO_EIGEN_HPP_
#define IDOCP_CSV_TO_EIGEN_HPP_

#include <string>
#include <fstream>
#include <vector>

#include "Eigen/Core"


namespace idocp {

class CSVToEigen {
public:
  explicit CSVToEigen(const std::string& file_name, const int dim);
  ~CSVToEigen();

  bool readCSVLine();

  Eigen::VectorXd get() const;

private:
  std::ifstream ifs_;
  std::string row_;
  int dim_;
  std::vector<double> vec_;

}; 

} // namespace idocp 

#endif // IDOCP_CSV_TO_EIGEN_HPP_ 
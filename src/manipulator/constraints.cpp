#include "idocp/manipulator/constraints.hpp"


namespace idocp {
namespace manipulator {

Constraints::Constraints(const Robot& robot)
  : ConstraintsInterface() {
}


Constraints::Constraints()
  : ConstraintsInterface() {
}


Constraints::~Constraints() {
}


bool Constraints::isFeasible(const Robot& robot, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                             const Eigen::VectorXd& u) {
}


void Constraints::setSlackAndDual(const Robot& robot, const double dtau, 
                                  const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& v, 
                                  const Eigen::VectorXd& a, 
                                  const Eigen::VectorXd& u) {
}


void Constraints::augmentDualResidual(const Robot& robot, const double dtau, 
                                      Eigen::VectorXd& Cu) {
}


void Constraints::augmentDualResidual(const Robot& robot, const double dtau, 
                                      Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                                      Eigen::VectorXd& Ca) {
}


void Constraints::condenseSlackAndDual(const Robot& robot, const double dtau, 
                                       const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const Eigen::VectorXd& a, 
                                       Eigen::MatrixXd& Cqq, 
                                       Eigen::MatrixXd& Cqv, 
                                       Eigen::MatrixXd& Cqa,  
                                       Eigen::MatrixXd& Cvv, 
                                       Eigen::MatrixXd& Cva,  
                                       Eigen::MatrixXd& Caa,  
                                       Eigen::VectorXd& Cq, 
                                       Eigen::VectorXd& Cv, 
                                       Eigen::VectorXd& Ca) {
}


void Constraints::condenseSlackAndDual(const Robot& robot, const double dtau, 
                                       const Eigen::VectorXd& u, 
                                       Eigen::MatrixXd& Cuu, 
                                       Eigen::VectorXd& Cu) {
}


void Constraints::computeSlackAndDualDirection(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& dq,
                                               const Eigen::VectorXd& dv, 
                                               const Eigen::VectorXd& da, 
                                               const Eigen::VectorXd& du) {
}


double Constraints::maxSlackStepSize() {
}


double Constraints::maxDualStepSize() {
}


void Constraints::updateSlack(const double step_size) {
}


void Constraints::updateDual(const double step_size) {
}


double Constraints::costSlackBarrier() {
}


double Constraints::costSlackBarrier(const double step_size) {
}


double Constraints::residualL1Nrom(const Robot& robot, const double dtau,
                                   const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v, 
                                   const Eigen::VectorXd& a, 
                                   const Eigen::VectorXd& u) {
}


double Constraints::residualSquaredNrom(const Robot& robot, const double dtau,
                                        const Eigen::VectorXd& q, 
                                        const Eigen::VectorXd& v, 
                                        const Eigen::VectorXd& a, 
                                        const Eigen::VectorXd& u) {
}

} // namespace manipulator
} // namespace idocp
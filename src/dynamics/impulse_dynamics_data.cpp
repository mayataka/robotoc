#include "robotoc/dynamics/impulse_dynamics_data.hpp"

namespace robotoc {

ImpulseDynamicsData::ImpulseDynamicsData(const Robot& robot) 
  : dImDddv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                         2*robot.dimv())), 
    MJtJinv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJtJinv_dImDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                                 2*robot.dimv())), 
    Qdvfqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                       2*robot.dimv())), 
    ImDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    MJtJinv_ImDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    ldvf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()),
    dimf_(0),
    dimvf_(robot.dimv()) {
}


ImpulseDynamicsData::ImpulseDynamicsData() 
  : dImDddv(),
    dImDCdqv_full_(), 
    MJtJinv_full_(), 
    MJtJinv_dImDCdqv_full_(), 
    Qdvfqv_full_(), 
    ImDC_full_(),
    MJtJinv_ImDC_full_(),
    ldvf_full_(),
    dimv_(0),
    dimf_(0),
    dimvf_(0) {
}

} // namespace robotoc 
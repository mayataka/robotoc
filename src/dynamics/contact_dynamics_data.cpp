#include "robotoc/dynamics/contact_dynamics_data.hpp"

namespace robotoc {

ContactDynamicsData::ContactDynamicsData(const Robot& robot) 
  : Qxu_passive(Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dim_passive())),
    Quu_passive_topRight(Eigen::MatrixXd::Zero(robot.dim_passive(), robot.dimu())),
    lu_passive(Eigen::VectorXd::Zero(robot.dim_passive())),
    dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dIDddv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dCda_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    dIDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        2*robot.dimv())),
    MJtJinv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJtJinv_dIDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                                2*robot.dimv())), 
    Qafqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                      2*robot.dimv())), 
    Qafu_full_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                          robot.dimv())), 
    IDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    MJtJinv_IDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    Phia_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    laf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    haf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    dimf_(0),
    dimvf_(robot.dimv()),
    dims_(0),
    dim_passive_(robot.dim_passive()),
    has_floating_base_(robot.hasFloatingBase()) {
}


ContactDynamicsData::ContactDynamicsData() 
  : Qxu_passive(),
    Quu_passive_topRight(),
    lu_passive(),
    dIDda(),
    dIDddv(),
    dIDCdqv_full_(),
    MJtJinv_full_(), 
    MJtJinv_dIDCdqv_full_(), 
    Qafqv_full_(), 
    Qafu_full_full_(), 
    IDC_full_(),
    MJtJinv_IDC_full_(),
    Phia_full_(),
    laf_full_(),
    haf_full_(),
    dimv_(0),
    dimu_(0),
    dimf_(0),
    dimvf_(0),
    dims_(0),
    dim_passive_(0),
    has_floating_base_(false) {
}

} // namespace robotoc 
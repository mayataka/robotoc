#ifndef IDOCP_KKT_COMPOSITION_HPP_
#define IDOCP_KKT_COMPOSITION_HPP_

#include <vector>
#include <memory>

#include "idocp/robot/robot.hpp"


namespace idocp {

class KKTComposition {
public:
  KKTComposition(const Robot& robot)
    : dimv_(robot.dimv()), 
      max_dimf_(robot.max_dimf()), 
      dimf_(robot.dimf()), 
      dim_passive_(robot.dim_passive()), 
      max_dimc_(robot.max_dimf()+robot.dim_passive()), 
      dimc_(robot.dimf()+robot.dim_passive()),
      max_dimKKT_(5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive()),
      dimKKT_(5*robot.dimv()+2*robot.dimf()+robot.dim_passive()),
      Fq_begin_(0),
      Fv_begin_(robot.dimv()),
      C_begin_(2*robot.dimv()),
      Qa_begin_(2*robot.dimv()+robot.dimf()+robot.dim_passive()),
      Qf_begin_(3*robot.dimv()+robot.dimf()+robot.dim_passive()),
      Qq_begin_(3*robot.dimv()+2*robot.dimf()+robot.dim_passive()),
      Qv_begin_(4*robot.dimv()+2*robot.dimf()+robot.dim_passive()) {
  }

  KKTComposition() 
    : dimv_(0), 
      max_dimf_(0), 
      dimf_(0), 
      dim_passive_(0), 
      max_dimc_(0), 
      dimc_(0),
      max_dimKKT_(0),
      dimKKT_(0),
      Fq_begin_(0),
      Fv_begin_(0),
      C_begin_(0),
      Qa_begin_(0),
      Qf_begin_(0),
      Qq_begin_(0),
      Qv_begin_(0) {
  }

  ~KKTComposition() {
  }

  KKTComposition(const KKTComposition&) = default;

  KKTComposition& operator=(const KKTComposition&) = default;
 
  KKTComposition(KKTComposition&&) noexcept = default;

  KKTComposition& operator=(KKTComposition&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    dimf_ = robot.dimf();
    dimc_ = robot.dimf() + robot.dim_passive();
    dimKKT_ = 5*robot.dimv() + 2*robot.dimf() + robot.dim_passive();
    Qa_begin_ = 2*robot.dimv() + robot.dimf() + robot.dim_passive();
    Qf_begin_ = 3*robot.dimv() + robot.dimf() + robot.dim_passive();
    Qq_begin_ = 3*robot.dimv() + 2*robot.dimf() + robot.dim_passive();
    Qv_begin_ = 4*robot.dimv() + 2*robot.dimf() + robot.dim_passive();
  }

  inline int max_dimKKT() const {
    return max_dimKKT_;
  }

  inline int dimKKT() const {
    return dimKKT_;
  }

  inline int Fq_begin() const {
    return Fq_begin_;
  }

  inline int Fq_size() const {
    return dimv_;
  }

  inline int Fv_begin() const {
    return Fv_begin_;
  }

  inline int Fv_size() const {
    return dimv_;
  }

  inline int C_begin() const {
    return C_begin_;
  }

  inline int C_size() const {
    return dimc_;
  }

  inline int Qa_begin() const {
    return Qa_begin_;
  }

  inline int Qa_size() const {
    return dimv_;
  }

  inline int Qf_begin() const {
    return Qf_begin_;
  }

  inline int Qf_size() const {
    return dimf_;
  }

  inline int Qq_begin() const {
    return Qq_begin_;
  }

  inline int Qq_size() const {
    return dimv_;
  }

  inline int Qv_begin() const {
    return Qv_begin_;
  }

  inline int Qv_size() const {
    return dimv_;
  }

  inline int Fx_begin() const {
    return Fq_begin_;
  }

  inline int Fx_size() const {
    return 2*dimv_;
  }

  inline int Qx_begin() const {
    return Qq_begin_;
  }

  inline int Qx_size() const {
    return 2*dimv_;
  }

  inline int Q_begin() const {
    return Qq_begin_;
  }

  inline int Q_size() const {
    return 3*dimv_+dimf_;
  }

private:
  int dimv_, max_dimf_, dimf_, dim_passive_, max_dimc_, dimc_,  
      max_dimKKT_, dimKKT_, Fq_begin_, Fv_begin_, C_begin_, 
      Qa_begin_, Qf_begin_, Qq_begin_, Qv_begin_, Qx_begin_;

};

} // namespace idocp 


#endif // IDOCP_KKT_COMPOSITION_HPP_
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/cost/time_varying_task_space_3d_cost.hpp"


namespace idocp {
namespace python {

class PyTimeVaryingTaskSpace3DRefBase : public TimeVaryingTaskSpace3DRefBase {
public:
  // Inherit the constructors
  using TimeVaryingTaskSpace3DRefBase::TimeVaryingTaskSpace3DRefBase;

  void update_q_3d_ref(const double t, Eigen::VectorXd& q_3d_ref) const override {
    PYBIND11_OVERRIDE_PURE(void, TimeVaryingTaskSpace3DRefBase, 
                           update_q_3d_ref, 
                           t, q_3d_ref);
  }

  bool isActive(const double t) const override {
    PYBIND11_OVERRIDE_PURE(bool, TimeVaryingTaskSpace3DRefBase, 
                           isActive, t);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(time_varying_task_space_3d_ref_base, m) {
  py::class_<TimeVaryingTaskSpace3DRefBase, 
             PyTimeVaryingTaskSpace3DRefBase,
             std::shared_ptr<TimeVaryingTaskSpace3DRefBase>>(m, "TimeVaryingTaskSpace3DRefBase")
    .def(py::init<>())
    .def("update_q_3d_ref", &TimeVaryingTaskSpace3DRefBase::update_q_3d_ref)
    .def("isActive", &TimeVaryingTaskSpace3DRefBase::isActive);
}

} // namespace python
} // namespace idocp
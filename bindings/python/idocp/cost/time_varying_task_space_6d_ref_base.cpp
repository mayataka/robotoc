#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/cost/time_varying_task_space_6d_cost.hpp"
#include "idocp/robot/se3.hpp"


namespace idocp {
namespace python {

class PyTimeVaryingTaskSpace6DRefBase : public TimeVaryingTaskSpace6DRefBase {
public:
  // Inherit the constructors
  using TimeVaryingTaskSpace6DRefBase::TimeVaryingTaskSpace6DRefBase;

  void update_SE3_ref(const double t, SE3& SE3_ref) const override {
    PYBIND11_OVERRIDE_PURE(void, TimeVaryingTaskSpace6DRefBase, 
                           update_SE3_ref, 
                           t, SE3_ref);
  }

  bool isActive(const double t) const override {
    PYBIND11_OVERRIDE_PURE(bool, TimeVaryingTaskSpace6DRefBase, 
                           isActive, t);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(time_varying_task_space_6d_ref_base, m) {
  py::class_<TimeVaryingTaskSpace6DRefBase, 
             PyTimeVaryingTaskSpace6DRefBase,
             std::shared_ptr<TimeVaryingTaskSpace6DRefBase>>(m, "TimeVaryingTaskSpace6DRefBase")
    .def(py::init<>())
    .def("update_SE3_ref", &TimeVaryingTaskSpace6DRefBase::update_SE3_ref)
    .def("isActive", &TimeVaryingTaskSpace6DRefBase::isActive);
}

} // namespace python
} // namespace idocp
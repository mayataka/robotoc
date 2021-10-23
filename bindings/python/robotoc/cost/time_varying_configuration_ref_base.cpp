#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/time_varying_configuration_space_cost.hpp"


namespace robotoc {
namespace python {

class PyTimeVaryingConfigurationRefBase : public TimeVaryingConfigurationRefBase {
public:
  // Inherit the constructors
  using TimeVaryingConfigurationRefBase::TimeVaryingConfigurationRefBase;

  void update_q_ref(const Robot& robot, const double t, 
                    Eigen::VectorXd& q_ref) const override {
    PYBIND11_OVERRIDE_PURE(void, TimeVaryingConfigurationRefBase, 
                           update_q_ref, 
                           robot, t, q_ref);
  }

  bool isActive(const double t) const override {
    PYBIND11_OVERRIDE_PURE(bool, TimeVaryingConfigurationRefBase, 
                           isActive, t);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(time_varying_configuration_ref_base, m) {
  py::class_<TimeVaryingConfigurationRefBase, 
             PyTimeVaryingConfigurationRefBase,
             std::shared_ptr<TimeVaryingConfigurationRefBase>>(m, "TimeVaryingConfigurationRefBase")
    .def(py::init<>())
    .def("update_q_ref", &TimeVaryingConfigurationRefBase::update_q_ref)
    .def("isActive", &TimeVaryingConfigurationRefBase::isActive);
}

} // namespace python
} // namespace robotoc
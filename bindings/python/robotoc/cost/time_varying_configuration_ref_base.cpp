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

  void update_q_ref(const Robot& robot, const GridInfo& grid_info,
                    Eigen::VectorXd& q_ref) const override {
    PYBIND11_OVERRIDE_PURE(void, TimeVaryingConfigurationRefBase, 
                           update_q_ref, 
                           robot, grid_info, q_ref);
  }

  bool isActive(const GridInfo& grid_info) const override {
    PYBIND11_OVERRIDE_PURE(bool, TimeVaryingConfigurationRefBase, 
                           isActive, grid_info);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(time_varying_configuration_ref_base, m) {
  py::class_<TimeVaryingConfigurationRefBase, 
             PyTimeVaryingConfigurationRefBase,
             std::shared_ptr<TimeVaryingConfigurationRefBase>>(m, "TimeVaryingConfigurationRefBase")
    .def(py::init<>())
    .def("update_q_ref", &TimeVaryingConfigurationRefBase::update_q_ref,
          py::arg("robot"), py::arg("grid_info"), py::arg("q_ref"))
    .def("isActive", &TimeVaryingConfigurationRefBase::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
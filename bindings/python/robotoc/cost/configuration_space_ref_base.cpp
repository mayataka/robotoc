#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/configuration_space_ref_base.hpp"


namespace robotoc {
namespace python {

class PyConfigurationSpaceRefBase : public ConfigurationSpaceRefBase {
public:
  // Inherit the constructors
  using ConfigurationSpaceRefBase::ConfigurationSpaceRefBase;

  void updateRef(const Robot& robot, const GridInfo& grid_info,
                 Eigen::VectorXd& q_ref) const override {
    PYBIND11_OVERRIDE_PURE(void, ConfigurationSpaceRefBase, 
                           updateRef, 
                           robot, grid_info, q_ref);
  }

  bool isActive(const GridInfo& grid_info) const override {
    PYBIND11_OVERRIDE_PURE(bool, ConfigurationSpaceRefBase, 
                           isActive, grid_info);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(configuration_space_ref_base, m) {
  py::class_<ConfigurationSpaceRefBase, 
             PyConfigurationSpaceRefBase,
             std::shared_ptr<ConfigurationSpaceRefBase>>(m, "ConfigurationSpaceRefBase")
    .def(py::init<>())
    .def("updateRef", &ConfigurationSpaceRefBase::updateRef,
          py::arg("robot"), py::arg("grid_info"), py::arg("q_ref"))
    .def("isActive", &ConfigurationSpaceRefBase::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
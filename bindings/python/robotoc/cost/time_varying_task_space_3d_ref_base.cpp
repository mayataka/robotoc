#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"


namespace robotoc {
namespace python {

class PyTimeVaryingTaskSpace3DRefBase : public TimeVaryingTaskSpace3DRefBase {
public:
  // Inherit the constructors
  using TimeVaryingTaskSpace3DRefBase::TimeVaryingTaskSpace3DRefBase;

  void update_x3d_ref(const GridInfo& grid_info, 
                      Eigen::VectorXd& x3d_ref) const override {
    PYBIND11_OVERRIDE_PURE(void, TimeVaryingTaskSpace3DRefBase, 
                           update_x3d_ref, 
                           grid_info, x3d_ref);
  }

  bool isActive(const GridInfo& grid_info) const override {
    PYBIND11_OVERRIDE_PURE(bool, TimeVaryingTaskSpace3DRefBase, 
                           isActive, grid_info);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(time_varying_task_space_3d_ref_base, m) {
  py::class_<TimeVaryingTaskSpace3DRefBase, 
             PyTimeVaryingTaskSpace3DRefBase,
             std::shared_ptr<TimeVaryingTaskSpace3DRefBase>>(m, "TimeVaryingTaskSpace3DRefBase")
    .def(py::init<>())
    .def("update_x3d_ref", &TimeVaryingTaskSpace3DRefBase::update_x3d_ref,
          py::arg("grid_info"), py::arg("x3d_ref"))
    .def("isActive", &TimeVaryingTaskSpace3DRefBase::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/task_space_3d_ref_base.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

class PyTaskSpace3DRefBase : public TaskSpace3DRefBase {
public:
  // Inherit the constructors
  using TaskSpace3DRefBase::TaskSpace3DRefBase;

  std::shared_ptr<TaskSpace3DRefBase> clone() const override {
    PYBIND11_OVERRIDE_PURE(std::shared_ptr<TaskSpace3DRefBase>, TaskSpace3DRefBase, 
                           clone, );
  }

  void updateRef(const GridInfo& grid_info, 
                 Eigen::VectorXd& ref) const override {
    PYBIND11_OVERRIDE_PURE(void, TaskSpace3DRefBase, 
                           updateRef, 
                           grid_info, ref);
  }

  bool isActive(const GridInfo& grid_info) const override {
    PYBIND11_OVERRIDE_PURE(bool, TaskSpace3DRefBase, 
                           isActive, 
                           grid_info);
  }
};

PYBIND11_MODULE(task_space_3d_ref_base, m) {
  py::class_<TaskSpace3DRefBase, 
             PyTaskSpace3DRefBase, 
             std::shared_ptr<TaskSpace3DRefBase>>(m, "TaskSpace3DRefBase")
    .def(py::init<>())
    .def("clone", &TaskSpace3DRefBase::clone)
    .def("updateRef", &TaskSpace3DRefBase::updateRef,
          py::arg("grid_info"), py::arg("ref"))
    .def("isActive", &TaskSpace3DRefBase::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
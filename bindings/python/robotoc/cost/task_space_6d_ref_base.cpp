#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/task_space_6d_ref_base.hpp"
#include "robotoc/robot/se3.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

class PyTaskSpace6DRefBase : public TaskSpace6DRefBase {
public:
  // Inherit the constructors
  using TaskSpace6DRefBase::TaskSpace6DRefBase;

  std::shared_ptr<TaskSpace6DRefBase> clone() const override {
    PYBIND11_OVERRIDE_PURE(std::shared_ptr<TaskSpace6DRefBase>, TaskSpace6DRefBase, 
                           clone, );
  }

  void updateRef(const GridInfo& grid_info, 
                 SE3& ref) const override {
    PYBIND11_OVERRIDE_PURE(void, TaskSpace6DRefBase, 
                           updateRef, 
                           grid_info, ref);
  }

  bool isActive(const GridInfo& grid_info) const override {
    PYBIND11_OVERRIDE_PURE(bool, TaskSpace6DRefBase, 
                           isActive, 
                           grid_info);
  }
};

PYBIND11_MODULE(task_space_6d_ref_base, m) {
  py::class_<TaskSpace6DRefBase, 
             PyTaskSpace6DRefBase, 
             std::shared_ptr<TaskSpace6DRefBase>>(m, "TaskSpace6DRefBase")
    .def(py::init<>())
    .def("clone", &TaskSpace6DRefBase::clone)
    .def("updateRef", &TaskSpace6DRefBase::updateRef,
          py::arg("grid_info"), py::arg("ref"))
    .def("isActive", &TaskSpace6DRefBase::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
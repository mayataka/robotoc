#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/com_ref_base.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

class PyCoMRefBase : public CoMRefBase {
public:
  // Inherit the constructors
  using CoMRefBase::CoMRefBase;

  void updateRef(const GridInfo& grid_info, 
                 Eigen::VectorXd& ref) const override {
    PYBIND11_OVERRIDE_PURE(void, CoMRefBase, 
                           updateRef, 
                           grid_info, ref);
  }

  bool isActive(const GridInfo& grid_info) const override {
    PYBIND11_OVERRIDE_PURE(bool, CoMRefBase, 
                           isActive, 
                           grid_info);
  }
};

PYBIND11_MODULE(com_ref_base, m) {
  py::class_<CoMRefBase, 
             PyCoMRefBase, 
             std::shared_ptr<CoMRefBase>>(m, "CoMRefBase")
    .def(py::init<>())
    .def("updateRef", &CoMRefBase::updateRef,
          py::arg("grid_info"), py::arg("ref"))
    .def("isActive", &CoMRefBase::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc
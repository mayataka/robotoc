#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/se3.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

SE3 createSE3FromQuaternion(const Eigen::Vector4d& quat, const Eigen::Vector3d& trans) {
  return SE3(Eigen::Quaterniond(quat.coeff(3), quat.coeff(0), quat.coeff(1), quat.coeff(2)), trans);
}

PYBIND11_MODULE(se3, m) {
  py::class_<SE3>(m, "SE3")
    .def(py::init<const Eigen::Matrix3d&, const Eigen::Vector3d&>(),
          py::arg("R"), py::arg("trans"))
    .def(py::init(&createSE3FromQuaternion),
          py::arg("quat_x_y_z_w"), py::arg("trans"))
    .def(py::init<const Eigen::Matrix4d&>(),
          py::arg("T"))
    .def(py::init<>())
    .def("clone", [](const SE3& self) {
       auto other = self;
       return other;
     })
    .def_property("R", [](const SE3& self) { return self.rotation(); },
                  [](SE3& self, const Eigen::Matrix3d& R) { self.rotation() = R; },
                  py::return_value_policy::copy)
    .def_property("trans", [](const SE3& self) { return self.translation(); },
                  [](SE3& self, const Eigen::Vector3d& trans) { self.translation() = trans; },
                  py::return_value_policy::copy)
    .def("__str__", [](const SE3& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc
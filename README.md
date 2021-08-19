## idocp - Inverse Dynamics based Optimal Control Problem solver for rigid body systems 

[![build](https://github.com/mayataka/idocp/workflows/build/badge.svg?branch=master)](https://github.com/mayataka/idocp/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/mayataka/idocp/branch/master/graph/badge.svg?token=UOWOF0XO51)](https://codecov.io/gh/mayataka/idocp)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/running_yoko.gif" width="600">

## Features for efficient optimal control for rigid body systems with contacts
- Direct multiple-shooting method based on inverse dynamics.
- Riccati recursion / Parallel Newton's method (ParNMPC) for solving the KKT systems.
- Efficient constraint handling method for pure-state equality constraints in Riccati recursion.
- Primal-dual interior point method for inequality constraints.
- Very fast computation of rigid body dynamics and its sensitivities thanks to [Pinocchio](https://github.com/stack-of-tasks/pinocchio).

## Requirements
- Ubuntu 18.04 or 20.04
- gcc (at least C++11 is required), CMake (at least version 3.1)
- Eigen3, [Pinocchio](https://stack-of-tasks.github.io/pinocchio/download.html)  , 
- Python3, NumPy (for Python binding)
- [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git), [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer) (optional to visualize the solution trajectory) 

## Installation 
1. Install the latest stable version of Eigen3 by 

```
sudo apt install libeigen3-dev
```

2. Install the latest stable version of Pinocchio by following the [instruction](https://stack-of-tasks.github.io/pinocchio/download.html)
3. Clone this repository and change directory as

```
git clone https://github.com/mayataka/idocp
cd idocp
```

4. Build and install `idocp` as

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release 
make install -j$(nproc)
```

5. If you want to visualize the solution trajectory of the OCP, first install [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git) and [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer), e.g., by
```
sudo apt update && sudo apt install robotpkg-py38-qt5-gepetto-viewer-corba -y
git clone https://github.com/stack-of-tasks/pinocchio-gepetto-viewer.git --recursive && cd pinocchio-gepetto-viewer
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make install
```
and change the CMake configuration of `idocp` as 
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_VIEWER=ON
```


6. If you do not want to install Python binding, change the CMake configuration as
```
cmake .. -DBUILD_PYTHON_INTERFACE=OFF
```


## Usage
### C++:
You can link your exectables to `idocp` by writing `CMakeLists.txt`, e.g., as
```
find_package(idocp REQUIRED)

add_executable(
    YOUR_EXECTABLE
    YOUR_EXECTABLE.cpp
)
target_link_libraries(
    YOUR_EXECTABLE
    PRIVATE
    idocp::idocp
)
target_include_directories(
    YOUR_EXECTABLE
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)
```

### Python:
Suppose that thePython version is 3.8. The Python binding will be then installed at `IDOCP_INSTALL_DIR/lib/python3.8/site-packages` where `IDOCP_INSTALL_DIR` is the install directory of `idocp` configured in CMake (e.g., by `-DCMAKE_INSTALL_PREFIX`).
To use the installed Python library, it is convenient to set the environment variable as

```
export PYTHONPATH=IDOCP_INSTALL_DIR/lib/python3.8/site-packages:$PYTHONPATH 
```

e.g., in ~/.bashrc. Note that if you use another Python version than `python3.8`, please adapt it.

## Solvers
The following three solvers are provided:
- `OCPSolver` : Solves the OCP for rigid-body systems (possibly with contacts) by using Riccati recursion.
- `UnconstrOCPSolver` : Solves the OCP for "unconstrained" rigid-body systems by using Riccati recursion.
- `UnconstrParNMPCSolver` : Solves the OCP for "unconstrained" rigid-body systems by using ParNMPC algorithm.

where "unconstrained" rigid-body systems are systems without any contacts or a floating-base.


## Examples
Examples are found in `examples` directory.
The following animations are the solution trajectory of the `UnconstrOCPSolver` for a robot manipulator iiwa14.

- Configuration-space and task-space optimal control (`iiwa14/config_space_ocp.cpp`, `iiwa14/task_space_ocp.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/config_ocp.gif" width="135"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/task_ocp.gif" width="135">


The following animations are the solution trajectory of the `OCPSolver` for a quadruped ANYmal (yellow arrows denote contact forces and blue polyhedrons denote linearized friction cone constraints).

- Walking, trotting gaits (`anymal/walking.cpp`, `anymal/trotting.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/walking.gif" width="250"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/trotting.gif" width="250">

- Pacing, bounding, jumping gaits (`anymal/pacing.cpp`, `anymal/bounding.cpp`, `anymal/jumping.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/pacing.gif" width="250"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/bounding.gif" width="250"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/jumping.gif" width="250">

- Running gait (`anymal/running.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/running.gif" width="520">


## Citing idocp

Citing `UnconstrOCPSolver` and `UnconstrParNMPCSolver`:
```
@inproceedings{katayama2021idocp,
  title={Efficient solution method based on inverse dynamics for optimal control problems of rigid body systems},
  author={S. Katayama and T. Ohtsuka},
  booktitle={{IEEE International Conference on Robotics and Automation (ICRA)}},
  year={2021}}
```

## Related publications
- S. Katayama and T. Ohtsuka, Lifted contact dynamics for efficient direct optimal control of rigid body systems with contacts, https://arxiv.org/abs/2108.01781, 2021
- S. Katayama and T. Ohtsuka, Efficient Riccati recursion for optimal control problems with pure-state equality constraints, https://arxiv.org/abs/2102.09731, 2021
- S. Katayama and T. Ohtsuka, Efficient solution method based on inverse dynamics for optimal control problems of rigid body systems, IEEE International Conference on Robotics and Automation (ICRA), 2021
- H. Deng and T. Ohtsuka, A parallel Newton-type method for nonlinear model predictive control, Automatica, Vol. 109, pp. 108560, 2019
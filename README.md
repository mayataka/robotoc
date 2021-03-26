## idocp - Inverse Dynamics based Optimal Control Problem solver for rigid body systems 

[![build](https://github.com/mayataka/idocp/workflows/build/badge.svg?branch=master)](https://github.com/mayataka/idocp/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/mayataka/idocp/branch/master/graph/badge.svg?token=UOWOF0XO51)](https://codecov.io/gh/mayataka/idocp)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/anymal_running_minimal.gif" width="600">

## Features for efficient optimal control for rigid body systems with contacts
- Direct multiple-shooting formulation based on inverse dynamics with novel condensing method dedicated for rigid body systems with contacts.
- Riccati recursion / Parallel Newton's method (ParNMPC) for solving the condensed KKT system.
- Efficient constraint handling method for pure-state equality constraints in Riccati recursion.
- Primal-dual interior point method for inequality constraints.
- Filter line-search method.
- Very fast computation of rigid body dynamics and its sensitivities thanks to [pinocchio](https://github.com/stack-of-tasks/pinocchio).

## Requirements
- Ubuntu 18.04 or 20.04
- gcc (at least C++11 is required), CMake (at least version 3.1)
- [Eigen3](https://stack-of-tasks.github.io/pinocchio/download.html)  
- [pinocchio](https://github.com/stack-of-tasks/pinocchio) (instruction for installation is found [here](https://stack-of-tasks.github.io/pinocchio/download.html))
- [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer), [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git) (Optional to visualize the solution trajectory) 

## Installation 
1. Install the latest stable version of Eigen3 by 

```
sudo apt install libeigen3-dev
```

2. Install the latest stable version of pinocchio by following the [instruction](https://stack-of-tasks.github.io/pinocchio/download.html)
3. Clone this repository and change directory as

```
git clone https://github.com/mayataka/idocp
cd idocp
```

5. Build and install `idocp` as

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DTESTING=False
make -j$(nproc)
make install
```

6. If you want to visualize the solution trajectory of the OCP, first install [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git) and [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer), e.g., by
```
sudo apt update && sudo apt install robotpkg-py38-qt5-gepetto-viewer-corba -y
git clone https://github.com/stack-of-tasks/pinocchio-gepetto-viewer.git --recursive && cd pinocchio-gepetto-viewer
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make install
```
and change the CMake configuration of `idocp` as 
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DTESTING=False -DBUILD_VIEWER=True
```

## Usage
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

The following four solvers are provided:
- `idocp::UnOCPSolver` : Solves the OCP for "unconstrained" rigid-body systems by using Riccati recursion.
- `idocp::UnParNMPCSolver` : Solves the OCP for "unconstrained" rigid-body systems by using ParNMPC algorithm.
- `idocp::OCPSolver` : Solves the OCP for rigid-body systems by using Riccati recursion.
- `idocp::ParNMPCSolver` : Solves the OCP for rigid-body systems by using ParNMPC algorithm.

where "unconstrained" rigid-body systems are systems without any contacts or a floating-base.


## Examples
Examples are found in `examples` directory.
The following animations are the solution trajectory of the `idocp::UnOCPSolver` for a robot manipulator iiwa14.

- Configuration-space optimal control (`iiwa14/config_space_ocp.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/iiwa14_config_ocp.gif" width="170">

- Task-space optimal control (`iiwa14/task_space_ocp.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/iiwa14_task_ocp.gif" width="170">

The following animations are the solution trajectory of the `idocp::OCPSolver` for a quadruped ANYmal.


- Trotting (`anymal/anymal_trotting.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/anymal_trotting.gif" width="300">

- Jumping (`anymal/anymal_jumping.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/anymal_jumping.gif" width="300">

- Running (`anymal/anymal_running.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/anymal_running.gif" width="600">


## MPC Simulation 
Simulation of the MPC of systems with rigid contacts are shown in [idocp-sim](https://github.com/mayataka/idocp-sim).


## Citing idocp

Citing `idocp::UnOCPSolver` and `idocp::UnParNMPCSolver`:
```
@inproceedings{katayama2021idocp,
  title={Efficient solution method based on inverse dynamics for optimal control problems of rigid body systems},
  author={S. Katayama and T. Ohtsuka},
  booktitle={{IEEE International Conference on Robotics and Automation (ICRA)}},
  year={2021}}
```

## Related publications
- S. Katayama and T. Ohtsuka, Efficient Riccati recursion for optimal control problems with pure-state equality constraints, https://arxiv.org/abs/2102.09731, 2021
- S. Katayama and T. Ohtsuka, Efficient solution method based on inverse dynamics for optimal control problems of rigid body systems, IEEE International Conference on Robotics and Automation (ICRA), 2021
- H. Deng and T. Ohtsuka, A parallel Newton-type method for nonlinear model predictive control, Automatica, Vol. 109, pp. 108560, 2019
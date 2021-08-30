## idocp - Inverse Dynamics based Optimal Control Problem solver for rigid body systems 

[![build](https://github.com/mayataka/idocp/workflows/build/badge.svg?branch=master)](https://github.com/mayataka/idocp/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/mayataka/idocp/branch/master/graph/badge.svg?token=UOWOF0XO51)](https://codecov.io/gh/mayataka/idocp)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/running_yoko.gif" width="530">

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
- [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git) and/or [meshcat-python](https://github.com/rdeits/meshcat-python) (optional to visualize the solution trajectory in Python) 
- [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer) (optional to visualize the solution trajectory in C++) 
- [PyBullet](https://pybullet.org/wordpress/) (optional to simulate MPC examples)

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
NOTE: if you want to maximize the performance, use CMake option
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DOPTIMIZE_FOR_NATIVE=ON
```

5. If you want to visualize the solution trajectory with Python, you have to install [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git) by
```
sudo apt update && sudo apt install robotpkg-py38-qt5-gepetto-viewer-corba -y
```
and/or [meshcat-python](https://github.com/rdeits/meshcat-python) by
```
pip install meshcat
```

6. If you want to visualize the solution trajectory with C++, in addition to [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git), you have to install [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer), e.g., by
```
git clone https://github.com/stack-of-tasks/pinocchio-gepetto-viewer.git --recursive && cd pinocchio-gepetto-viewer
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make install
```
and add the CMake option for `idocp` as 
```
cmake .. -DBUILD_VIEWER=ON
```

7. If you do not want to install the Python bindings, change the CMake configuration as
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
Suppose that the Python version is 3.8. The Python bindings will then be installed at `IDOCP_INSTALL_DIR/lib/python3.8/site-packages` where `IDOCP_INSTALL_DIR` is the install directory of `idocp` configured in CMake (e.g., by `-DCMAKE_INSTALL_PREFIX`).
To use the installed Python library, it is convenient to set the environment variable as

```
export PYTHONPATH=IDOCP_INSTALL_DIR/lib/python3.8/site-packages:$PYTHONPATH 
```

e.g., in `~/.bashrc`. Note that if you use another Python version than `python3.8`, please adapt it.

## Solvers 
The following three solvers are provided:
- `OCPSolver` : Solves the OCP for rigid-body systems (possibly with contacts) by using Riccati recursion.
- `UnconstrOCPSolver` : Solves the OCP for "unconstrained" rigid-body systems by using Riccati recursion.
- `UnconstrParNMPCSolver` : Solves the OCP for "unconstrained" rigid-body systems by using ParNMPC algorithm.

where "unconstrained" rigid-body systems are systems without any contacts or a floating-base.

## Optimal control examples
Examples of these solvers are found in `examples` directory.
The following animations are the solution trajectory of the `UnconstrOCPSolver` for a robot manipulator iiwa14.

- Configuration-space and task-space optimal control (`iiwa14/config_space_ocp.cpp`, `iiwa14/task_space_ocp.cpp`, or `iiwa14/python/config_space_ocp.py`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/config_ocp.gif" width="115"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/task_ocp.gif" width="115">


The following animations are the solution trajectory of the `OCPSolver` for a quadruped ANYmal (yellow arrows denote contact forces and blue polyhedrons denote linearized friction cone constraints).

- Walking, trotting gaits (`anymal/walking.cpp`, `anymal/trotting.cpp`, or `anymal/python/walking.py`, `anymal/python/trotting.py`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/walking.gif" width="215"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/trotting.gif" width="215">

- Pacing, bounding, jumping gaits (`anymal/pacing.cpp`, `anymal/bounding.cpp`, `anymal/jumping.cpp`, or `anymal/python/pacing.py`, `anymal/python/bounding.py`, `anymal/python/jumping.py`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/pacing.gif" width="215"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/bounding.gif" width="215"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/jumping.gif" width="215">

- Running gait (`anymal/running.cpp`)

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/running.gif" width="500">


## Whole-body MPC examples
The following two example implementations of whole-body MPC are provided:
- `MPCQuadrupedalWalking` : MPC with `OCPSolver` for the walking gait of quadrupedal robots.
- `MPCQuadrupedalTrotting` : MPC with `OCPSolver` for the trotting gait of quadrupedal robots.

You can run the simulations of these MPC with `anymal/mpc/walking.py` and `anymal/mpc/trotting.py` (you need to install [PyBullet](https://pybullet.org/wordpress/), e.g., by `pip install pybullet`).

<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/mpc_walking.gif" width="300"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/idocp/images/mpc_trotting.gif" width="300">


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
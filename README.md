## robotoc - efficient ROBOT Optimal Control solvers  

[![build](https://github.com/mayataka/robotoc/workflows/build/badge.svg?branch=master)](https://github.com/mayataka/robotoc/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/mayataka/robotoc/branch/master/graph/badge.svg?token=UOWOF0XO51)](https://codecov.io/gh/mayataka/robotoc)

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/forward_c.gif" width="140"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/side_c.gif" width="140"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/curve_c.gif" width="140"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/rotation_c.gif" width="140"> 

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/longitudinal_mpc_x05.gif" width="140"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/lateral_mpc_x05.gif" width="140"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/back_mpc_x05.gif" width="140"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/rotational_mpc_x05.gif" width="140"> 

## Features for efficient optimal control of robotic systems
- Direct multiple-shooting method based on the **lifted contact dynamics** / **inverse dynamics**.
- Riccati recursion algorithm for **switching time optimization (STO)** problems and efficient **pure-state equality constraint handling**.
- **Primal-dual interior point method** for inequality constraints.
- Very fast computation of rigid body dynamics and its sensitivities thanks to [Pinocchio](https://github.com/stack-of-tasks/pinocchio).

## Requirements
- Ubuntu 22.04, 20.04, and 18.04 (possibly Mac OSX, but not well tested)
- gcc (at least C++11 is required), CMake (at least version 3.11)
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
3. Clone this repository and change the directory as

```
git clone https://github.com/mayataka/robotoc
cd robotoc 
```

4. Build and install `robotoc` as

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
and add the CMake option for `robotoc` as 
```
cmake .. -DBUILD_VIEWER=ON
```

7. If you do not want to install the Python bindings, change the CMake configuration as
```
cmake .. -DBUILD_PYTHON_INTERFACE=OFF
```

8. In OSX, explicitly set g++ as the complier. First, find the path of g++ as 
```
ls -l /usr/local/bin | grep g++
```
Then set the full path in the cmake as 
```
cmake .. -DCMAKE_CXX_COMPILER=FULL_PATH_TO_GPLUSPLUS
```


## Usage
### C++:
You can link your executables to `robotoc` by writing `CMakeLists.txt`, e.g., as
```
find_package(robotoc REQUIRED)

add_executable(
    YOUR_EXECTABLE
    YOUR_EXECTABLE.cpp
)
target_link_libraries(
    YOUR_EXECTABLE
    PRIVATE
    robotoc::robotoc
)
```

### Python:
Suppose that the Python version is 3.8. The Python bindings will then be installed at `ROBOTOC_INSTALL_DIR/lib/python3.8/site-packages` where `ROBOTOC_INSTALL_DIR` is the install directory of `robotoc` configured in CMake (e.g., by `-DCMAKE_INSTALL_PREFIX`).
To use the installed Python library, it is convenient to set the environment variable as

```
export PYTHONPATH=ROBOTOC_INSTALL_DIR/lib/python3.8/site-packages:$PYTHONPATH 
```

e.g., in `~/.bashrc`. Note that if you use another Python version than `python3.8`, please adapt it.

## Solvers 
The following three solvers are provided:
- `OCPSolver` : Solves the OCP for rigid-body systems (possibly with contacts) by using Riccati recursion. Can optimize the switching times and the trajectories simultaneously. 
- `UnconstrOCPSolver` : Solves the OCP for "unconstrained" rigid-body systems by using Riccati recursion.
- `UnconstrParNMPCSolver` : Solves the OCP for "unconstrained" rigid-body systems by using ParNMPC algorithm.

where "unconstrained" rigid-body systems are systems without any contacts or a floating-base.

## Examples
Examples of these solvers are found in `examples` directory.
Further explanations are found at https://mayataka.github.io/robotoc/page_examples.html.

### Switching time optimization (STO) examples
- `OCPSolver` can solve the switching time optimization (STO) problem, which optimizes the trajectory and the contact timings simultaneously.
- The following videos display the solution trajectory of the STO problems (`anymal/python/jumping_sto.py` and `icub/python/jumping_sto.py`).

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/jumping_sto.gif" width="240"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/icub.gif" width="230">

### Whole-body MPC examples
- The following example implementations of whole-body MPC are provided:
  - `MPCCrawl` : MPC with `OCPSolver` for the crawl gait of quadrupedal robots.
  - `MPCTrot` : MPC with `OCPSolver` for the trot gait of quadrupedal robots.
  - `MPCPace` : MPC with `OCPSolver` for the pace gait of quadrupedal robots.
  - `MPCFlyingTrot` : MPC with `OCPSolver` for the flying trot gait of quadrupedal robots.
  - `MPCJump` : MPC with `OCPSolver` for the jump motion of quadrupedal or bipedal robots.
  - `MPCBipedWalk` : MPC with `OCPSolver` for the walking motion of bipedal robots.
- You can find the simulations of these MPC at `a1/mpc`, `anymal/mpc`, and `icub/mpc`.
- You need [PyBullet](https://pybullet.org/wordpress/) for the MPC simulations. (easy to install, e.g., by `pip install pybullet`)
- The following videos display the simulation results of quadrupedal walking on rough terrain (`a1/mpc/trot_terrain.py`) and bipedal walking (`icub/mpc/walk.py`).

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/a1_trot_terrain.gif" width="270"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/icub_walk.gif" width="250">


## Documentation
More detailed documentation is available at https://mayataka.github.io/robotoc/.

## Tutorial
Our tutorials on robotoc is available at [robotoc_tutorial](https://github.com/mayataka/robotoc_tutorial).

## Citing robotoc
- Citing the STO algorithm of `OCPSolver`:
```
@misc{katayama2021sto,
  title={Structure-exploiting {N}ewton-type method for optimal control of switched systems}, 
  author={Sotaro Katayama and Toshiyuki Ohtsuka},
  url={arXiv:2112.07232},
  eprint={2112.07232},
  archivePrefix={arXiv}
  year={2021}}
```

- Citing `OCPSolver` without the switching time optimization (STO):
```
@inproceedings{katayama2022lifted,
  title={Lifted contact dynamics for efficient optimal control of rigid body systems with contacts}, 
  author={Sotaro Katayama and Toshiyuki Ohtsuka},
  booktitle={{2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS 2022) (to appear)}},
  year={2022}}
```

- Citing `UnconstrOCPSolver` and `UnconstrParNMPCSolver` (the repository name was `idocp` in this paper (https://github.com/mayataka/idocp)):
```
@inproceedings{katayama2021idocp,
  title={Efficient solution method based on inverse dynamics for optimal control problems of rigid body systems},
  author={Sotaro Katayama and Toshiyuki Ohtsuka},
  booktitle={{2021 IEEE International Conference on Robotics and Automation (ICRA)}},
  pages={2070--2076},
  year={2021}}
```

## Related publications
- S. Katayama and T. Ohtsuka, "Whole-body model predictive control with rigid contacts via online switching time optimization," 2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS 2022) (to appear), https://arxiv.org/abs/2203.00997, 2022
- S. Katayama and T. Ohtsuka, "Structure-exploiting Newton-type method for optimal control of switched systems," https://arxiv.org/abs/2112.07232, 2021
- S. Katayama and T. Ohtsuka, Lifted contact dynamics for efficient optimal control of rigid body systems with contacts, 2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (to appear), https://arxiv.org/abs/2108.01781, 2022
- S. Katayama and T. Ohtsuka, Efficient Riccati recursion for optimal control problems with pure-state equality constraints, 2022 American Control Conference (ACC), pp. 3579-3586, 2022
- S. Katayama and T. Ohtsuka, Efficient solution method based on inverse dynamics for optimal control problems of rigid body systems, 2021 IEEE International Conference on Robotics and Automation (ICRA), pp.2070-2076, 2021
- H. Deng and T. Ohtsuka, A parallel Newton-type method for nonlinear model predictive control, Automatica, Vol. 109, pp. 108560, 2019

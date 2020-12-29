## idocp - Inverse Dynamics based Optimal Control Problem solver for rigid body systems 

[![build](https://github.com/mayataka/idocp/workflows/build/badge.svg?branch=master)](https://github.com/mayataka/idocp/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/mayataka/idocp/branch/master/graph/badge.svg?token=UOWOF0XO51)](https://codecov.io/gh/mayataka/idocp)

## Features for efficient optimal control 
- Solves the optimal control problem (OCP) for rigid body systems with contacts based on inverse dynamics and event-driven scheme.
- Riccati recursion / Parallel Newton's method (ParNMPC) for solving the KKT system.
- Primal-dual interior point method for inequality constraints.
- Filter line-search method.
- Very fast computation of rigid body dynamics and its sensitivities thanks to [pinocchio](https://github.com/stack-of-tasks/pinocchio).

## Requirements
- Ubuntu 18.04 or 20.04
- gcc, CMake
- [Eigen3](https://stack-of-tasks.github.io/pinocchio/download.html)  
- [pinocchio](https://github.com/stack-of-tasks/pinocchio) (instruction for installation is found [here](https://stack-of-tasks.github.io/pinocchio/download.html))
- [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer), [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git) (Optional to visualization of the solution trajectory) 

## Installation 
1. Install latest stable version of Eigen3 by 

```
sudo apt install libeigen3-dev
```

2. Install latest stable version of pinocchio by following the [instruction](https://stack-of-tasks.github.io/pinocchio/download.html)
3. Clone this repository and change directory as

```
git clone https://github.com/mayataka/idocp
cd idocp
```

5. Build and install `idocp` as

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DTESTING=False
make -j$(nproc)
sudo make install
```

6. If you want to visualize the solution trajectory of the OCP, first install [gepetto-viewer-corba](https://github.com/Gepetto/gepetto-viewer-corba.git) and [pinocchio-gepetto-viewer](https://github.com/stack-of-tasks/pinocchio-gepetto-viewer), e.g., by
```
sudo apt update && sudo apt install robotpkg-py27-qt4-gepetto-viewer-corba 
git clone https://github.com/stack-of-tasks/pinocchio-gepetto-viewer.git --recursive
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
sudo ln -s /usr/lib/x86_64-linux-gnu/libSM.so.6 /usr/lib/x86_64-linux-gnu/libSM.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libICE.so.6 /usr/lib/x86_64-linux-gnu/libICE.so
sudo make install
```
and configure the build of `idocp` as 
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DTESTING=False -DBUILD_VIEWER=True
```

## Usage
You can link your exectables to `idocp` by writing `CMakeLists.txt` as
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

## MPC Simulation 
Simulation of the MPC of systems with rigid contacts are provided in [idocp-sim](https://github.com/mayataka/idocp-sim).


## Related publications

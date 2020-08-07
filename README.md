## idocp - Inverse Dynamics based Optimal Control Problem solver for rigid body systems 

[![Build Status](https://travis-ci.com/mayataka/IDOCP.svg?token=fusqwLK1c8Q529AAxFz6&branch=master)](https://travis-ci.com/mayataka/IDOCP)
[![codecov](https://codecov.io/gh/mayataka/IDOCP/branch/master/graph/badge.svg?token=UOWOF0XO51)](https://codecov.io/gh/mayataka/IDOCP)

## Features for efficient optimal control 
- Solves the optimal control problem for rigid body systems based on inverse dynamics.
- Sparsity-exploiting Riccati recursion for computing the Newton's direction.
- Primal-dual interior point method for inequality constraints.
- Filter line-search method.
- Very fast computation of rigid body dynamics and its sensitivities thanks to [pinocchio](https://github.com/stack-of-tasks/pinocchio).

## Requirements
- Ubuntu 18.04 
- gcc
- [pinocchio](https://github.com/stack-of-tasks/pinocchio) (instruction for installation is found [here](https://stack-of-tasks.github.io/pinocchio/download.html))
- [Eigen3](https://stack-of-tasks.github.io/pinocchio/download.html)  

## Installation 
1. Install Eigen3 by `sudo apt install libeigen3-dev`.

2. Install pinocchio by following the [instruction](https://stack-of-tasks.github.io/pinocchio/download.html).
3. Clone this repository.
4. Change current directory by `cd idocp`.
5. Build and install idocp as

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
sudo make install
```


## Related publications


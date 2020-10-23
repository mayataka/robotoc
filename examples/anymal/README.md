## Simulation of MPC

Class `idocp::QuadrupedSimulator` provides quadruped simulation for MPC utilizing RaiSim.
You can build exectables as (assume ${RAISIM_LOCAL_BUILD_DIR} is ~/raisim_build)
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DSIM_MPC=ON -DCMAKE_PREFIX_PATH=~/raisim_build
make
```
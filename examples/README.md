## Examples
Examples of these solvers are found in `examples` directory.
Further explanations are found at https://mayataka.github.io/robotoc/page_examples.html.

### Optimal control example with fixed contact timings

- Configuration-space and task-space optimal control for a robot manipulator iiwa14 using `UnconstrOCPSolver` (`iiwa14/config_space_ocp.cpp`, `iiwa14/task_space_ocp.cpp`, or `iiwa14/python/config_space_ocp.py`):

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/config_ocp.gif" width="100"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/task_ocp.gif" width="100">


- Crawl, trot, pacing, and bounding gaits of quadruped ANYmal using `OCPSolver` (yellow arrows denote contact forces and blue polyhedrons denote linearized friction cone constraints) with fixed contact timings (e.g., `anymal/crawl.cpp` or `anymal/python/crawl.py`):

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/walking.gif" width="180"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/trotting.gif" width="180"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/pacing.gif" width="180"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/bounding.gif" width="180">

### Switching time optimization (STO) examples
-  `OCPSolver` for the switching time optimization (STO) problem, which optimizes the trajectory and the contact timings simultaneously (`anymal/python/jumping_sto.py` and `icub/python/jumping_sto.py`):

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/jumping_sto.gif" width="250"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/icub.gif" width="230">

### Whole-body MPC examples
- The following example implementations of whole-body MPC are provided:
  - `MPCCrawl` : MPC with `OCPSolver` for the crawl gait of quadrupedal robots.
  - `MPCTrot` : MPC with `OCPSolver` for the trot gait of quadrupedal robots.
  - `MPCPace` : MPC with `OCPSolver` for the pace gait of quadrupedal robots.
  - `MPCFlyingTrot` : MPC with `OCPSolver` for the flying trot gait of quadrupedal robots.
  - `MPCJump` : MPC with `OCPSolver` for the jump motion of quadrupedal or bipedal robots.
  - `MPCBipedWalk` : MPC with `OCPSolver` for the walking motion of bipedal robots.
- You can find the simulations of these MPC at `a1/mpc`, `anymal/mpc`, and `icub/mpc` (you need to install [PyBullet](https://pybullet.org/wordpress/), e.g., by `pip install pybullet`):

<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/a1_trotting.gif" width="250"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/anymal_crawling.gif" width="250"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/robotoc/images/icub_walking.gif" width="240">
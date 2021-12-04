import time


def cpu_time(ocp_solver, t, q, v, num_iteration):
    start_clock = time.time()
    for i in range(num_iteration):
        ocp_solver.update_solution(t, q, v)
    end_clock = time.time()
    print('---------- OCP benchmark : CPU time ----------')
    print('total CPU time: {:.6g}'.format(1e03*(end_clock-start_clock)) + '[ms]')
    print('CPU time per update: {:.6g}'.format(1e03*(end_clock-start_clock)/num_iteration) + '[ms]')
    print('-----------------------------------')


def convergence(ocp_solver, t, q, v, num_iteration, logger=None):
    print('---------- OCP benchmark : Convergence ----------')
    print('Initial KKT error = {:.6g}'.format(ocp_solver.KKT_error(t, q, v)))
    if logger is not None:
        logger.take_log(ocp_solver)
    for i in range(num_iteration):
        ocp_solver.update_solution(t, q, v)
        print('KKT error after iteration ' + str(i+1) + ' = {:.6g}'.format(ocp_solver.KKT_error(t, q, v)))
        if logger is not None:
            logger.take_log(ocp_solver)
    print('-----------------------------------')


def convergence_sto(ocp_solver, t, q, v, num_iteration, dt_tol_mesh, kkt_tol_mesh, logger=None):
    print('---------- OCP benchmark : Convergence ----------')
    print('Initial KKT error = {:.6g}'.format(ocp_solver.KKT_error(t, q, v)))
    if logger is not None:
        logger.take_log(ocp_solver)
    for i in range(num_iteration):
        ocp_solver.update_solution(t, q, v)
        if ocp_solver.get_OCP_discretization().dt_max() > dt_tol_mesh and ocp_solver.KKT_error() < kkt_tol_mesh:
            print('Mesh refinement is carried out!')
            ocp_solver.mesh_refinement(t)
        print('KKT error after iteration ' + str(i+1) + ' = {:.6g}'.format(ocp_solver.KKT_error(t, q, v)))
        if logger is not None:
            logger.take_log(ocp_solver)
    print('-----------------------------------')

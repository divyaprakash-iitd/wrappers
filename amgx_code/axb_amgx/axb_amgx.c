#define MAX_MSG_LEN 4096

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "axb_amgx_headers.h"
#include <amgx_c.h>

//#include "cuda_runtime.h" 
/* CUDA error macro */
#define CUDA_SAFE_CALL(call) do {                                 \
  cudaError_t err = call;                                         \
  if(cudaSuccess != err) {                                        \
    fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
            __FILE__, __LINE__, cudaGetErrorString( err) );       \
    exit(EXIT_FAILURE);                                           \
  } } while (0)

/* standard or dynamically load library */
#ifdef AMGX_DYNAMIC_LOADING
#include "amgx_capi.h"
#else
#include "amgx_c.h"
#endif

/* print error message and exit */
void errAndExit(const char *err)
{
    printf("%s\n", err);
    fflush(stdout);
    exit(1);
}

/* print callback (could be customized) */
void print_callback(const char *msg, int length)
{
    printf("%s", msg);
}

// GLOBAL VARIABLES
//input matrix and rhs/solution
int n;
int bsize_x;
int bsize_y;
int sol_size;
int sol_bsize;
int N; int nnz; int block_dimx; int block_dimy;

//library handles
AMGX_Mode mode;
AMGX_config_handle cfg;
AMGX_resources_handle rsrc;
AMGX_matrix_handle A;
AMGX_vector_handle b, x;
AMGX_solver_handle solver;
//status handling
AMGX_SOLVE_STATUS status;

//int main(int argc, const char **argv)
int initialize_amgx(int *crs_data, double *data, int *col_ind, int *row_ptr, double *rhs, double *sol) 
{
    
    //printf("The size of crs_data : %d\n", sizeof(crs_data));
    //printf("The size of val is: %d\n", sizeof(data));
    //printf("The size of col_ind is: %d\n", sizeof(col_ind));
    //printf("The size of row_ptr is: %d\n", sizeof(row_ptr));
    //printf("The size of rhs is: %d\n", sizeof(rhs));
    //printf("The size of sol is: %d\n", sizeof(sol));
    
    //input matrix and rhs/solution
    n = 0;
    bsize_x = 0;
    bsize_y = 0;
    sol_size = 0;
    sol_bsize = 0;


    /* load the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
    void *lib_handle = NULL;
    //open the library
#ifdef _WIN32
    lib_handle = amgx_libopen("amgxsh.dll");
#else
    lib_handle = amgx_libopen("libamgxsh.so");
#endif

    if (lib_handle == NULL)
    {
        errAndExit("ERROR: can not load the library");
    }

    //load all the routines
    if (amgx_liblink_all(lib_handle) == 0)
    {
        amgx_libclose(lib_handle);
        errAndExit("ERROR: corrupted library loaded\n");
    }

#endif
    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    /* system */
    AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
    AMGX_SAFE_CALL(AMGX_install_signal_handler());


    /* get mode */
    mode = AMGX_mode_dDDI;

    /* Specify configuration */
    AMGX_SAFE_CALL(AMGX_config_create(&cfg,
    "config_version=2,\
    solver=FGMRES,\
    gmres_n_restart=20,\
    max_iters=10000,\
    norm=L2,\
    convergence=RELATIVE_INI_CORE,\
    monitor_residual=1,\
    tolerance=1e-4,\
    preconditioner(amg_solver)=AMG,\
    amg_solver:algorithm=CLASSICAL,\
    amg_solver:max_iters=2,\
    amg_solver:presweeps=1,\
    amg_solver:postsweeps=1,\
    amg_solver:cycle=V,\
    print_solve_stats=0,\
    print_grid_stats=0,\
    obtain_timings=0"));

    ///* Specify configuration */
    //AMGX_SAFE_CALL(AMGX_config_create(&cfg,
    //"config_version=2,\
    //solver(my_solver)=FGMRES,\
    //my_solver:gmres_n_restart=20,\
    //my_solver:max_iters=10000,\
    //my_solver:norm=L2,\
    //my_solver:convergence=RELATIVE_INI_CORE,\
    //my_solver:monitor_residual=1,\
    //my_solver:tolerance=1e-5,\
    //my_solver:preconditioner(amg_solver)=AMG,\
    //amg_solver:algorithm=CLASSICAL,\
    //amg_solver:max_iters=2,\
    //amg_solver:presweeps=1,\
    //amg_solver:postsweeps=1,\
    //amg_solver:cycle=V,\
    //my_solver:print_solve_stats=1,\
    //my_solver:print_grid_stats=1,\
    //my_solver:obtain_timings=1"));

    //AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version=2, solver(amg_solver)=AMG"));
    //AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version=2, solver(amg_solver)=AMG, amg_solver:print_solve_stats=0,amg_solver: monitor_residual=1, amg_solver:max_iters=1000000, amg_solver:tolerance=0.00000000000001, amg_solver:norm=L2, amg_solver:store_res_history=0, amg_solver:obtain_timings=0"));
    //AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version=2, solver(amg_solver)=AMG, amg_solver:print_solve_stats=1,amg_solver: monitor_residual=1, amg_solver:max_iters=1000000, amg_solver:tolerance=0.0000000001, amg_solver:norm=L2, amg_solver:store_res_history=1, amg_solver:obtain_timings=1")); 
    //AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version= 2, use_scalar_norm= 1,solver= BLOCK_JACOBI,print_solve_stats= 1, obtain_timings= 1,monitor_residual= 1, convergence= RELATIVE_INI_CORE, tolerance= 1e-14, norm= L2,max_iters=10000"));

   //AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version= 2, determinism_flag= 1, print_grid_stats= 1, max_uncolored_percentage= 0.15, algorithm= AGGREGATION, obtain_timings= 1, solver= AMG, smoother= MULTICOLOR_GS, print_solve_stats= 1, presweeps= 1, symmetric_GS= 1, selector= SIZE_2, coarsest_sweeps= 2, max_iters= 1000, monitor_residual= 1, postsweeps= 1, max_levels= 1000, matrix_coloring_scheme= MIN_MAX, tolerance= 0.1, norm= L1, cycle= V"));

  //AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version= 2, solver= AMG, smoother= MULTICOLOR_GS, max_iters= 1000, tolerance= 1e-6, norm= L2"));
    
    //AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version=2, convergence=ABSOLUTE, max_iters=1e6, monitor_residual=1, norm=L2, solver(my_solver)=AMG, tolerance=1e-10, print_solve_stats=1, print_grid_stats=1,")); 


    //AMGX_SAFE_CALL(AMGX_config_create(&cfg,"config_version=2,solver=FGMRES,gmres_n_restart=20,max_iters=100,norm=L2,convergence=RELATIVE_INI_CORE,monitor_residual=1,tolerance=1e-4,preconditioner(amg_solver)=AMG,amg_solver:algorithm=CLASSICAL,amg_solver:max_iters=2,amg_solver:presweeps=1,amg_solver:postsweeps=1,amg_solver:cycle=V,print_solve_stats=1,print_grid_stats=1,obtain_timings=1"));

   //AMGX_SAFE_CALL(AMGX_config_create(&cfg,"config_version=2, preconditioner(my_solver)=AMG, my_solver:error_scaling= 0, my_solver:print_grid_stats=1, my_solver:max_uncolored_percentage=0.05, my_solver:algorithm=AGGREGATION,  my_solver:smoother=MULTICOLOR_DILU, my_solver:presweeps=0, my_solver:selector=SIZE_2, my_solver:coarse_solver=DENSE_LU_SOLVER, my_solver:max_iters=1, my_solver:postsweeps=3, my_solver:min_coarse_rows=32, my_solver:relaxation_factor=0.75, my_solver:max_levels=100, my_solver:matrix_coloring_scheme=PARALLEL_GREEDY, my_solver:cycle= V, use_scalar_norm=1, solver=FGMRES, my_solver:print_solve_stats=1, obtain_timings=1, max_iters=100, my_solver:monitor_residual=1, gmres_n_restart=10, convergence=RELATIVE_INI_CORE, tolerance=1e-10, norm=L2"));
    /* create resources, matrix, vector and solver */
    AMGX_resources_create_simple(&rsrc, cfg);
    AMGX_matrix_create(&A, rsrc, mode);
    AMGX_vector_create(&x, rsrc, mode);
    AMGX_vector_create(&b, rsrc, mode);
    AMGX_solver_create(&solver, rsrc, mode, cfg);

    /* read the input system: matrix [and rhs & solution]
       Please refer to AMGX_read_system description in the AMGX_Reference.pdf
       manual for details on how to specify the rhs and the solution inside
       the input file. If these are not specified than rhs=[1,...,1]^T and
       (initial guess) sol=[0,...,0]^T. */

    /* Input your own matrix */

    /*double data[] = {-2, 1, 1, -2, 1, 1, -2, 1, 1, -2 };
    int col_ind[] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    int row_ptr[] = {0, 2, 5, 8, 10};*/

    // int N = 4; int nnz = 10; int block_dimx = 1; int block_dimy = 1;
    
    // int N; int nnz; int block_dimx; int block_dimy;
    N = crs_data[0]; nnz = crs_data[1]; block_dimx = crs_data[2]; block_dimy = crs_data[3];
    
    AMGX_matrix_upload_all(A, N, nnz, block_dimx, block_dimy, row_ptr, col_ind, data, 0);

    
    AMGX_matrix_get_size(A, &n, &bsize_x, &bsize_y);

    /* Input your RHS vector */
    //double rhs[] = {-300, 0, 0, -100};
    //AMGX_pin_memory(rhs);
    AMGX_vector_upload(b, N, 1, rhs);
    
    AMGX_vector_get_size(x, &sol_size, &sol_bsize);
    /* Input your initial guess, x */
    //AMGX_vector_set_zero(x, n, bsize_x);
    AMGX_vector_upload(x, N, 1, sol);


    /* solver setup */
    AMGX_solver_setup(solver, A);
    /* solver solve */
    AMGX_solver_solve(solver, b, x);
    /* example of how to change parameters between non-linear iterations */
    //AMGX_config_add_parameters(&cfg, "config_version=2, default:tolerance=1e-12");
    
    //AMGX_solver_solve(solver, b, x);
    AMGX_solver_get_status(solver, &status);
    /* example of how to print the residual history */
    //int nit;
    //double res;
    //AMGX_solver_get_iterations_number(solver, &nit);
    //for (int i=0; i<nit; i++) {
    //  printf("residual from iteration %d=", i);
    //  for (int j=0; j<bsize_y; j++) {
    //    AMGX_solver_get_iteration_residual(solver, i, j, &res);
    //    printf("%f ", (float)(res));
    //  }
    //  printf("\n");
    //}

    /* example of how to write the linear system to the output */
    AMGX_vector_download(x,sol);
    //AMGX_write_system(A, b, x, "output.system.mtx");
    /* destroy resources, matrix, vector and solver */
    //AMGX_solver_destroy(solver);
    //AMGX_vector_destroy(x);
    //AMGX_vector_destroy(b);
    //AMGX_matrix_destroy(A);
    //AMGX_resources_destroy(rsrc);
    /* destroy config (need to use AMGX_SAFE_CALL after this point) */
    //AMGX_SAFE_CALL(AMGX_config_destroy(cfg));
    /* shutdown and exit */
    //AMGX_SAFE_CALL(AMGX_finalize_plugins());
    //AMGX_SAFE_CALL(AMGX_finalize());
    /* close the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
    //amgx_libclose(lib_handle);
#endif
    //CUDA_SAFE_CALL(cudaDeviceReset());
    //return status;
}

int solveamg(double *rhs, double *sol) 
{
    AMGX_vector_upload(b, N, 1, rhs);
    
    // AMGX_vector_get_size(x, &sol_size, &sol_bsize);
    /* Input your initial guess, x */
    //AMGX_vector_set_zero(x, n, bsize_x);
    AMGX_vector_upload(x, N, 1, sol);


    /* solver setup */
    //AMGX_solver_setup(solver, A);
    /* solver solve */
    AMGX_solver_solve(solver, b, x);
    AMGX_solver_get_status(solver, &status);

    /* example of how to write the linear system to the output */
    AMGX_vector_download(x,sol);
}


int destroy_amgx() 
{
    /* destroy resources, matrix, vector and solver */
    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(x);
    AMGX_vector_destroy(b);
    AMGX_matrix_destroy(A);
    AMGX_resources_destroy(rsrc);
    /* destroy config (need to use AMGX_SAFE_CALL after this point) */
    AMGX_SAFE_CALL(AMGX_config_destroy(cfg));
    /* shutdown and exit */
    AMGX_SAFE_CALL(AMGX_finalize_plugins());
    AMGX_SAFE_CALL(AMGX_finalize());
    /* close the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
    amgx_libclose(lib_handle);
#endif
    //CUDA_SAFE_CALL(cudaDeviceReset());
    //return status;
}

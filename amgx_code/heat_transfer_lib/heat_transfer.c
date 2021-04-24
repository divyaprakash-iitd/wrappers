#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "amgaxb_headers.h"
#include "amgx_c.h"

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


//int main(int argc, const char **argv)
void solveAMG(double data[], int col_ind[], int row_ptr[], double rhs[], double xi[]) 
{
    //input matrix and rhs/solution
    int n = 0;
    int bsize_x = 0;
    int bsize_y = 0;
    int sol_size = 0;
    int sol_bsize = 0;
    //library handles
    AMGX_Mode mode;
    AMGX_config_handle cfg;
    AMGX_resources_handle rsrc;
    AMGX_matrix_handle A;
    AMGX_vector_handle b, x;
    AMGX_solver_handle solver;
    //status handling
    AMGX_SOLVE_STATUS status;


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
    AMGX_SAFE_CALL(AMGX_config_create(&cfg, "config_version=2, solver(amg_solver)=AMG"));

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

    int N = 4; int nnz = 10; int block_dimx = 1; int block_dimy = 1;
    AMGX_matrix_upload_all(A, N, nnz, block_dimx, block_dimy, row_ptr, col_ind, data, 0);

    
    AMGX_matrix_get_size(A, &n, &bsize_x, &bsize_y);

    /* Input your RHS vector */
    //double rhs[] = {-300, 0, 0, -100};
    //AMGX_pin_memory(rhs);
    AMGX_vector_upload(b, 4, 1, rhs);
    
    AMGX_vector_get_size(x, &sol_size, &sol_bsize);
    /* Input your initial guess, x */
    AMGX_vector_set_zero(x, n, bsize_x);


    /* solver setup */
    AMGX_solver_setup(solver, A);
    /* solver solve */
    AMGX_solver_solve(solver, b, x);
    /* example of how to change parameters between non-linear iterations */
    //AMGX_config_add_parameters(&cfg, "config_version=2, default:tolerance=1e-12");
    
    //AMGX_solver_solve(solver, b, x);
    AMGX_solver_get_status(solver, &status);
    /* example of how to write the linear system to the output */
    

    AMGX_write_system(A, b, x, "output.system.mtx");
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

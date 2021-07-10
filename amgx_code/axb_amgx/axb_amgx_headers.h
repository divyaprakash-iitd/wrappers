int initialize_amgx(int *crs_data, double *data, int *col_ind, int *row_ptr, double *rhs, double *sol); 
int solveamg(int *crs_data, double *rhs, double *sol);
int destroy_amgx(); 

void errAndExit(const char *err);
void print_callback(const char *msg, int length);
int solveamg(int *crs_data, double *data, int *col_ind, int *row_ptr, double *rhs, double *sol);

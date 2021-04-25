#include <stdio.h>
#include "amgaxb_headers.h"

int main () 
{
    int i;
    double data[] = {-2, 1, 1, -2, 1, 1, -2, 1, 1, -2 };
    int col_ind[] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    int row_ptr[] = {0, 2, 5, 8, 10}; 

    double rhs[] = {-300, 0, 0, -100};
    double sol[4];    

    solveAMG(data, col_ind, row_ptr, rhs, sol);
    for (i=0;i<4;i++)
    {
        printf("%lf\n",sol[i]);
    }
    return 0; 
}

#include<stdio.h>
#include "my_headers.h"

int factorial(int n)
{   
    int f, i;
    f = 1;
    for (i = 1; i <=n; i++)
    {
        f = f*i;   
    }
    return f;
}


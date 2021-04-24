program solveaxb
    implicit none
    ! Declare variables
    integer :: nnz, col, row, n_x
    real(8), dimension(:), allocatable :: mtx, rhs, xinit
    integer, dimension(:), allocatable :: col_ind, row_ptr

    nnz = 10; col = 10; row = 5; n_x = 4;
    allocate(mtx(nnz),col_ind(col),row_ptr(row),rhs(n_x),xinit(n_x))

    mtx = (/-2, 1, 1, -2, 1, 1, -2, 1, 1, -2/)
    col_ind = (/0, 1, 0, 1, 2, 1, 2, 3, 2, 3/)
    row_ptr = (/0, 2, 5, 8, 10/)

    rhs = (/-300, 0, 0, -100/)
    xinit = (/0, 0, 0, 0/)
   
    ! To-do (using AMGX library) 
 
    !AMGX_matrix_upload_all(A, N, nnz, block_dimx, block_dimy, row_ptr, col_ind, data, 0);
    !AMGX_vector_upload(b, 4, 1, rhs);
    ! AMGX_solver_solve(solver, b, x);
    ! AMGX_vector_download(vector, data);

    print *, mtx
    print *, col_ind
    print *, row_ptr
    print *, rhs
    print *, xinit
        
end program solveaxb

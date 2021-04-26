program solveaxb
    use iso_c_binding, only: c_int, c_double, c_loc
    use ftn_c
    !implicit none
    ! Declare variables
    integer :: nnz, col, row, b_size
    real (c_double), dimension(:), allocatable, target :: mtx, rhs, sol
    integer (c_int), dimension(:), allocatable, target :: col_ind, row_ptr

    nnz = 10; col = 10; row = 5; b_size = 4;
    allocate(mtx(nnz),col_ind(col),row_ptr(row),rhs(b_size),sol(b_size))

    mtx = (/-2, 1, 1, -2, 1, 1, -2, 1, 1, -2/)
    col_ind = (/0, 1, 0, 1, 2, 1, 2, 3, 2, 3/)
    row_ptr = (/0, 2, 5, 8, 10/)

    rhs = (/-300, 0, 0, -100/)
    sol = (/-300, 0, 0, -100/)
    
    call solveamg(c_loc(mtx), c_loc(col_ind), c_loc(row_ptr),c_loc(rhs),c_loc(sol))  
 
   ! print *, mtx
   ! print *, col_ind
   ! print *, row_ptr
   ! print *, rhs
        
end program solveaxb

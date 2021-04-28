program solveaxb
    use iso_c_binding, only: c_int, c_double, c_loc
    use ftn_c
    implicit none
    ! Declare variables
    integer :: nnz, col, row, b_size, i
    real (c_double), dimension(:), allocatable, target :: datam, rhs, sol
    integer (c_int), dimension(:), allocatable, target :: col_ind, row_ptr, crs_data

    nnz = 10; col = 10; row = 5; b_size = 4;
    allocate(datam(nnz),col_ind(col),row_ptr(row),rhs(b_size),sol(b_size),crs_data(b_size))

    datam = (/-2, 1, 1, -2, 1, 1, -2, 1, 1, -2/)
    col_ind = (/0, 1, 0, 1, 2, 1, 2, 3, 2, 3/)
    row_ptr = (/0, 2, 5, 8, 10/)

    rhs = (/-300, 0, 0, -100/)
    
    sol = (/0, 0, 0, 0/)
    
    ! crs data = (N, nnz, block_dimx, block_dimy)
    crs_data = (/4, 10, 1, 1/)
    
    write(*,*) solveamg(c_loc(crs_data), c_loc(datam), c_loc(col_ind), c_loc(row_ptr),c_loc(rhs),c_loc(sol)) 

    write(*,100) ( real(sol(i)), i=1,4 ) 
    100 format (4(F14.7))
 
   ! print *, datam
   ! print *, col_ind
   ! print *, row_ptr
   ! print *, rhs
        
end program solveaxb

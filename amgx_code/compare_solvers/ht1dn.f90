program ht1dn
    use iso_c_binding, only: c_int, c_double, c_loc
    use ftn_c
    implicit none
    ! Declare variables
    integer :: N ! Number of unknowns
    integer :: nnz, n_row_ptr, row, n_row, b_size, i, id
    real (c_double), dimension(:), allocatable, target :: datam, rhs, sol
    integer (c_int), dimension(:), allocatable, target :: col_ind, row_ptr, crs_data

    N = 1e6
    nnz = 2*2 + 3*(N-2);

    n_row = N+1; b_size = N;
    allocate(datam(nnz),col_ind(nnz),row_ptr(n_row),rhs(b_size),sol(b_size),crs_data(b_size))
    
    id = 1
    do row = 1,N    
        if (row.eq.1) then
            datam(id) = -2; col_ind(id) = 0; id = id + 1
            datam(id) =  1; col_ind(id) = 1; id = id + 1
        else if (row.eq.N) then
            datam(id) =  1; col_ind(id) = row - 3 + 1; id = id + 1
            datam(id) = -2; col_ind(id) = row - 3 + 2; id = id + 1
        else
            datam(id) =  1; col_ind(id) = row - 3 + 1; id = id + 1
            datam(id) = -2; col_ind(id) = row - 3 + 2; id = id + 1
            datam(id) =  1; col_ind(id) = row - 3 + 3; id = id + 1
        end if
    end do

    do n_row_ptr = 2,N
        row_ptr(n_row_ptr) = 2 + 3*(n_row_ptr-2) 
    end do 

    row_ptr(1) = 0
    row_ptr(N+1) = nnz 

    rhs = 0.0d0
    rhs(1) = -300
    rhs(N) = -100
    
    sol = 0.0d0
    
    ! crs data = (N, nnz, block_dimx, block_dimy)
    crs_data = (/N, nnz, 1, 1/)
    
    write(*,*) solveamg(c_loc(crs_data), c_loc(datam), c_loc(col_ind), c_loc(row_ptr),c_loc(rhs),c_loc(sol)) 

    !write(*,100) ( real(sol(i)), i=1,N ) 
    !100 format (10000(F14.7))
 
    !print *, datam
    !print *, col_ind
    !print *, row_ptr
    !print *, rhs
        
end program ht1dn

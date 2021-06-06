module dataht
    implicit none

    ! Heat Transfer Parameters
    real(8), parameter :: eGen    = 5000   ! Uniform Heat Generation           [W/m^3]
    real(8), parameter :: k       = 2      ! Thermal Conductivity              [w/m/C]
    real(8), parameter :: hT      = 70     ! Heat Transfer Coefficient, Top    [W/m^2/K]
    real(8), parameter :: hB      = 10     ! Heat Transfer Coefficient, Bottom [W/m^2/K]
    real(8), parameter :: Tinf    = 25     ! Ambient Temperature               [C]
    real(8), parameter :: delta   = 0.01    ! Mesh Size                         [m]
    real(8), parameter :: qf      = 500    ! Heat Flux                         [W/m^2]
    real(8), parameter :: Tright  = 45     ! Right Boundary Temperature        [C]
    real(8), parameter :: Lx      = 1.5    ! X-Dimension of Plate              [m]
    real(8), parameter :: Ly      = 1      ! Y-Dimension of Plate              [m]
 
    ! Grid Parameters
    integer :: Ni    ! Number of Grid Points (x-direction)
    integer :: Nj    ! Number of Grid Points (y-direction)
    integer :: N     ! Number of Unknowns

    ! Declare Matrices and Vectors
    real(8), dimension(:), allocatable      :: A    ! Coefficient Matrix
    real(8), dimension(:), allocatable      :: b    ! RHS Vector
    real(8), dimension(:), allocatable      :: T    ! Unknown (Temperature Distribution)

    ! Declare Loop Indices variables 
    integer(8) :: i, j, row
    integer(8) :: LEFT, RIGHT, TOP, BOTTOM, CENTER

    ! Declare Function Variable
    !integer(8) :: id
    
    contains
    subroutine init() 
        ! Calculate Number of Grid Points and Total Number of Unknowns
        Ni = ceiling(Lx/delta) + 1
        Nj = ceiling(Ly/delta) + 1
        N  = (Ni-1)*Nj 

        ! Allocate Matrices and Vectors
        allocate(A(N))
        allocate(b(N))
        allocate(T(N))
        
        ! Initialize  Matrices and Vectors
        A = 0.0d0
        b = 0.0d0
        T = 0.0d0
    end subroutine
    
    subroutine gauss_siedel(n,A,b,T)
    implicit none
        
        ! Calling parameters
        integer, intent(in)                  :: n
        real(8), intent(in), dimension(n,n)     :: A
        real(8), intent(in), dimension(n)       :: b
        real(8), intent(inout), dimension(n)    :: T

        ! Local variables
        real(8)                     :: S, error, tol
        real(8)                     :: w
        integer(8)                  :: iter, maxIter, i, j
        real(8), dimension(n)       :: xold
        
        error   = 1.0
        tol     = 1e-10
        w       = 1.8     ! Relaxation Factor 
        maxIter = 10000
        iter    = 0

        do while ((error > tol).and.(iter < maxIter))
            xold = T
            do i = 1,n
                ! Calculate S
                S = 0.0d0
                do j = 1,n
                    if (i.ne.j) then
                        S = S + A(i,j)*T(j)
                    end if     
                end do
                T(i) = (1-w)*T(i) + w/a(i,i)*(b(i) - S)
            end do
            error = norm2(T-xold)
            iter = iter + 1
        end do
        if (iter == maxIter) then
            print *, "Did not converge!"
        end if
        print *, "error: ", error
        print *, "iter: ", iter
    end subroutine gauss_siedel
    
    function id(i,j)
    implicit none
        ! Dummy arguments
        integer(8) :: id

        ! Local variables
        integer(8)  :: i,j

        id = (j-1)*(Ni-1) + i     
        !print *,Ni
    end function id      
end module dataht

program ht2d
    use dataht
   
    !------------Fortran to C interface-------------!
    use iso_c_binding, only: c_int, c_double, c_loc
    use ftn_c
    !-----------------------------------------------!
    
    implicit none
   
    !-----------------------------------------------Variables for AMGX solver--------------------------------------------------!
    integer :: NNZ ! Number of non-zero (NNZ) elements in the coefficient matrix (C)
    integer :: NZN ! Variable used to iterate through the NNZ elements
    integer :: flag! Flag used to mark the first NNZ element of each row (used by row_ptr array)
    integer :: m   ! Index used to iterate through the columns of a row of C

    ! Compressed Row Storage (CRS) format matrix data
    real (c_double), dimension(:), allocatable, target :: val       ! Stores the values of the NNZ of C 
    integer (c_int), dimension(:), allocatable, target :: col_ind   ! Stores the column indices of the elements in val
    integer (c_int), dimension(:), allocatable, target :: row_ptr   ! Stores the locations in the val array that start a row
    
    real (c_double), dimension(:), allocatable, target :: rhs       ! Stores the RHS vector
    real (c_double), dimension(:), allocatable, target :: sol       ! Stores the solution vector
    integer (c_int), dimension(:), allocatable, target :: crs_data  ! Stores the CRS data needed by the AMGX solver
    !--------------------------------------------------------------------------------------------------------------------------!

    ! Declare Time Variables
    real :: start, finish

    ! Miscellaneous Variables
    integer :: ret ! Stores the return value of AMGX function
 
    ! Initialize variables
    call init()
    
    ! Count the Number of Non-Zero Elements (NNZ)
    NNZ = 0
    NNZ = (Ni-3)*(Nj-2)*5                   ! Internal Nodes
    NNZ = NNZ + (2*(Ni-3) + 2*(Nj-2))*4     ! Boundary Nodes
    NNZ = NNZ + 4*3                         ! Corner Nodes
    
    ! Allocate arrays for CRS matrix
    allocate(val(NNZ),col_ind(NNZ),row_ptr(N+1))
    row_ptr(N+1) = NNZ + 1 
    
    ! Fill up the Coefficient Matrix
    row = 1
    NNZ = 1
    do j = 1,Nj
        do i = 1,Ni-1
            A = 0.0d0  ! Initialize the coefficient matrix with zeros

            LEFT    = i-1
            RIGHT   = i+1
            BOTTOM  = j-1
            TOP     = j+1

            ! Internal Nodes
            if (LEFT /=0 .and. RIGHT/=Ni .and. TOP/=Nj+1 .and. BOTTOM/=0) then 
                A(id(LEFT,j))   =  1
                A(id(RIGHT,j))  =  1  
                A(id(i,BOTTOM)) =  1  
                A(id(i,TOP))    =  1 
                A(id(i,j))      = -4   ! Center
                b(row)          = -eGen*delta**2/k
            end if
            
            ! Boundary Nodes
            ! Left
            if (LEFT == 0 .and. TOP /=Nj+1 .and. BOTTOM /= 0) then
                A(id(RIGHT,j))  =  2  
                A(id(i,BOTTOM)) =  1  
                A(id(i,TOP))    =  1 
                A(id(i,j))      = -4   ! Center
                b(row)          = -eGen*delta**2/k - 2*qf*delta/k
            end if
            
            ! Right
            if (RIGHT == Ni .and. TOP /= Nj+1 .and. BOTTOM /= 0) then
                A(id(LEFT,j))   =  1
                A(id(i,BOTTOM)) =  1  
                A(id(i,TOP))    =  1 
                A(id(i,j))      = -4   ! Center
                b(row)          = -eGen*delta**2/k -Tright
            end if
            
            
            ! Top
            if (TOP == Nj + 1 .and. LEFT /=0 .and. RIGHT /=Ni) then
                A(id(LEFT,j))   =  1
                A(id(RIGHT,j))  =  1  
                A(id(i,BOTTOM)) =  2
                A(id(i,j))      = -4 - 2*delta*hT/k
                b(row)          = -eGen*delta**2/k - 2*delta*hT/k*Tinf
            end if
            
            ! Bottom
            if (BOTTOM == 0 .and. LEFT /=0 .and. RIGHT /=Ni) then
                A(id(LEFT,j))   =  1
                A(id(RIGHT,j))  =  1
                A(id(i,TOP))    =  2
                A(id(i,j))      = -4 - 2*delta*hB/k
                b(row)          = -eGen*delta**2/k - 2*delta*hB/k*Tinf
            end if
            
            ! Corner Nodes
            if (TOP == Nj+1 .and. LEFT == 0) then ! Top-Left
                A(id(RIGHT,j))  =  2  
                A(id(i,BOTTOM)) =  2   
                A(id(i,j))      = -4 - 2*delta*hT/k;   ! Center
                b(row)          = -eGen*delta**2/k - 2*qf*delta/k - 2*delta*hT/k*Tinf
            end if
            
            if (TOP == Nj+1 .and. RIGHT == Ni) then ! Top-Right
                A(id(LEFT,j))   =  1
                A(id(i,BOTTOM)) =  2  
                A(id(i,j))      = -4-2*delta*hT/k   ! Center
                b(row)          = -eGen*delta**2/k - 2*delta*hT/k*Tinf - Tright
            end if
            
            if (BOTTOM == 0 .and. LEFT == 0) then ! Bottom-Left
                A(id(RIGHT,j))  =  2   
                A(id(i,TOP))    =  2 
                A(id(i,j))      = -4 - 2*delta*hB/k   ! Center
                b(row)          = -eGen*delta**2/k - 2*qf*delta/k - 2*delta*hB/k*Tinf
            end if
                 
            if (BOTTOM == 0 .and. RIGHT == Ni) then ! Bottom-Right
                A(id(LEFT,j))   =  1
                A(id(i,TOP))    =  2  
                A(id(i,j))      = -4 - 2*delta*hB/k   ! Center
                b(row)          = -eGen*delta**2/k - 2*delta*hB/k*Tinf - Tright
            end if
                  
            ! Inspect the row for non-zero elements               
            flag = 0
            do m = 1,N
                if (abs(A(m)) > 1e-14) then
                    if (flag == 0) then
                        row_ptr(row) = NNZ - 1
                        flag = 1
                    end if
                
                    val(NNZ) = A(m)
                    col_ind(NNZ) = m - 1
                    NNZ = NNZ + 1
                end if
            end do

            row = row + 1 
        end do
    end do
    
    ! Solve using the AMGX Solver
    crs_data = (/N, NNZ-1, 1, 1/)
    allocate(sol(N))
    sol = 1
    val = val
    rhs = b

    call cpu_time(start) 
    ret = solveamg(c_loc(crs_data), c_loc(val), c_loc(col_ind), c_loc(row_ptr),c_loc(rhs),c_loc(sol))
    call cpu_time(finish)
    print *, 'Ni: ', Ni
    print *, 'Nj: ', Nj
    
    !open(unit=2, file='A.txt', ACTION="write", STATUS="replace")
    !do i=1,N
    !    write(2, '(*(F14.7))')( real(A(i,j)) ,j=1,N)
    !end do
    
    !open(unit=3, file='b.txt', ACTION="write", STATUS="replace")
    !do i=1,N
    !    write(3, '(*(F14.7))') real(b(i))
    !end do

    print *, 'N: ', N
   
    !call cpu_time(start) 
    !call gauss_siedel(n,A,b,T)
    !call cpu_time(finish)

    print '("Time: ",f6.3," seconds.")',finish-start
    
    open(unit=11, file='T_AMGX.txt', ACTION="write", STATUS="replace")
    do i=1,N
        write(11, '(*(F14.7))') real(sol(i))
    end do
end program ht2d


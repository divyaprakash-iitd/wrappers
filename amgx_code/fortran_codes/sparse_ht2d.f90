module dataht
    implicit none

    ! Heat Transfer Parameters
    real(8), parameter :: eGen    = 5000   ! Uniform Heat Generation           [W/m^3]
    real(8), parameter :: k       = 2      ! Thermal Conductivity              [w/m/C]
    real(8), parameter :: hT      = 70     ! Heat Transfer Coefficient, Top    [W/m^2/K]
    real(8), parameter :: hB      = 10     ! Heat Transfer Coefficient, Bottom [W/m^2/K]
    real(8), parameter :: Tinf    = 25     ! Ambient Temperature               [C]
    real(8), parameter :: delta   = 0.1    ! Mesh Size                         [m]
    real(8), parameter :: qf      = 500    ! Heat Flux                         [W/m^2]
    real(8), parameter :: Tright  = 45     ! Right Boundary Temperature        [C]
    real(8), parameter :: Lx      = 1.5    ! X-Dimension of Plate              [m]
    real(8), parameter :: Ly      = 1      ! Y-Dimension of Plate              [m]

    ! Grid Parameters
    integer(8) :: Ni    ! Number of Grid Points (x-direction)
    integer(8) :: Nj    ! Number of Grid Points (y-direction)

    ! Declare Matrices and Vectors
    real(8), dimension(:,:,:), allocatable  :: A    ! Coefficient Matrix
    real(8), dimension(:,:), allocatable    :: b    ! RHS Vector
    real(8), dimension(:,:), allocatable    :: T    ! Unknown (Temperature Distribution)
    real(8), dimension(:,:), allocatable    :: TOUT ! Output
    ! Declare Loop Indices variables 

    integer(8) :: i, j, row
    integer(8) :: LEFT, RIGHT, UP, DOWN, CENTER
    
    ! Coefficients
    integer, parameter :: L  = 1
    integer, parameter :: R  = 2
    integer, parameter :: U  = 3
    integer, parameter :: D  = 4
    integer, parameter :: C  = 5
    integer, parameter :: N  = 5
    
    contains
    subroutine init() 
        ! Calculate Number of Grid Points and Total Number of Unknowns
        Ni = ceiling(Lx/delta) + 1
        Nj = ceiling(Ly/delta) + 1

        ! Allocate Matrices and Vectors
        allocate(A(Ni,Nj,N))
        allocate(b(Ni,Nj))
        allocate(T(0:Ni+1,0:Nj+1))
        allocate(TOUT(Ni,Nj))
        
        ! Initialize  Matrices and Vectors
        A       = 0.0d0
        b       = 0.0d0
        T       = 0.0d0
        TOUT    = 0.0d0
        TOUT(Ni,:) = Tright 

    end subroutine
    
end module dataht

program ht2d
    use dataht
    implicit none
   
    ! Declare Time Variables
    real :: start, finish
    
    ! Define constants for Gauss-Seidel(SOR)
    real, parameter     :: OMEGA        = 1.8
    real, parameter     :: TOLERANCE    = 1E-10
    integer, parameter  :: MAXITER      = 10000
    integer             :: ITER
    real(8)             :: TTEMP
    
    ! Initialize variables
    call init()
    
    ! Fill up the Coefficient Matrix
    row = 1
    do j = 1,Nj
        do i = 1,Ni-1

            LEFT    = i-1
            RIGHT   = i+1
            DOWN    = j-1
            UP      = j+1

            ! Internal Nodes
            if (LEFT /=0 .and. RIGHT/=Ni .and. UP/=Nj+1 .and. DOWN/=0) then 
                A(i,j,L) =  1
                A(i,j,R) =  1  
                A(i,j,D) =  1  
                A(i,j,U) =  1 
                A(i,j,C) = -4   ! Center
                b(i,j)   = -eGen*delta**2/k
            end if
            
            ! Boundary Nodes
            ! Left
            if (LEFT == 0 .and. UP /=Nj+1 .and. DOWN /= 0) then
                A(i,j,R) =  2  
                A(i,j,D) =  1  
                A(i,j,U) =  1 
                A(i,j,C) = -4   ! Center
                b(i,j)   = -eGen*delta**2/k - 2*qf*delta/k
            end if
            
            ! Right
            if (RIGHT == Ni .and. UP /= Nj+1 .and. DOWN /= 0) then
                A(i,j,L) =  1
                A(i,j,D) =  1  
                A(i,j,U) =  1 
                A(i,j,C) = -4   ! Center
                b(i,j)   = -eGen*delta**2/k -Tright
            end if
            
            
            ! Top
            if (UP == Nj + 1 .and. LEFT /=0 .and. RIGHT /=Ni) then
                A(i,j,L) =  1
                A(i,j,R) =  1  
                A(i,j,D) =  2
                A(i,j,C) = -4 - 2*delta*hT/k
                b(i,j)   = -eGen*delta**2/k - 2*delta*hT/k*Tinf
            end if
            
            ! Bottom
            if (DOWN == 0 .and. LEFT /=0 .and. RIGHT /=Ni) then
                A(i,j,L)  =  1
                A(i,j,R)  =  1
                A(i,j,U)  =  2
                A(i,j,C)  = -4 - 2*delta*hB/k
                b(i,j)    = -eGen*delta**2/k - 2*delta*hB/k*Tinf
            end if
            
            ! Corner Nodes
            if (UP == Nj+1 .and. LEFT == 0) then ! Top-Left
                A(i,j,R) =  2  
                A(i,j,D) =  2   
                A(i,j,C) = -4 - 2*delta*hT/k;   ! Center
                b(i,j)   = -eGen*delta**2/k - 2*qf*delta/k - 2*delta*hT/k*Tinf
            end if
            
            if (UP == Nj+1 .and. RIGHT == Ni) then ! Top-Right
                A(i,j,L) =  1
                A(i,j,D) =  2  
                A(i,j,C) = -4-2*delta*hT/k   ! Center
                b(i,j)   = -eGen*delta**2/k - 2*delta*hT/k*Tinf - Tright
            end if
            
            if (DOWN == 0 .and. LEFT == 0) then ! Bottom-Left
                A(i,j,R) =  2   
                A(i,j,U) =  2 
                A(i,j,C) = -4 - 2*delta*hB/k   ! Center
                b(i,j)   = -eGen*delta**2/k - 2*qf*delta/k - 2*delta*hB/k*Tinf
            end if
                 
            if (DOWN == 0 .and. RIGHT == Ni) then ! Bottom-Right
                A(i,j,L) =  1
                A(i,j,U) =  2  
                A(i,j,C) = -4 - 2*delta*hB/k   ! Center
                b(i,j)   = -eGen*delta**2/k - 2*delta*hB/k*Tinf - Tright
            end if
                   
            row = row + 1 
        end do
    end do
    
    ! Gauss-Seidel Method (Succesive Over-Relaxation)
    ITER = 1
    do while ((ITER < MAXITER)) !.and.(resnorm > TOLERANCE))
        do j = 1,Nj
            do i = 1,Ni-1

                LEFT    = i-1
                RIGHT   = i+1
                DOWN    = j-1
                UP      = j+1

                ! Multiply the coefficients with the corresponding temperature
                ! value and the source terms
                   TTemp  =               b(i,j)    - &
                            (A(i,j,L) * T(LEFT,j)   + &
                             A(i,j,R) * T(RIGHT,j)  + &
                             A(i,j,U) * T(i,UP)     + &
                             A(i,j,D) * T(i,DOWN))
                   TTemp = TTemp/A(i,j,C)
                T(i,j) = (1-OMEGA)*T(i,j) + OMEGA*TTemp
             end do
          end do
         ITER = ITER + 1
     end do
  print *,'iter = ', iter 
    !open(unit=2, file='A.txt', ACTION="write", STATUS="replace")
    !do i=1,N
    !    write(2, '(*(F14.7))')( real(A(i,j)) ,j=1,N)
    !end do
    
    !open(unit=3, file='b.txt', ACTION="write", STATUS="replace")
    !do i=1,N
    !    write(3, '(*(F14.7))') real(b(i))
    !end do

    !call cpu_time(start) 
    !call gauss_siedel(n,A,b,T)
    !call cpu_time(finish)

    !print '("Time: ",f6.3," seconds.")',finish-start
   
    print *,'Ni = ', Ni 
    print *,'Nj = ', Nj 
    
    TOUT(1:Ni-1,:) = T(1:Ni-1,1:Nj)

    open(unit=11, file='T_GS_SPARSE.txt', ACTION="write", STATUS="replace")
    do j=1,Nj
        write(11, '(*(F14.7))') (real(TOUT(i,j)), i=1,Ni)
    end do

end program ht2d


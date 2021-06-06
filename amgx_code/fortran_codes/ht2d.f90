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
    integer(8) :: N     ! Number of Unknowns

    ! Declare Matrices and Vectors
    real(8), dimension(:,:), allocatable    :: A    ! Coefficient Matrix
    real(8), dimension(:), allocatable      :: b    ! RHS Vector
    real(8), dimension(:), allocatable      :: T    ! Unknown (Temperature Distribution)

    ! Declare Loop Indices variables 
    integer(8) :: i, j, row
    integer(8) :: LEFT, RIGHT, TOP, BOTTOM, CENTER

    ! Declare Function Variable
    !integer(8) :: id
   
    ! Variables to run the gs_sor subroutine
    real(8) :: eps   = 1e-6
    real(8) :: omega = 0.9
    integer, parameter :: iter  = 1000
    
    contains
    subroutine init() 
        ! Calculate Number of Grid Points and Total Number of Unknowns
        Ni = ceiling(Lx/delta) + 1
        Nj = ceiling(Ly/delta) + 1
        N  = (Ni-1)*Nj 

        ! Allocate Matrices and Vectors
        allocate(A(N,N))
        allocate(b(N))
        allocate(T(N))
        
        ! Initialize  Matrices and Vectors
        A = 0.0d0
        b = 0.0d0
        T = 0.0d0
    end subroutine
    
    !subroutine gauss_siedel()
    !implicit none
    !    
    !    ! Calling parameters
    !    !integer, intent(in) :: n
    !    !real(8), intent(in), dimension(n,n)  :: A
    !    !real(8), intent(in), dimension(n) :: b
    !    !real(8), intent(inout), dimension(n) :: T

    !    ! Local variables
    !    real(8) :: S, error, tol
    !    real(8) :: w
    !    integer :: iter, maxIter, i, j
    !    real(8), dimension(n) :: Told
    !    
    !    error = 1.0
    !    tol = 0.5
    !    w = 0.1     ! Relaxation Factor 
    !    maxIter = 1
    !    iter = 0
    !    do while ((error > tol).or.(iter < maxIter))
    !        Told = T
    !        do i = 1,n
    !            ! Calculate S
    !            S = 0.0d0
    !            do j = 1,n
    !                if (i.ne.j) then
    !                    S = S + A(i,j)*T(j)
    !                    print *, S    
    !                end if     
    !            end do
    !            T(i) = (1-w)*T(i) + w/a(i,i)*(b(i) - S)
    !        end do
    !        error = norm2(T-Told)
    !        iter = iter + 1
    !    end do
    !    if (iter == maxIter) then
    !        print *, "Did not converge!"
    !    end if
    !    print *, "error: "
    !    print *, error
    !end subroutine gauss_siedel

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
    implicit none
    
    ! Initialize variables
    call init()
    !print *,N
    ! Fill up the Coefficient Matrix
    row = 1
    do j = 1,Nj
        do i = 1,Ni-1
            
            LEFT    = i-1
            RIGHT   = i+1
            BOTTOM  = j-1
            TOP     = j+1

            ! Internal Nodes
            if (LEFT /=0 .and. RIGHT/=Ni .and. TOP/=Nj+1 .and. BOTTOM/=0) then 
                A(row,id(LEFT,j))   =  1
                A(row,id(RIGHT,j))  =  1  
                A(row,id(i,BOTTOM)) =  1  
                A(row,id(i,TOP))    =  1 
                A(row,id(i,j))      = -4   ! Center
                b(row)              = -eGen*delta**2/k
            end if
            
            ! Boundary Nodes
            ! Left
            if (LEFT == 0 .and. TOP /=Nj+1 .and. BOTTOM /= 0) then
                A(row,id(RIGHT,j))  =  2  
                A(row,id(i,BOTTOM)) =  1  
                A(row,id(i,TOP))    =  1 
                A(row,id(i,j))      = -4   ! Center
                b(row)              = -eGen*delta**2/k - 2*qf*delta/k
            end if
            
            ! Right
            if (RIGHT == Ni .and. TOP /= Nj+1 .and. BOTTOM /= 0) then
                A(row,id(LEFT,j))   =  1
                A(row,id(i,BOTTOM)) =  1  
                A(row,id(i,TOP))    =  1 
                A(row,id(i,j))      = -4   ! Center
                b(row)              = -eGen*delta**2/k -Tright
            end if
            
            
            ! Top
            if (TOP == Nj + 1 .and. LEFT /=0 .and. RIGHT /=Ni) then
                A(row,id(LEFT,j))   =  1
                A(row,id(RIGHT,j))  =  1  
                A(row,id(i,BOTTOM)) =  2
                A(row,id(i,j))      = -4 - 2*delta*hT/k
                b(row)              = -eGen*delta**2/k - 2*delta*hT/k*Tinf
            end if
            
            ! Bottom
            if (BOTTOM == 0 .and. LEFT /=0 .and. RIGHT /=Ni) then
                A(row,id(LEFT,j))   =  1
                A(row,id(RIGHT,j))  =  1
                A(row,id(i,TOP))    =  2
                A(row,id(i,j))      = -4 - 2*delta*hB/k
                b(row)              = -eGen*delta**2/k - 2*delta*hB/k*Tinf
            end if
            
            ! Corner Nodes
            if (TOP == Nj+1 .and. LEFT == 0) then ! Top-Left
                A(row,id(RIGHT,j))  =  2  
                A(row,id(i,BOTTOM)) =  2   
                A(row,id(i,j))      = -4 - 2*delta*hT/k;   ! Center
                b(row)              = -eGen*delta**2/k - 2*qf*delta/k - 2*delta*hT/k*Tinf
            end if
            
            if (TOP == Nj+1 .and. RIGHT == Ni) then ! Top-Right
                A(row,id(LEFT,j))   =  1
                A(row,id(i,BOTTOM)) =  2  
                A(row,id(i,j))      = -4-2*delta*hT/k   ! Center
                b(row)              = -eGen*delta**2/k - 2*delta*hT/k*Tinf - Tright
            end if
            
            if (BOTTOM == 0 .and. LEFT == 0) then ! Bottom-Left
                A(row,id(RIGHT,j))  =  2   
                A(row,id(i,TOP))    =  2 
                A(row,id(i,j))      = -4 - 2*delta*hB/k   ! Center
                b(row)              = -eGen*delta**2/k - 2*qf*delta/k - 2*delta*hB/k*Tinf
            end if
                 
            if (BOTTOM == 0 .and. RIGHT == Ni) then ! Bottom-Right
                A(row,id(LEFT,j))   =  1
                A(row,id(i,TOP))    =  2  
                A(row,id(i,j))      = -4 - 2*delta*hB/k   ! Center
                b(row)              = -eGen*delta**2/k - 2*delta*hB/k*Tinf - Tright
            end if
                   
            row = row + 1 
        end do
    end do
   
    print *, size(A,1)
 
    open(unit=2, file='A.txt', ACTION="write", STATUS="replace")
    do i=1,N
        write(2, '(27225F14.7)')( real(A(i,j)) ,j=1,N)
    end do
    
    open(unit=3, file='b.txt', ACTION="write", STATUS="replace")
    do i=1,N
        write(3, '(165F14.7)') real(b(i))
    end do

    print *, 'Ni = ', Ni
    print *, 'Nj = ', Nj
    
    call gauss_siedel(N,A,b,T)
    ! call gs_sor(A,b,T,omega,eps,n,iter)
    
    open(unit=11, file='T.txt', ACTION="write", STATUS="replace")
    do i=1,N
        write(11, '(165F14.7)') real(T(i))
    end do
end program ht2d

subroutine gauss_siedel(n,A,b,x)
implicit none
    
    ! Calling parameters
    integer(8), intent(in) :: n
    real(8), intent(in), dimension(n,n)  :: A
    real(8), intent(in), dimension(n) :: b
    real(8), intent(inout), dimension(n) :: x

    ! Local variables
    real(8) :: S, error, tol
    real(8) :: w
    integer :: iter, maxIter, i, j
    real(8), dimension(n) :: xold
    
    error = 1.0
    tol = 1e-10
    w = 1.2     ! Relaxation Factor 
    maxIter = 1000
    iter = 0
    do while ((error > tol).and.(iter < maxIter))
        xold = x
        do i = 1,n
            ! Calculate S
            S = 0.0d0
            do j = 1,n
                if (i.ne.j) then
                    S = S + A(i,j)*x(j)
                end if     
            end do
            x(i) = (1-w)*x(i) + w/a(i,i)*(b(i) - S)
        end do
        error = norm2(x-xold)
        iter = iter + 1
    end do
    if (iter == maxIter) then
        print *, "Did not converge!"
    end if
    print *, "error: "
    print *, error
end subroutine gauss_siedel
subroutine gs_sor(a,b,x,omega,eps,n,iter)
!==========================================================
! Solutions to a system of linear equations A*x=b
! Method: The successive-over-relaxation (SOR)
! Alex G. (November 2009)
!----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - array of the right hand coefficients b
! x(n)   - solutions (initial guess)
! n      - number of equations (size of matrix A)
! omega  - the over-ralaxation factor
! eps    - convergence tolerance 
! output ...
! x(n)   - solutions
! iter   - number of iterations to achieve the tolerance
! coments ...
! kmax   - max number of allowed iterations
!==========================================================
implicit none 
integer, parameter::kmax=1000000
integer n
double precision a(n,n), b(n), x(n)
double precision c, omega, eps, delta, conv, sum
integer i, j, k, iter, flag

! check if the system is diagonally dominant
flag = 0
do i=1,n
  sum = 0.0
  do j=1,n
    if(i == j) cycle
    sum = sum+abs(a(i,j))
  end do
  if(abs(a(i,i)) < sum) flag = flag+1
end do
if(flag >0) write(*,*) 'The system is NOT diagonally dominant'    

do k=1,kmax
  conv = 0.0
  do i=1,n
    delta = b(i)
    do j=1,n
      delta = delta - a(i,j)*x(j)
    end do
    x(i) = x(i)+omega*delta/a(i,i)
    if(abs(delta) > conv) conv=abs(delta)
  end do
  if(conv < eps) exit
end do
iter = k
!if(k == kmax) write (*,*)'The system failed to converge'
if(k > kmax) write (*,*)'The system failed to converge'
end subroutine gs_sor

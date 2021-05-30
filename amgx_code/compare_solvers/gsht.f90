program gsht
    implicit none

    ! Declare variables
    integer, parameter :: N         = 1E6 ! Number of unknowns
    integer, parameter :: NN        = 3
    integer, parameter :: LEFT      = 1
    integer, parameter :: CENTER    = 2
    integer, parameter :: RIGHT     = 3
    integer, parameter :: MAXITER   = 1000
    
    integer :: i, ITER
    real(8) ::TL, TR, TTemp            
    real(8), dimension(N,NN) :: A
    real(8), dimension(N) :: T, B, R
    real(8) :: TOLERANCE, ER
    
    real :: start, finish

    ! Boundary conditions
    TL = 300
    TR = 100
    
    ! Initialize vectors
    T = 0.0d0
    R = 0.0d0
    B = 0.0d0 
    
    ! Fill up coefficient matrix    
    A(2:N-1,LEFT)   =  1
    A(2:N-1,RIGHT)  =  1
    A(:,CENTER)     = -2

    A(1,LEFT)  = 0
    A(1,RIGHT) = 1
    A(N,LEFT)  = 1 
    A(N,RIGHT) = 0

    B(1) = -TL
    B(N) = -TR

    
    call cpu_time(start)    
    ! Stopping criteria
    TOLERANCE = 1E-14
    ER = 1.0d0
    ITER = 0
    do while ((ITER < MAXITER).and.(ER > TOLERANCE))
        do i = 1,N
            if (i == 1) then
                TTemp = B(i) - A(i,LEFT)*TL - A(i,RIGHT)*T(i+1)
            else if (i == N) then
                TTemp = B(i) - A(i,LEFT)*T(i-1) - A(i,RIGHT)*TR
            else
                TTemp = B(i) - A(i,LEFT)*T(i-1) - A(i,RIGHT)*T(i+1)
            end if
            T(i) = 1/A(i,CENTER)*TTemp
        end do
        
        ! Calculate Residual: R = AT - B
        do i = 1,N
            if (i == 1) then
                R(i) = A(i,LEFT)*TL + A(i,CENTER)*T(i) + A(i,RIGHT)*T(i+1) - B(i)
            elseif (i == N) then
                R(i) = A(i,LEFT)*T(i-1) + A(i,CENTER)*T(i) + A(i,RIGHT)*TR - B(i)
            else
                R(i) = A(i,LEFT)*T(i-1) + A(i,CENTER)*T(i) + A(i,RIGHT)*T(i+1) - B(i)
            end if
        end do
        ITER = ITER + 1
        ER = norm2(R)
    end do
    call cpu_time(finish)    
    
    print '("Time = ",f6.3," seconds.")',finish-start
    print *, T(N)
end program gsht

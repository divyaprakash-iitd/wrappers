program coefficient_matrix
    integer, parameter :: Ni = 3
    integer, parameter :: Nj = 3
    integer, parameter :: Nk = 3
    
    integer(4) :: id, i, j, k, row ! Ni, Nj, Nk, row
    real(8),  dimension(Ni*Nj*Nk,Ni*Nj*Nk) :: C ! Coefficient Matrix
    real(8), dimension(Ni*Nj*Nk) :: deltayp, deltayv
    real(8) :: dx, dy, invdeltax2, invdeltaz2
    real(8), dimension(Ni*Nj*Nk) :: b, x

    dx = 0.10; dz = 0.1d0
    row = 1
   
    deltayp = 0.1d0
    deltayv = 0.1d0
     
    invdeltax2 = 1.0/dx**2
    invdeltaz2 = 1.0/dz**2    
        

    ! Coefficient Matrix
    C = 0.0d0
    do k =1,Nk
        do j = 1,Nj
            do i = 1,Ni
                ! x-direction
                C(row, idx(i-1,j,k,Ni,Nj,Nk)) = invdeltax2
                C(row, idx(i+1,j,k,Ni,Nj,Nk)) = invdeltax2

                ! z-direction
                C(row, idx(i,j,k-1,Ni,Nj,Nk)) = invdeltaz2
                C(row, idx(i,j,k+1,Ni,Nj,Nk)) = invdeltaz2

                ! y-direction
                if (j.eq.1) then
                    C(row, idx(i,j+1,k,Ni,Nj,Nk)) = 1/deltayv(j+1)/deltayp(j)
                    ! Center Node
                    C(row, idx(i,j,k,Ni,Nj,Nk))   = -2*invdeltax2 - 2*invdeltaz2 -1/deltayv(j+1)/deltayp(j)
                
                else if (j.eq.Nj) then
                    C(row, idx(i,j-1,k,Ni,Nj,Nk)) = 1/deltayv(j)/deltayp(j)
                    ! Center Node
                    C(row, idx(i,j,k,Ni,Nj,Nk))   = -2*invdeltax2 - 2*invdeltaz2 -1/deltayv(j)/deltayp(j)

                else
                    C(row, idx(i,j-1,k,Ni,Nj,Nk)) = 1/deltayv(j)/deltayp(j)
                    C(row, idx(i,j+1,k,Ni,Nj,Nk)) = 1/deltayv(j+1)/deltayp(j)
                    ! Center Node
                    C(row, idx(i,j,k,Ni,Nj,Nk))   = -2*invdeltax2 - 2*invdeltaz2 -1/deltayp(j)*(1/deltayv(j+1) + 1/deltayv(j))
                end if

                row = row + 1;
            end do
        end do
    end do
    
    !i =1; j = 1; k = 1
    !id = idx(i,j,k,Ni,Nj,Nk)
    !print *, "Transformed id: "
    !print *, id
    !print *, "invdeltax2: "
    !print *, invdeltax2
    
    b = 1.0
   ! x = gauss_siedel(Ni*Nj*Nk,C,b)    

    open(unit=2, file='graph1.txt', ACTION="write", STATUS="replace")
    do i=1,Ni*Nj*Nk
        write(2, '(27F14.7)')( real(C(i,j)) ,j=1,Ni*Nj*Nk)
    end do

    !print *, sum(C)
end program coefficient_matrix

function idx(i,j,k,Ni,Nj,Nk)
implicit none
    
    ! dummy arguments
    integer(4) :: idx

    ! local variables
    integer(4)  :: i,j,k,Ni,Nj,Nk
    
    ! Periodicity in x
    if (i.LT.1) then
        i = Ni
    else if (i.GT.Ni) then
        i = 1
    end if

    ! Periodicity in z
    if (k.LT.1) then
        k = Nk
    else if (k.GT.Nk) then
        k = 1
    end if

    idx = (k-1)*(Ni*Nj) + (j-1)*Ni + i     
end function idx      

function gauss_siedel(n,A,b)
!implicit none
    integer :: n, maxIter, i, j, iter
    real(8), dimension(n,n) :: A
    real(8), dimension(n) :: b, x, gauss_siedel
    real(8) :: S
    real(8) :: w

    x = 1.0 ! Initial Guess
    w = 1.2     ! Relaxation Factor 
    maxIter = 1000
    do iter = 1,maxIter
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
    end do
    gauss_siedel = x
end function gauss_siedel

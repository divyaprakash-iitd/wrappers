module dataht
    implicit none

    ! Heat Transfer Parameters
    real(8), parameter :: eGen    = 5000   ! Uniform Heat Generation           [W/m^3]
    real(8), parameter :: k       = 2      ! Thermal Conductivity              [w/m/C]
    real(8), parameter :: hT      = 70     ! Heat Transfer Coefficient, Top    [W/m^2/K]
    real(8), parameter :: hB      = 10     ! Heat Transfer Coefficient, Bottom [W/m^2/K]
    real(8), parameter :: Tinf    = 25     ! Ambient Temperature               [C]
    real(8), parameter :: delta   = 0.1   ! Mesh Size                         [m]
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

    contains
    subroutine init() 
        ! Calculate Number of Grid Points and Total Number of Unknowns
        Ni = ceiling(Lx/delta) + 1
        Nj = ceiling(Ly/delta) + 1
        N  = (Ni-1)*Nj 

        ! Allocate Matrices and Vectors
        allocate(A(Ni,Nj))
        allocate(b(N))
        allocate(T(N))
        
        ! Initialize  Matrices and Vectors
        A = 0.0d0
        b = 0.0d0
        T = 0.0d0
    end subroutine
    
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
    print *,Ni
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
end program ht2d



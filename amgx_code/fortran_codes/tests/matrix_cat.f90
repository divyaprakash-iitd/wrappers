program matrix
implicit none
real,dimension (3,3) :: mat1,mat2
real,dimension(18) :: mat3
integer i

mat1=reshape( (/1,2,3,4,5,6,7,8,9/),(/3,3/))
mat2=reshape( (/1,2,3,4,5,6,7,8,9/),(/3,3/))
mat3=[mat1,mat2]

print*, shape([mat1,mat2])  !check shape of concatenated array
!display
do i=1,18,1
write(*,10) mat3(i)
10 format(F10.4)
end do

end program

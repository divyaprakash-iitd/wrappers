program factorial_f
    use iso_c_binding, only: c_int
    use ftn_c
    
    integer (c_int) :: N
    N = 10
    write(*,*) factorial(N)
    
end program factorial_f

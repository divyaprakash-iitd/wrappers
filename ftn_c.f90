module ftn_c
    interface
        integer (c_int) function factorial(N) bind(C,'factorial')
            use iso_c_binding
            implicit none
            integer(c_int), value :: N
            end function factorial
    end interface    
end module ftn_c

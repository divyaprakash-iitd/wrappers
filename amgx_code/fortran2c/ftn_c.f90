module ftn_c
    interface
        integer (c_int) function solveamg(datam, col_ind, row_ptr, rhs, sol) bind(c, name='solveamg')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: datam
            type (c_ptr), value :: col_ind
            type (c_ptr), value :: row_ptr
            type (c_ptr), value :: rhs
            type (c_ptr), value :: sol
        end function solveamg
    end interface    
end module ftn_c

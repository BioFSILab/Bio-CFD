module ftn_c
    interface
        integer (c_int) function initialize_amgx() bind(c, name='initialize_amgx')
            use iso_c_binding
            implicit none

        end function initialize_amgx
        
        integer (c_int) function solveamg(crs_data, datam, col_ind, row_ptr, rhs,sol,nit) bind(c, name='solveamg')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: crs_data
            type (c_ptr), value :: datam
            type (c_ptr), value :: col_ind
            type (c_ptr), value :: row_ptr
            type (c_ptr), value :: rhs
            type (c_ptr), value :: sol   
            type (c_ptr), value :: nit

        end function solveamg
        

        integer (c_int) function solveamg1(rhs,sol,nit1) bind(c, name='solveamg1')
            use iso_c_binding
            implicit none
            
            type (c_ptr), value :: rhs
            type (c_ptr), value :: sol   
            type (c_ptr), value :: nit1

        end function solveamg1

        integer (c_int) function destroy_amgx() bind(c, name='destroy_amgx')
            use iso_c_binding
            implicit none
        end function destroy_amgx
    end interface    
end module ftn_c

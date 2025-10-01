MODULE biocfd_interfaces
       use, intrinsic :: iso_fortran_env, only: dp => real64, int64
       IMPLICIT NONE
       private

       public :: Interfaces

        type Interfaces

        INTEGER(int64), ALLOCATABLE, DIMENSION (:,:) :: px_interface_det, ux_interface_det, &
                                                        vx_interface_det, wx_interface_det, &
                                                        pz_interface_det, uz_interface_det, &
                                                        vz_interface_det, wz_interface_det, &
                                                        py_interface_det, uy_interface_det, &
                                                        vy_interface_det, wy_interface_det
        REAL(dp) :: a_blk, b_blk, xintf_start, xintf_end, yintf_start, yintf_end, &
                    zintf_start, zintf_end
        REAL(dp) :: xintf_st_new, xintf_en_new, yintf_st_new, yintf_en_new
        REAL(dp) :: zintf_st_new, zintf_en_new
        INTEGER(int64) :: a_msh, b_msh, a_intf, b_intf
        INTEGER(int64) :: counterxu,counteryu,counterxp,counteryp,counterxv,counteryv
        INTEGER(int64) :: counterxw,counteryw
        INTEGER(int64) :: counterzu,counterzp,counterzv, counterzw
        end type Interfaces

END MODULE biocfd_interfaces


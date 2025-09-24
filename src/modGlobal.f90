! allow(missing-accessibility-statement) - TODO: Attempt to remove in future
MODULE global
       use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int32, int64
       use iso_c_binding,only :c_int,c_double,c_loc,c_ptr
       use biocfd_blocks,only : Blocks
       IMPLICIT NONE
       CHARACTER (LEN = 128) :: line
       CHARACTER (LEN = 3)   :: char_f
       INTEGER               :: istart
       INTEGER (int64)   :: itamax, pcItaMax,amgxita,      &
                                ita, ita1,ita2,    &
                                inor,blk_start, coarse_flcnt_check
       REAL (dp)        ::      dt_order,  &
                                omega,omega1,omega2,omega3,omega4,  &
                                dxmin, &
                                freq, &
                                u0, &
                                epsi, re, rev, &
                                alpha, &
                                xfact,deltat,totime, totalTime, dfinish, dstart, &
                                pi, al, uc

       REAL (dp)       :: alpha_m, theta_m, alpha_m1, theta_m1, mu_f, rho_f, l_c, u_tip, disp

        INTEGER (int64) ::nblocks, intflines

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

        type(Blocks),allocatable ::block(:)
        type(Interfaces),allocatable ::intfr(:)

       REAL (dp)    :: phase_angle,a0y,aoa,piv_pt, aoa1,aoa2,ang_theta, alpha_t, theta_t
       INTEGER (int64)   :: surGeoPoints
END MODULE global


! allow(missing-accessibility-statement) - TODO: Attempt to remove in future
MODULE global
       use, intrinsic :: iso_fortran_env, only: dp => real64, int64
       use biocfd_blocks,only : Blocks
       use biocfd_interface_type, only: Interface_t
       IMPLICIT NONE
       INTEGER (int64)   ::     ita, ita1,   &
                                inor,blk_start, coarse_flcnt_check
       REAL (dp)        ::      dt_order,  &
                                omega,omega1,omega2,omega3,omega4,  &
                                dxmin, &
                                freq, &
                                u0, &
                                epsi, re, &
                                alpha, &
                                deltat, totime, totalTime, &
                                uc
       real(dp), parameter :: pi = 4._dp * atan(1._dp)

       REAL (dp)       ::  theta_m, alpha_m1, theta_m1, mu_f, rho_f, l_c, u_tip, disp

        INTEGER (int64) ::nblocks, intflines

        type(Blocks),allocatable ::block(:)
        type(Interface_t), allocatable :: intfr(:)

       REAL (dp)    :: phase_angle,aoa,piv_pt, aoa1,aoa2
END MODULE global


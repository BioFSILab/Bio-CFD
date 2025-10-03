! allow(missing-accessibility-statement) - TODO: Attempt to remove in future
MODULE global
       use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int32, int64
       use biocfd_blocks,only : Blocks
       use biocfd_interfaces, only : Interfaces
       IMPLICIT NONE
       CHARACTER (LEN = 128) :: line
       CHARACTER (LEN = 3)   :: char_f
       INTEGER               :: istart
       INTEGER (int64)   :: itamax, pcItaMax,      &
                                ita, ita1,ita2,    &
                                inor,blk_start, coarse_flcnt_check
       REAL (dp)        ::      dt_order,  &
                                omega,omega1,omega2,omega3,omega4,  &
                                dxmin, &
                                freq, &
                                u0, &
                                epsi, re, rev, &
                                alpha, &
                                xfact,deltat,totime, totalTime, &
                                pi, uc

       REAL (dp)       :: alpha_m, theta_m, alpha_m1, theta_m1, mu_f, rho_f, l_c, u_tip, disp

        INTEGER (int64) ::nblocks, intflines

        type(Blocks),allocatable ::block(:)
        type(Interfaces),allocatable ::intfr(:)

       REAL (dp)    :: phase_angle,aoa,piv_pt, aoa1,aoa2,ang_theta
       INTEGER (int64)   :: surGeoPoints
END MODULE global


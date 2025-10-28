! allow(missing-accessibility-statement) - TODO: Attempt to remove in future
MODULE global
       use, intrinsic :: iso_fortran_env, only: dp => real64, int64
       use biocfd_block_type, only: Block_t
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

        INTEGER (int64) ::nblocks, intflines

        type(Block_t),allocatable ::block(:)
        type(Interface_t), allocatable :: intfr(:)

END MODULE global


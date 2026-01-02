! allow(missing-accessibility-statement) - TODO: Attempt to remove in future
MODULE global
       use, intrinsic :: iso_fortran_env, only: dp => real64
       use biocfd_block_type, only: Block_t
       use biocfd_interface_type, only: Interface_t
       IMPLICIT NONE
       real(dp), parameter :: pi = 4._dp * atan(1._dp)

        type(Block_t),allocatable :: block(:)
        type(Interface_t), allocatable :: intfr(:)

END MODULE global

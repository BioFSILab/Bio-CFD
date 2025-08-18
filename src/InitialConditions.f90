module biocfd_initial_conditions
  use iso_fortran_env, only : dp => real64
  use global, only : block, ita, ita1, nblocks, totime, uc
  use biocfd_blocks, only : Blocks
  implicit none
  private

  public :: initialConditions

  contains
      SUBROUTINE initialConditions(blk)

        type(Blocks), intent(inout) :: blk
       INTEGER::  i, j, k, n
        WRITE(*,*) 'Enter initialcondtitions'
       !cell variables
        blk%u = uc
        blk%ut = uc
        blk%v = 0._dp
        blk%vt = 0._dp
        blk%w = 0._dp
        blk%wt = 0._dp
        blk%p = 0._dp
        blk%cell_pr=0
        ! SUMWSS = 0._dp            !SUMWSS global real array(nsurf)
        ! SIGNWSS = 0._dp           !SIGNWSS global real array(nsurf)
        ita = 0
        ita1 = 0
        totime = 0.
           DO n = 1, blk%fluidCellCount
              i = blk%fluidIndexPtr(n, 1)
              j = blk%fluidIndexPtr(n, 2)
              k = blk%fluidIndexPtr(n, 3)
              blk% u(i,j,k)  = uc  !396.33054782262406 !116.236233
              blk%v(i,j,k)  = 0.
              blk%w(i,j,k)  = 0._dp
              blk%ut(i,j,k) = uc
              blk%vt(i,j,k) = 0._dp
              blk%wt(i,j,k) = 0._dp
          END DO

        print*, 'initial'

      END SUBROUTINE initialConditions

end module biocfd_initial_conditions

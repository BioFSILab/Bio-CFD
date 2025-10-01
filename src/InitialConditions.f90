module biocfd_initial_conditions
  use, intrinsic :: iso_fortran_env, only : dp => real64
  use global, only : block, ita, ita1, nblocks, totime, uc
  implicit none
  private

  public :: initialConditions

  contains
      SUBROUTINE initialConditions

       INTEGER::  i, j, k, n, f
        WRITE(*,*) 'Enter initialcondtitions'
       !cell variables
        DO i=1,nblocks
        block(i)%u = uc
         block(i)%ut = uc
         block(i)%v = 0._dp
         block(i)%vt = 0._dp
        block(i)%w = 0._dp
        block(i)%wt = 0._dp
        block(i)%p = 0._dp
        block(i)%cell_pr=0
        end do
        ! SUMWSS = 0._dp            !SUMWSS global real array(nsurf)
        ! SIGNWSS = 0._dp           !SIGNWSS global real array(nsurf)
        ita = 0
        ita1 = 0
        totime = 0.
        DO f=1,nblocks
           DO n = 1, block(f)%fluidCellCount
              i = block(f)%fluidIndexPtr(n, 1)
              j = block(f)%fluidIndexPtr(n, 2)
              k = block(f)%fluidIndexPtr(n, 3)
              block(f)% u(i,j,k)  = uc  !396.33054782262406 !116.236233
              block(f)%v(i,j,k)  = 0.
              block(f)%w(i,j,k)  = 0._dp
              block(f)%ut(i,j,k) = uc
              block(f)%vt(i,j,k) = 0._dp
              block(f)%wt(i,j,k) = 0._dp
          END DO
        END DO

        print*, 'initial'

      END SUBROUTINE initialConditions

end module biocfd_initial_conditions

module biocfd_last_conditions
  use iso_fortran_env, only : dp => real64
  use global, only : block, ita, ita1, ita2, nblocks, re, totime
  use biocfd_blocks, only : Blocks
  implicit none
  private

  public :: lastConditions

  contains
      SUBROUTINE lastConditions(blk)

        type(Blocks), intent(inout) :: blk
       INTEGER::  i, j, k, g
        CHARACTER(len=150) :: filename3
        WRITE(*,*) 'Enter lastcondtitions'

        ! SUMWSS = 0._dp            !SUMWSS global real array(nsurf)
        ! SIGNWSS = 0._dp           !SIGNWSS global real array(nsurf)
        ita = 0
        ita1 = 0
        ita2 = 0

       WRITE(filename3,3) g, re
  3     FORMAT('out/aorta_chkpt.',i3.3,'.',f6.1,".dat")
       OPEN (1, FILE=filename3, FORM='formatted')
       DO k = 1, blk%nz+2
       DO j = 1, blk%ny+2
       DO i = 1, blk%nx+2
         READ(1,*) blk%u(i,j,k), blk%v(i,j,k), blk%w(i,j,k), blk%p(i,j,k), &
                   totime, ita, ita1
      END DO
      END DO
      END DO
       CLOSE(1)

       DO k = 1, blk%nz+2
       DO j = 1, blk%ny+2
       DO i = 1, blk%nx+2
           blk%ut(i,j,k) =  blk%u(i,j,k)
           blk%vt(i,j,k) =  blk%v(i,j,k)
           blk%wt(i,j,k) =  blk%w(i,j,k)

        END DO
        END DO
        END DO

      END SUBROUTINE lastConditions
end module biocfd_last_conditions

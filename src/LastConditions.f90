module biocfd_last_conditions
  use iso_fortran_env, only : dp => real64
  use global, only : block, ita, ita1, ita2, nblocks, re, totime
  implicit none
  private

  public :: lastConditions

  contains
      SUBROUTINE lastConditions

       INTEGER::  i, j, k, g
        CHARACTER(len=150) :: filename3
        WRITE(*,*) 'Enter lastcondtitions'

        ! SUMWSS = 0._dp            !SUMWSS global real array(nsurf)
        ! SIGNWSS = 0._dp           !SIGNWSS global real array(nsurf)
        ita = 0
        ita1 = 0
        ita2 = 0


        DO g=1,nblocks
       WRITE(filename3,3) g, re
  3     FORMAT('out/aorta_chkpt.',i3.3,'.',f6.1,".dat")
       OPEN (1, FILE=filename3, FORM='formatted')
       DO k = 1, block(g)%nz+2
       DO j = 1, block(g)%ny+2
       DO i = 1, block(g)%nx+2
         READ(1,*) block(g)%u(i,j,k), block(g)%v(i,j,k), block(g)%w(i,j,k), block(g)%p(i,j,k), &
                   totime, ita, ita1
      END DO
      END DO
      END DO
       CLOSE(1)
        END DO


       DO g=1,nblocks
       DO k = 1, block(g)%nz+2
       DO j = 1, block(g)%ny+2
       DO i = 1, block(g)%nx+2
           block(g)%ut(i,j,k) =  block(g)%u(i,j,k)
           block(g)%vt(i,j,k) =  block(g)%v(i,j,k)
           block(g)%wt(i,j,k) =  block(g)%w(i,j,k)

        END DO
        END DO
        END DO
        END DO

      END SUBROUTINE lastConditions
end module biocfd_last_conditions

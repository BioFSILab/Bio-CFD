module biocfd_last_conditions
  use, intrinsic :: iso_fortran_env, only : dp => real64, int64
  use global, only: ita, ita1, totime
  use biocfd_blocks, only: Blocks
  implicit none
  private

  public :: lastConditions

contains
  SUBROUTINE lastConditions(blk, id, re)

    !> The Block that we want to setup
    type(Blocks), intent(inout) :: blk
    !> The ID number of the block, usually 1 for the coarse block and
    !> >=2 for the fine blocks
    integer(int64), intent(in) :: id
    real(dp), intent(in) :: re
    INTEGER ::  i, j, k
    CHARACTER(len=150) :: filename3

    WRITE(filename3,3) id, re
3   FORMAT('out/aorta_chkpt.',i3.3,'.',f6.1,".dat")
    OPEN (1, FILE=filename3, FORM='formatted')
    DO k = 1, blk%nz+2
       DO j = 1, blk%ny+2
          DO i = 1, blk%nx+2
             ! TODO: totime, ita, and ita1 are all global variables,
             ! but they will always just be overwritten by the last
             ! block which is read. Do we want to have some validation
             ! here that things are working as one would expect?
             READ(1,*) blk%u(i,j,k), blk%v(i,j,k), blk%w(i,j,k), blk%p(i,j,k), &
                  totime, ita, ita1
          END DO
       END DO
    END DO
    CLOSE(1)

    DO k = 1, blk%nz+2
       DO j = 1, blk%ny+2
          DO i = 1, blk%nx+2
             blk%ut(i,j,k) = blk%u(i,j,k)
             blk%vt(i,j,k) = blk%v(i,j,k)
             blk%wt(i,j,k) = blk%w(i,j,k)
          END DO
       END DO
    END DO

  END SUBROUTINE lastConditions
end module biocfd_last_conditions

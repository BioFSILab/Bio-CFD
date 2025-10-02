module biocfd_initial_conditions
  use, intrinsic :: iso_fortran_env, only : dp => real64
  use biocfd_blocks, only: Blocks
  implicit none
  private

  public :: initialConditions

contains
  SUBROUTINE initialConditions(blk, uc)

    type(Blocks), intent(inout) :: blk
    real(dp), intent(in) :: uc
    integer :: i, j, k, n
    !cell variables
    blk%u = uc
    blk%ut = uc
    blk%v = 0._dp
    blk%vt = 0._dp
    blk%w = 0._dp
    blk%wt = 0._dp
    blk%p = 0._dp
    blk%cell_pr=0
    DO n = 1, blk%fluidCellCount
       i = blk%fluidIndexPtr(n, 1)
       j = blk%fluidIndexPtr(n, 2)
       k = blk%fluidIndexPtr(n, 3)
       blk%u(i,j,k) = uc  !396.33054782262406 !116.236233
       blk%v(i,j,k) = 0._dp
       blk%w(i,j,k) = 0._dp
       blk%ut(i,j,k) = uc
       blk%vt(i,j,k) = 0._dp
       blk%wt(i,j,k) = 0._dp
    END DO

  END SUBROUTINE initialConditions

end module biocfd_initial_conditions

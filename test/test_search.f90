module test_search
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use testdrive, only : error_type, unittest_type, new_unittest, check
  implicit none
  private

  public :: collect_search

contains

  !> Collect all exported unit tests
  subroutine collect_search(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("find_dist_node", test_find_dist_node)]
  end subroutine collect_search


  subroutine test_find_dist_node(error)
    use biocfd_blocks, only : Blocks
    use biocfd_search, only : findDistnode
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    integer(int64) :: expected
    integer :: ibnodes
    type(Blocks) :: blk

    ! Setup the required block variables
    ibnodes = 6
    blk%ibnodes = ibnodes

    ! Although not required by fortran>=2003 Lfortran requires arrays
    ! to be allocated before they are assigned to see
    ! https://github.com/lfortran/lfortran/issues/2946
    allocate(&
         blk%ibNodeId(ibnodes), &
         blk%xnode(ibnodes), &
         blk%ynode(ibnodes), &
         blk%znode(ibnodes) &
    )

    blk%ibNodeId = [51, 51, 51, 51, 51, 51]
    blk%xnode = [0, -7, -4, 1, -10, 8]
    blk%ynode = [0, -7, -4, 1, -10, 8]
    blk%znode = [0, -7, -4, 1, -10, 8]

    call findDistnode(blk)
    expected = 6
    call check(error, blk%mk, expected)
    if (allocated(error)) return
  end subroutine test_find_dist_node
end module test_search

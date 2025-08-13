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
    use biocfd_search, only : findDistnode
    use global, only : block, blk_start, nblocks
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    integer(int64) :: expected
    integer :: ibnodes

    ! Setup the required block variables
    blk_start = 1
    nblocks = 1
    ibnodes = 6
    allocate(block(nblocks))
    block(1)%ibnodes = ibnodes

    ! Although not required by fortran>=2003 Lfortran requires arrays
    ! to be allocated before they are assigned to see
    ! https://github.com/lfortran/lfortran/issues/2946
    allocate(&
         block(1)%ibNodeId(ibnodes), &
         block(1)%xnode(ibnodes), &
         block(1)%ynode(ibnodes), &
         block(1)%znode(ibnodes) &
    )
    
      
    block(1)%ibNodeId = [51, 51, 51, 51, 51, 51]
    block(1)%xnode = [0, -7, -4, 1, -10, 8]
    block(1)%ynode = [0, -7, -4, 1, -10, 8]
    block(1)%znode = [0, -7, -4, 1, -10, 8]

    call findDistnode
    expected = 6
    call check(error, block(1)%mk, expected)
    if (allocated(error)) return
  end subroutine test_find_dist_node
end module test_search

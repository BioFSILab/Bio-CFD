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

    testsuite = [&
      new_unittest("find_dist_node", test_find_dist_node), &
      new_unittest("shift_surface_nodes_initial", test_shift_surface_nodes_initial_simple) &
    ]
  end subroutine collect_search


  subroutine test_find_dist_node(error)
    use biocfd_search, only : findDistnode
    use global, only : block, blk_start, nblocks
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    integer(int64) :: expected

    ! Setup the required block variables
    if (allocated(block)) deallocate(block)
    blk_start = 1
    nblocks = 1
    allocate(block(nblocks))
    block(1)%ibNodes = 6
    block(1)%ibNodeId = [51, 51, 51, 51, 51, 51]
    block(1)%xnode = [0, -7, -4, 1, -10, 8]
    block(1)%ynode = [0, -7, -4, 1, -10, 8]
    block(1)%znode = [0, -7, -4, 1, -10, 8]

    call findDistnode
    expected = 6
    call check(error, block(1)%mk, expected)
    if (allocated(error)) return
  end subroutine test_find_dist_node

  !> Test the shiftSurfaceNodesInitial subroutine but only for the
  !> branch where the node ID is not equal to 51 or 52. This is the
  !> case where the expected output just becomes e.g. xnode + xshift
  subroutine test_shift_surface_nodes_initial_simple(error)
    use global, only : block, blk_start, nblocks
    use biocfd_search, only : shiftSurfaceNodesInitial
    type(error_type), allocatable, intent(out) :: error
    integer :: i

    ! Setup the required block variables
    if (allocated(block)) deallocate(block)
    blk_start = 1
    nblocks = 1
    allocate(block(nblocks))
    block(1)%ibNodes = 1
    ! As mentioned above, when node ID is 50, the code branch is
    ! relatively simple
    block(1)%ibNodeId = [50]
    block(1)%xnode = [1]
    block(1)%ynode = [2]
    block(1)%znode = [3]

    ! Set block shift up
    block(1)%xshift = 1
    block(1)%yshift = 2
    block(1)%zshift = 3

    call shiftSurfaceNodesInitial

    call check(error, all(block(1)%xnode1 == block(1)%xnode + block(1)%xshift))
    if (allocated(error)) return

    call check(error, all(block(1)%ynode1 == block(1)%ynode + block(1)%yshift))
    if (allocated(error)) return

    call check(error, all(block(1)%znode1 == block(1)%znode + block(1)%zshift))
    if (allocated(error)) return

  end subroutine test_shift_surface_nodes_initial_simple
end module test_search

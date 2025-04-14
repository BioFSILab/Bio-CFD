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
      new_unittest("shift_surface_nodes_initial", test_shift_surface_nodes_initial_simple), &
      new_unittest("compute_surface_norm", test_compute_surface_norm) &
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

    do i=1, size(block(1)%ibNodeId)
      call check(error, block(1)%xnode1(i), block(1)%xnode(i) + block(1)%xshift)
      if (allocated(error)) return

      call check(error, block(1)%ynode1(i), block(1)%ynode(i) + block(1)%yshift)
      if (allocated(error)) return

      call check(error, block(1)%znode1(i), block(1)%znode(i) + block(1)%zshift)
      if (allocated(error)) return
    end do

  end subroutine test_shift_surface_nodes_initial_simple

  !> Tests the computeSurfaceNorm subroutine, which computes the
  !> surface normals for triangles formed by three points
  subroutine test_compute_surface_norm(error)
    use global, only : block, blk_start, nblocks, inor
    use biocfd_search, only : computeSurfaceNorm
    type(error_type), allocatable, intent(out) :: error

   ! Setup the required block variables
    if (allocated(block)) deallocate(block)
    blk_start = 1
    nblocks = 1
    allocate(block(nblocks))
    block(1)%ibElems = 1

    block(1)%ibElP1 = [1]
    block(1)%ibElP2 = [2]
    block(1)%ibElP3 = [3]

    ! Make a triangle of points in the z=0 plane
    block(1)%xnode1 = [0, 1, 0]
    block(1)%ynode1 = [0, 0, 1]
    block(1)%znode1 = [0, 0, 0]

    inor = 1
    call computeSurfaceNorm

    ! The avg point should be [1/3, 1/3, 0]
    call check(error, block(1)%xcent(1), 1._dp/3)
    if (allocated(error)) return
    call check(error, block(1)%ycent(1), 1._dp/3)
    if (allocated(error)) return
    call check(error, block(1)%zcent(1), 0._dp)
    if (allocated(error)) return

    ! Cross-product of input vectors should give 0, 0, +1
    call check(error, block(1)%cosAlpha(1), 0._dp)
    if (allocated(error)) return
    call check(error, block(1)%cosBeta(1), 0._dp)
    if (allocated(error)) return
    call check(error, block(1)%cosGamma(1), 1._dp)
    if (allocated(error)) return

    deallocate(&
      block(1)%xcent, &
      block(1)%ycent, &
      block(1)%zcent, &
      block(1)%cosAlpha, &
      block(1)%cosBeta, &
      block(1)%cosGamma, &
      block(1)%alpha3, &
      block(1)%beta3, &
      block(1)%gamma3)

    ! Now swap x and y
    block(1)%xnode1 = [0, 0, 1]
    block(1)%ynode1 = [0, 1, 0]
    block(1)%znode1 = [0, 0, 0]

    call computeSurfaceNorm

    ! The avg point should be [1/3, 1/3, 0]
    call check(error, block(1)%xcent(1), 1._dp/3)
    if (allocated(error)) return
    call check(error, block(1)%ycent(1), 1._dp/3)
    if (allocated(error)) return
    call check(error, block(1)%zcent(1), 0._dp)
    if (allocated(error)) return

    ! Cross-product of input vectors should now give 0, 0, -1
    call check(error, block(1)%cosAlpha(1), 0._dp)
    if (allocated(error)) return
    call check(error, block(1)%cosBeta(1), 0._dp)
    if (allocated(error)) return
    call check(error, block(1)%cosGamma(1), -1._dp)
    if (allocated(error)) return


  end subroutine test_compute_surface_norm
end module test_search

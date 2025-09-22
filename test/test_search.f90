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

    testsuite = [ &
         new_unittest("find_dist_node", test_find_dist_node), &
         new_unittest("cellCount_solid_coarse", test_cell_count_solid_coarse) &
         ]
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

  subroutine test_cell_count_solid_coarse(error)
    use biocfd_blocks, only: Blocks
    use biocfd_search, only: cellCount_solid_coarse
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    type(Blocks) :: blk
    integer(int64) :: expected

    blk%nx = 5
    blk%ny = 5
    blk%nz = 5

    ! See biocfd_allocate_arrays (this is just a simpler way)
    allocate(blk%cell(blk%nx+2, blk%ny+2, blk%nz+2))
    blk%cell = 0

    call cellCount_solid_coarse(blk)

    expected = blk%nx * blk%ny * blk%nz
    ! Check that all cells are identified as fluid cells
    call check(error, blk%fluidCellCount, expected)
    if (allocated(error)) return
    ! Check that the sum of red + black cells == fluid cells
    call check(error, blk%redCellCount + blk%blackCellCount, blk%fluidCellCount)
    if (allocated(error)) return

  end subroutine test_cell_count_solid_coarse

end module test_search

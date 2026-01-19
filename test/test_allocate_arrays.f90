module test_allocate
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use testdrive, only : error_type, unittest_type, new_unittest, check
  implicit none
  private

  public :: collect_allocate_arrays

contains

  !> Collect all exported unit tests
  subroutine collect_allocate_arrays(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("allocate_arrays", test_allocate_arrays) &
         ]
  end subroutine collect_allocate_arrays

  subroutine test_allocate_arrays(error)
    use biocfd_allocate_arrays, only: allocateArrays
    use biocfd_block_type, only: Block_t
    ! Error handling
    type(error_type), allocatable, intent(out) :: error

    type(Block_t) :: blk

    blk%nx = 5
    blk%ny = 5
    blk%nz = 5
    blk%ibElems = 10

    call allocateArrays(blk)

    ! Check that several (but not all) of the arrays we expect to be
    ! allocated are
    call check(error, allocated(blk%u))
    if (allocated(error)) return
    call check(error, allocated(blk%ut))
    if (allocated(error)) return
    call check(error, allocated(blk%v))
    if (allocated(error)) return
    call check(error, allocated(blk%vt))
    if (allocated(error)) return
    call check(error, allocated(blk%w))
    if (allocated(error)) return
    call check(error, allocated(blk%wt))
    if (allocated(error)) return


  end subroutine test_allocate_arrays



end module test_allocate

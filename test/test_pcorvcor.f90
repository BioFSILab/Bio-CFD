module test_pcor_vcor
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use testdrive, only : error_type, unittest_type, new_unittest, check
  implicit none
  private

  public :: collect_pcor_vcor

contains

  !> Collect all exported unit tests
  subroutine collect_pcor_vcor(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("test_update_velocity_newv", test_update_velocity_newv) &
         ]
  end subroutine collect_pcor_vcor


  !> Test that when we call updateVelocity_newv the contents of u, v,
  !> and w, are correctly copied to ut, vt, and wt, when move_check is
  !> 1
  subroutine test_update_velocity_newv(error)
    use biocfd_blocks, only: Blocks
    use biocfd_pcor_vcor, only: updateVelocity_newv
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    type(Blocks) :: blk
    integer(int64) :: expected

    blk%nx = 5
    blk%ny = 5
    blk%nz = 5

    ! See biocfd_allocate_arrays (this is just a simpler way)
    allocate(blk%u(blk%nx+2, blk%ny+2, blk%nz+2))
    allocate(blk%v(blk%nx+2, blk%ny+2, blk%nz+2))
    allocate(blk%w(blk%nx+2, blk%ny+2, blk%nz+2))

    allocate(blk%ut(blk%nx+2, blk%ny+2, blk%nz+2))
    allocate(blk%vt(blk%nx+2, blk%ny+2, blk%nz+2))
    allocate(blk%wt(blk%nx+2, blk%ny+2, blk%nz+2))

    blk%u = 1
    blk%v = 2
    blk%w = 3

    blk%ut = 4
    blk%vt = 5
    blk%wt = 6

    ! Set move check to 1 so that the subroutine does something
    blk%move_check = 1

    ! Before calling updateVelocity_newv the arrays should be different
    call check(error, any(blk%u == blk%ut), .false.)
    if (allocated(error)) return
    call check(error, any(blk%v == blk%vt), .false.)
    if (allocated(error)) return
    call check(error, any(blk%w == blk%wt), .false.)
    if (allocated(error)) return

    call updateVelocity_newv(blk)

    ! After calling updatingVelocity_newv the arrays should be the same
    call check(error, all(blk%u == blk%ut), .true.)
    if (allocated(error)) return
    call check(error, all(blk%v == blk%vt), .true.)
    if (allocated(error)) return
    call check(error, all(blk%w == blk%wt), .true.)
    if (allocated(error)) return

  end subroutine test_update_velocity_newv

end module test_pcor_vcor

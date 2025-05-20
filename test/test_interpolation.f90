module test_interpolation
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use testdrive, only : error_type, unittest_type, new_unittest, check
  implicit none
  private

  public :: collect_interpolation

contains

  !> Collect all exported unit tests
  subroutine collect_interpolation(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("linear_interpolation", test_linear_interpolation), &
                 new_unittest("bilinear_interpolation", test_bilinear_interpolation)]
  end subroutine collect_interpolation

  subroutine test_linear_interpolation(error)
    use biocfd_interpolation, only: linear_interpolation
    ! Error handling
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: result

    result = linear_interpolation(1.5_dp, 2._dp, 1._dp, 2._dp, 1._dp)
    call check(error, result, 1.5_dp)
    if (allocated(error)) return

    ! Actually extrapolation
    result = linear_interpolation(3._dp, 2._dp, 1._dp, 2._dp, 1._dp)
    call check(error, result, 3._dp)
    if (allocated(error)) return

    ! Flat distribution
    result = linear_interpolation(50._dp, 100._dp, 1._dp, 1._dp, 1._dp)
    call check(error, result, 1._dp)
    if (allocated(error)) return

    ! Negative slope
    result = linear_interpolation(0.5_dp, 1._dp, 0._dp, 0._dp, 10._dp)
    call check(error, result, 5._dp)
    if (allocated(error)) return

  end subroutine test_linear_interpolation


  subroutine test_bilinear_interpolation(error)
    use biocfd_interpolation, only: bilinear_interpolation
    ! Error handling
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: result

    ! Flat distribution
    result = bilinear_interpolation(0.5_dp, 0.5_dp, 1._dp, 0._dp, 1._dp, 0._dp, &
                                    [1._dp, 1._dp, 1._dp, 1._dp])
    call check(error, result, 1._dp)
    if (allocated(error)) return

     ! Actually extrapolation
    result = bilinear_interpolation(2._dp, 0.5_dp, 1._dp, 0._dp, 1._dp, 0._dp, &
                                    [1._dp, 2._dp, 3._dp, 4._dp])
    call check(error, result, 4._dp)
    if (allocated(error)) return

  end subroutine test_bilinear_interpolation

end module test_interpolation

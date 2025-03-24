module test_demo
  use testdrive, only : error_type, unittest_type, new_unittest, check
  implicit none
  private

  public :: collect_demo

contains

  !> Collect all exported unit tests
  subroutine collect_demo(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("substitute", test_substitute)]
  end subroutine collect_demo

  !> Check substitution of a single line
  subroutine test_substitute(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    call check(error, "This is a valid example", "This is a valid example")
  end subroutine test_substitute
end module test_demo

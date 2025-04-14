module test_search
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
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    call check(error, "This is a valid example", "This is a valid example")
  end subroutine test_find_dist_node
end module test_search

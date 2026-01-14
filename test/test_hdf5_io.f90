#if USE_HDF5 == 1
module test_hdf5_io
  use, intrinsic :: iso_fortran_env, only: int64, dp =>real64
  use testdrive, only : error_type, unittest_type, new_unittest, check
  use biocfd_hdf5_io, only: hdf5_write_real, hdf5_write_int, hdf5_read_real_3d, &
       hdf5_read_real_2d, hdf5_read_real_1d, hdf5_read_int_3d, hdf5_read_int_2d, &
       hdf5_read_int_1d, hdf5_read_real_scalar, hdf5_read_int_scalar
  use hdf5, only: hid_t, hsize_t, h5f_acc_rdonly_f, h5t_std_i64le, &
       h5t_native_double, h5aread_f, h5dread_f, h5aclose_f, h5aopen_name_f, &
       h5dclose_f, h5close_f, h5dget_space_f, h5dopen_f, h5fclose_f, &
       h5fopen_f, h5gclose_f, h5gopen_f, h5open_f, h5sclose_f, &
       h5sget_simple_extent_dims_f
  implicit none
  private

  public :: collect_hdf5

contains

  !> Collect all exported unit tests
  subroutine collect_hdf5(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("hdf5_write_real_3d", test_hdf5_write_real_3d), &
         new_unittest("hdf5_write_real_2d", test_hdf5_write_real_2d), &
         new_unittest("hdf5_write_real_1d", test_hdf5_write_real_1d), &
         new_unittest("hdf5_write_real_scalar", test_hdf5_write_real_scalar), &
         new_unittest("hdf5_write_int_3d", test_hdf5_write_int_3d), &
         new_unittest("hdf5_write_int_2d", test_hdf5_write_int_2d), &
         new_unittest("hdf5_write_int_1d", test_hdf5_write_int_1d), &
         new_unittest("hdf5_write_int_scalar", test_hdf5_write_int_scalar) &
         ]
  end subroutine collect_hdf5

  subroutine test_hdf5_write_real_3d(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j, k
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "3d_real_array"

    real(dp), allocatable :: array3d_read(:,:,:)
    real(dp), allocatable :: array3d(:,:,:)

    logical :: arrays_equal

    allocate(array3d(4,3,2))
    do k = 1, 2
      do j = 1, 3
        do i = 1, 4
           array3d(i,j,k) = real(i + 10*j + 100*k, dp)
        end do
      end do
    end do

    call hdf5_write_real(filename=filename, array_input_3d=array3d, key=key, group=group)

    call hdf5_read_real_3d(filename=filename, group=group, key=key, output=array3d_read)

    arrays_equal = all(array3d_read == array3d)

    call check(error, .true. , arrays_equal)
    if (allocated(error)) return
  end subroutine test_hdf5_write_real_3d

  subroutine test_hdf5_write_real_2d(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    integer :: j, k
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "2d_real_array"

    real(dp), allocatable :: array2d_read(:,:)
    real(dp), allocatable :: array2d(:,:)

    logical :: arrays_equal

    allocate(array2d(3,2))
    do k = 1, 2
      do j = 1, 3
           array2d(j,k) = real(10*j + 100*k, dp)
      end do
    end do

    call hdf5_write_real(filename=filename, array_input_2d=array2d, key=key, group=group)

    call hdf5_read_real_2d(filename=filename, group=group, key=key, output=array2d_read)

    arrays_equal = all(array2d_read == array2d)

    call check(error, .true. , arrays_equal)
    if (allocated(error)) return
  end subroutine test_hdf5_write_real_2d

  subroutine test_hdf5_write_real_1d(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    integer :: k
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "1d_real_array"

    real(dp), allocatable :: array1d_read(:)
    real(dp), allocatable :: array1d(:)

    logical :: arrays_equal

    allocate(array1d(2))
    do k = 1, 2
           array1d(k) = real(100*k, dp)
    end do

    call hdf5_write_real(filename=filename, array_input_1d=array1d, key=key, group=group)

    call hdf5_read_real_1d(filename=filename, group=group, key=key, output=array1d_read)

    arrays_equal = all(array1d_read == array1d)

    call check(error, .true. , arrays_equal)
    if (allocated(error)) return
  end subroutine test_hdf5_write_real_1d

  subroutine test_hdf5_write_real_scalar(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: adims
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "real_scalar"

    real(dp) :: scalar_read
    real(dp) :: scalar

    adims=1
    scalar = 20.5_dp

    call hdf5_write_real(filename=filename, scalar_input=scalar, key=key, group=group)

    call hdf5_read_real_scalar(filename=filename, group=group, key=key, output=scalar_read)

    call check(error, scalar , scalar_read)
    if (allocated(error)) return
  end subroutine test_hdf5_write_real_scalar

  subroutine test_hdf5_write_int_3d(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j, k
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "3d_integer_array"

    integer(int64), allocatable :: array3d_read(:,:,:)
    integer(int64), allocatable :: array3d(:,:,:)

    logical :: arrays_equal

    allocate(array3d(4,3,2))
    do k = 1, 2
      do j = 1, 3
        do i = 1, 4
           array3d(i,j,k) = i + 10*j + 100*k
        end do
      end do
    end do

    call hdf5_write_int(filename=filename, array_input_3d=array3d, key=key, group=group)

    call hdf5_read_int_3d(filename=filename, group=group, key=key, output=array3d_read)

    arrays_equal = all(array3d_read == array3d)

    call check(error, .true. , arrays_equal)
    if (allocated(error)) return
  end subroutine test_hdf5_write_int_3d

  subroutine test_hdf5_write_int_2d(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    integer :: j, k
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "2d_integer_array"

    integer(int64), allocatable :: array2d_read(:,:)
    integer(int64), allocatable :: array2d(:,:)

    logical :: arrays_equal

    allocate(array2d(3,2))
    do k = 1, 2
      do j = 1, 3
           array2d(j,k) = 10*j + 100*k
      end do
    end do

    call hdf5_write_int(filename=filename, array_input_2d=array2d, key=key, group=group)

    call hdf5_read_int_2d(filename=filename, group=group, key=key, output=array2d_read)

    arrays_equal = all(array2d_read == array2d)

    call check(error, .true. , arrays_equal)
    if (allocated(error)) return
  end subroutine test_hdf5_write_int_2d

  subroutine test_hdf5_write_int_1d(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    integer :: k
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "1d_integer_array"

    integer(int64), allocatable :: array1d_read(:)
    integer(int64), allocatable :: array1d(:)

    logical :: arrays_equal

    allocate(array1d(2))
    do k = 1, 2
           array1d(k) = 100*k
    end do

    call hdf5_write_int(filename=filename, array_input_1d=array1d, key=key, group=group)

    call hdf5_read_int_1d(filename=filename, group=group, key=key, output=array1d_read)

    arrays_equal = all(array1d_read == array1d)

    call check(error, .true. , arrays_equal)
    if (allocated(error)) return
  end subroutine test_hdf5_write_int_1d

  subroutine test_hdf5_write_int_scalar(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: adims
    character(len=*), parameter :: filename = "test_output.h5"
    character(len=*), parameter :: group = "/test_group"
    character(len=*), parameter :: key = "integer_scalar"

    integer(int64) :: scalar_read
    integer(int64) :: scalar

    adims=1
    scalar = 20

    call hdf5_write_int(filename=filename, scalar_input=scalar, key=key, group=group)

    call hdf5_read_int_scalar(filename=filename, group=group, key=key, output=scalar_read)

    call check(error, scalar , scalar_read)
    if (allocated(error)) return
  end subroutine test_hdf5_write_int_scalar

end module test_hdf5_io
#endif

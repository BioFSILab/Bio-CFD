! FPM resolves which modules need to be compiled before any
! preprocessing so we need this guard here to make sure this file is
! effectively empty
#if USE_HDF5 == 1
module biocfd_hdf5_io
  use, intrinsic :: iso_fortran_env, only: dp => real64, int8,int32,int64
  use hdf5, only: hsize_t, hid_t, h5open_f, h5fopen_f, h5fcreate_f, h5lexists_f, &
       h5gcreate_f, h5gopen_f, h5screate_f, h5screate_simple_f, h5dcreate_f, &
       h5dwrite_f, h5dclose_f, h5sclose_f, h5gclose_f, h5fclose_f, h5close_f, &
       h5F_acc_rdwr_f, H5S_SCALAR_F, H5T_STD_I64LE, h5f_acc_excl_f, h5t_native_double, &
       h5acreate_f, h5awrite_f, size_t, h5aclose_f, h5sclose_f, H5F_ACC_RDONLY_F, &
       h5dopen_f, h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
       h5aopen_name_f, h5aread_f
  implicit none

  private

  public :: hdf5_write_real, hdf5_write_int, hdf5_read_real_3d, hdf5_read_real_2d, &
       hdf5_read_real_1d, hdf5_read_int_3d, hdf5_read_int_2d, &
       hdf5_read_int_1d, hdf5_read_real_scalar, hdf5_read_int_scalar
contains

  subroutine hdf5_read_real_3d(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer(hsize_t), dimension(3) :: dims, maxdims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    real(dp), allocatable, intent(out) :: output(:,:,:)

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5dopen_f(file_id, trim(group)//"/"//trim(key), dset_id, error_hdf5)

    call h5dget_space_f(dset_id, dspace_id, error_hdf5)

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error_hdf5)

    allocate(output(dims(1), dims(2), dims(3)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, output, dims, error_hdf5)

    call h5dclose_f(dset_id, error_hdf5)
    call h5sclose_f(dspace_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_real_3d

  subroutine hdf5_read_int_3d(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer(hsize_t), dimension(3) :: dims, maxdims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    integer(int64), allocatable, intent(out) :: output(:,:,:)

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5dopen_f(file_id, trim(group)//"/"//trim(key), dset_id, error_hdf5)

    call h5dget_space_f(dset_id, dspace_id, error_hdf5)

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error_hdf5)

    allocate(output(dims(1), dims(2), dims(3)))

    call h5dread_f(dset_id, H5T_STD_I64LE, output, dims, error_hdf5)

    call h5dclose_f(dset_id, error_hdf5)
    call h5sclose_f(dspace_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_int_3d

  subroutine hdf5_read_real_2d(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer(hsize_t), dimension(2) :: dims, maxdims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    real(dp), allocatable, intent(out) :: output(:,:)

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5dopen_f(file_id, trim(group)//"/"//trim(key), dset_id, error_hdf5)

    call h5dget_space_f(dset_id, dspace_id, error_hdf5)

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error_hdf5)

    allocate(output(dims(1), dims(2)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, output, dims, error_hdf5)

    call h5dclose_f(dset_id, error_hdf5)
    call h5sclose_f(dspace_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_real_2d

  subroutine hdf5_read_int_2d(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer(hsize_t), dimension(2) :: dims, maxdims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    integer(int64), allocatable, intent(out) :: output(:,:)

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5dopen_f(file_id, trim(group)//"/"//trim(key), dset_id, error_hdf5)

    call h5dget_space_f(dset_id, dspace_id, error_hdf5)

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error_hdf5)

    allocate(output(dims(1), dims(2)))

    call h5dread_f(dset_id, H5T_STD_I64LE, output, dims, error_hdf5)

    call h5dclose_f(dset_id, error_hdf5)
    call h5sclose_f(dspace_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_int_2d

  subroutine hdf5_read_real_1d(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer(hsize_t), dimension(1) :: dims, maxdims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    real(dp), allocatable, intent(out) :: output(:)

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5dopen_f(file_id, trim(group)//"/"//trim(key), dset_id, error_hdf5)

    call h5dget_space_f(dset_id, dspace_id, error_hdf5)

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error_hdf5)

    allocate(output(dims(1)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, output, dims, error_hdf5)

    call h5dclose_f(dset_id, error_hdf5)
    call h5sclose_f(dspace_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_real_1d

  subroutine hdf5_read_int_1d(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer(hsize_t), dimension(1) :: dims, maxdims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    integer(int64), allocatable, intent(out) :: output(:)

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5dopen_f(file_id, trim(group)//"/"//trim(key), dset_id, error_hdf5)

    call h5dget_space_f(dset_id, dspace_id, error_hdf5)

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error_hdf5)

    allocate(output(dims(1)))

    call h5dread_f(dset_id, H5T_STD_I64LE, output, dims, error_hdf5)

    call h5dclose_f(dset_id, error_hdf5)
    call h5sclose_f(dspace_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_int_1d

  subroutine hdf5_read_real_scalar(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, group_id, attr_id, attr_space_id
    INTEGER(HSIZE_T), DIMENSION(1) :: adims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    real(dp), intent(out) :: output

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5gopen_f(file_id, group, group_id, error_hdf5)

    call h5aopen_name_f(group_id, key, attr_id, error_hdf5)

    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, output, adims, error_hdf5)

    call h5aclose_f(attr_id, error_hdf5)
    call h5gclose_f(group_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_real_scalar

  subroutine hdf5_read_int_scalar(filename, group, key, output)
    integer :: error_hdf5
    integer(hid_t) :: file_id, group_id, attr_id, attr_space_id
    INTEGER(HSIZE_T), DIMENSION(1) :: adims
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: key

    integer(int64), intent(out) :: output

    call h5open_f(error_hdf5)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error_hdf5)

    call h5gopen_f(file_id, group, group_id, error_hdf5)

    call h5aopen_name_f(group_id, key, attr_id, error_hdf5)

    call h5aread_f(attr_id, H5T_STD_I64LE, output, adims, error_hdf5)

    call h5aclose_f(attr_id, error_hdf5)
    call h5gclose_f(group_id, error_hdf5)
    call h5fclose_f(file_id, error_hdf5)
    call h5close_f(error_hdf5)

  end subroutine hdf5_read_int_scalar

  subroutine hdf5_write_real(filename,scalar_input,&
                             array_input_1d,array_input_2d,array_input_3d,key,group)

    character(len=*), intent(in) :: filename
    real (dp), optional, intent(in) :: array_input_1d(:)
    real (dp), optional, intent(in) :: array_input_2d(:,:)
    real (dp), optional, intent(in) :: array_input_3d(:,:,:)
    real (dp), optional, intent (in) :: scalar_input
    character(len=*), intent(in) :: key, group
    character(len=20) :: dataset_location
    integer(hsize_t),allocatable :: data_dims(:)
    integer(hid_t) :: file_id, dspace_id, dset_id, group_id
    integer(hid_t) :: attr_id, aspace_id
    integer(size_t) :: adims(1)
    integer (int32) :: space_rank
    integer (int8) :: arguments_present
    integer :: error
    logical :: file_exists, dataset_exists, group_exists
    arguments_present=0
    if (present(scalar_input)) then
      arguments_present = arguments_present + 1
      adims = [0]
    end if

    if (present(array_input_1d)) then
      if(arguments_present /= 0) then
        error stop "More than one argument present hdf5_write"
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(1))
      data_dims(1) = size(array_input_1d,1)
      space_rank = 1
    end if

    if (present(array_input_2d)) then
      if(arguments_present /= 0) then
        error stop "More than one argument present hdf5_write"
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(2))
      data_dims(1) = size(array_input_2d,1)
      data_dims(2) = size(array_input_2d,2)
      space_rank = 2
    end if

    if (present(array_input_3d)) then
      if(arguments_present /= 0) then
        error stop "More than one argument present hdf5_write"
      end if
      allocate(data_dims(3))
      data_dims(1) = size(array_input_3d,1)
      data_dims(2) = size(array_input_3d,2)
      data_dims(3) = size(array_input_3d,3)
      space_rank = 3
    end if

    ! Initialize Fortran interface.
    call h5open_f(error)

    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
      ! Open an existing file.
      call h5fopen_f(trim(filename), h5F_acc_rdwr_f, file_id, error)
    else
      ! Create file requested
      call h5fcreate_f(trim(filename), h5f_acc_excl_f, file_id, error)
    end if

    ! Check if the group exists
    call h5lexists_f(file_id, group, group_exists, error)

    if (.not. group_exists) then
      ! Create a group
      call h5gcreate_f(file_id, group, group_id, error)
    else
      ! Open the existing group
      call h5gopen_f(file_id, group, group_id, error)
    end if

    if (present(scalar_input)) then
      ! Create a scalar dataspace
      call h5screate_f(H5S_SCALAR_F, aspace_id, error)
    else
      ! Open dataspace
      call h5screate_simple_f(space_rank,data_dims,dspace_id,error)
   end if
    dataset_location = group//"/"//key
    ! Check if the dataset exists
    call h5lexists_f(file_id, trim(dataset_location), dataset_exists, error)

    if (.not. dataset_exists) then
      if (present(scalar_input)) then
        ! Create attribute attached to the group
        call h5acreate_f(group_id, key, h5t_native_double, aspace_id, attr_id, error)
      else
        ! Create dataset if it doesn't exist already
         call h5dcreate_f(group_id,key,h5t_native_double,dspace_id,dset_id,error)
      end if
    end if

    if (present(scalar_input)) then
       ! Write to dataset
       call h5awrite_f(attr_id,h5t_native_double, scalar_input, adims, error)
    end if

    if (present(array_input_1d)) then
      ! Write to dataset
      call h5dwrite_f(dset_id,h5t_native_double,array_input_1d,data_dims,error)
    end if

    if (present(array_input_2d)) then
      ! Write to dataset
      call h5dwrite_f(dset_id,h5t_native_double,array_input_2d,data_dims,error)
    end if

    if (present(array_input_3d)) then
      ! Write to dataset
      call h5dwrite_f(dset_id,h5t_native_double,array_input_3d,data_dims,error)
    end if

    if (present(scalar_input)) then
       call h5aclose_f(attr_id, error)
       call h5sclose_f(aspace_id, error)
    else
       ! Close dataset
       call h5dclose_f(dset_id,error)

       ! Close dataspace
       call h5sclose_f(dspace_id, error)
    end if

    ! Close the group
    call h5gclose_f(group_id, error)

    ! Close the file.
    call h5fclose_f(file_id, error)

    ! Close Fortran interface.
    call h5close_f(error)

  end subroutine hdf5_write_real

  subroutine hdf5_write_int(filename,scalar_input,&
                            array_input_1d,array_input_2d,array_input_3d,key,group)

    character(len=*), intent(in) :: filename
    integer (int64), optional, intent(in) :: array_input_1d(:)
    integer (int64), optional, intent(in) :: array_input_2d(:,:)
    integer (int64), optional, intent(in) :: array_input_3d(:,:,:)
    integer (int64), optional, intent(in) :: scalar_input
    character(len=*), intent(in) :: key, group
    character(len=20) :: dataset_location
    integer(hsize_t),allocatable :: data_dims(:)
    integer(hid_t) :: file_id, dspace_id, dset_id, group_id
    integer(hid_t) :: attr_id, aspace_id
    integer(size_t) :: adims(1)
    integer (int32) :: space_rank
    integer (int8) :: arguments_present
    integer :: error
    logical :: file_exists, dataset_exists, group_exists
    arguments_present=0
    if (present(scalar_input)) then
      arguments_present = arguments_present + 1
      adims = [0]
    end if

    if (present(array_input_1d)) then
       if(arguments_present /= 0) then
        error stop "More than one argument present hdf5_write"
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(1))
      data_dims(1) = size(array_input_1d,1)
      space_rank = 1
    end if

    if (present(array_input_2d)) then
      if(arguments_present /= 0) then
        error stop "More than one argument present hdf5_write"
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(2))
      data_dims(1) = size(array_input_2d,1)
      data_dims(2) = size(array_input_2d,2)
      space_rank = 2
    end if

    if (present(array_input_3d)) then
      if(arguments_present /= 0) then
        error stop "More than one argument present hdf5_write"
      end if
      allocate(data_dims(3))
      data_dims(1) = size(array_input_3d,1)
      data_dims(2) = size(array_input_3d,2)
      data_dims(3) = size(array_input_3d,3)
      space_rank = 3
    end if

    ! Initialize Fortran interface.
    call h5open_f(error)

    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
      ! Open an existing file.
      call h5fopen_f(trim(filename), h5F_acc_rdwr_f, file_id, error)
    else
      ! Create file requested
      call h5fcreate_f(trim(filename), h5f_acc_excl_f, file_id, error)
    end if

    ! Check if the group exists
    call h5lexists_f(file_id, group, group_exists, error)

    if (.not. group_exists) then
      ! Create a group
      call h5gcreate_f(file_id, group, group_id, error)
    else
      ! Open the existing group
      call h5gopen_f(file_id, group, group_id, error)
   end if

    if (present(scalar_input)) then
      ! Create a scalar dataspace
      call h5screate_f(H5S_SCALAR_F, aspace_id, error)
    else
      ! Open dataspace
      call h5screate_simple_f(space_rank,data_dims,dspace_id,error)
    end if
     dataset_location = group//"/"//key
    ! Check if the dataset exists
     call h5lexists_f(file_id, trim(dataset_location), dataset_exists, error)

    if (.not. dataset_exists) then
      if (present(scalar_input)) then
        ! Create attribute attached to the group
        call h5acreate_f(group_id, key, H5T_STD_I64LE, aspace_id, attr_id, error)
      else
        ! Create dataset if it doesn't exist already
         call h5dcreate_f(group_id,key,H5T_STD_I64LE,dspace_id,dset_id,error)
      end if
    end if

    if (present(scalar_input)) then
       ! Write to dataset
       call h5awrite_f(attr_id,H5T_STD_I64LE, scalar_input, adims, error)
    end if

    if (present(array_input_1d)) then
      ! Write to dataset
      call h5dwrite_f(dset_id,H5T_STD_I64LE,array_input_1d,data_dims,error)
    end if

    if (present(array_input_2d)) then
      ! Write to dataset
      call h5dwrite_f(dset_id,H5T_STD_I64LE,array_input_2d,data_dims,error)
    end if

    if (present(array_input_3d)) then
      ! Write to dataset
       call h5dwrite_f(dset_id,H5T_STD_I64LE,array_input_3d,data_dims,error)
    end if

    if (present(scalar_input)) then
       call h5aclose_f(attr_id, error)
       call h5sclose_f(aspace_id, error)
    else
       ! Close dataset
       call h5dclose_f(dset_id,error)

       ! Close dataspace
       call h5sclose_f(dspace_id, error)
    end if

    ! Close the group
    call h5gclose_f(group_id, error)

    ! Close the file.
    call h5fclose_f(file_id, error)

    ! Close Fortran interface.
    call h5close_f(error)

  end subroutine hdf5_write_int
end module biocfd_hdf5_io
#endif

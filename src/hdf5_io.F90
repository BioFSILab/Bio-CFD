! FPM resolves which modules need to be compiled before any
! preprocessing so we need this guard here to make sure this file is
! effectively empty
#if USE_HDF5 == 1
module biocfd_hdf5_io
  use, intrinsic :: iso_fortran_env, only: dp => real64, int8,int16,int32,int64
  use hdf5, only: hsize_t, hid_t, h5open_f, h5fopen_f, h5fcreate_f, h5lexists_f, &
       h5gcreate_f, h5gopen_f, h5screate_f, h5screate_simple_f, h5dcreate_f, &
       h5dwrite_f, h5dclose_f, h5sclose_f, h5gclose_f, h5fclose_f, h5close_f, &
       h5F_acc_rdwr_f, H5S_SCALAR_F, H5T_STD_I64LE, h5f_acc_excl_f, h5t_native_double
  implicit none

  private

  public :: hdf5_write_real, hdf5_write_int
contains

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
    integer (int32) :: space_rank
    integer (int8) ::arguments_present
    integer :: error
    logical :: file_exists, dataset_exists, group_exists
    arguments_present=0
    if (present(scalar_input)) then
      arguments_present = arguments_present + 1
    end if

    if (present(array_input_1d)) then
      if(arguments_present /= 0 ) then
        error stop 'More than one argument present hdf5_write'
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(1))
      data_dims(1) = size(array_input_1d,1)
      space_rank = 1
    end if

    if (present(array_input_2d)) then
      if(arguments_present /= 0 ) then
        error stop 'More than one argument present hdf5_write'
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(2))
      data_dims(1) = size(array_input_2d,1)
      data_dims(2) = size(array_input_2d,2)
      space_rank = 2
    end if

    if (present(array_input_3d)) then
      if(arguments_present /= 0 ) then
        error stop 'More than one argument present hdf5_write'
      end if
      allocate(data_dims(3))
      data_dims(1) = size(array_input_3d,1)
      data_dims(2) = size(array_input_3d,2)
      data_dims(3) = size(array_input_3d,3)
      space_rank = 3
    end if

    ! Initialize Fortran interface.
    call h5open_f(error)

    inquire(file=trim(filename), exist=file_exists )

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
      call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    else
      ! Open dataspace
      call h5screate_simple_f(space_rank,data_dims,dspace_id,error)
   end if
    dataset_location = group//"/"//key
    ! Check if the dataset exists
    call h5lexists_f(file_id, trim(dataset_location), dataset_exists, error)

    if (.not. dataset_exists) then
      ! Create dataset if it doesn't exist already
      call h5dcreate_f(group_id,key,h5t_native_double,dspace_id,dset_id,error)
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

    ! Close dataset
    call h5dclose_f(dset_id,error)

    ! Close dataspace
    call h5sclose_f(dspace_id, error)

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
    integer (int32) :: space_rank
    integer (int8) :: arguments_present
    integer :: error
    logical :: file_exists, dataset_exists, group_exists
    arguments_present=0
    if (present(scalar_input)) then
      arguments_present = arguments_present + 1
    end if

    if (present(array_input_1d)) then
       if(arguments_present /= 0 ) then
        error stop 'More than one argument present hdf5_write'
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(1))
      data_dims(1) = size(array_input_1d,1)
      space_rank = 1
    end if

    if (present(array_input_2d)) then
      if(arguments_present /= 0 ) then
        error stop 'More than one argument present hdf5_write'
      end if
      arguments_present = arguments_present + 1
      allocate(data_dims(2))
      data_dims(1) = size(array_input_2d,1)
      data_dims(2) = size(array_input_2d,2)
      space_rank = 2
    end if

    if (present(array_input_3d)) then
      if(arguments_present /= 0 ) then
        error stop 'More than one argument present hdf5_write'
      end if
      allocate(data_dims(3))
      data_dims(1) = size(array_input_3d,1)
      data_dims(2) = size(array_input_3d,2)
      data_dims(3) = size(array_input_3d,3)
      space_rank = 3
    end if

    ! Initialize Fortran interface.
    call h5open_f(error)

    inquire(file=trim(filename), exist=file_exists )

    if (file_exists) then
      ! Open an existing file.
      call h5fopen_f(trim(filename), h5F_acc_rdwr_f, file_id, error)
    else
      ! Create file requested
      call h5fcreate_f(trim(filename), h5f_acc_rdwr_f, file_id, error)
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
      call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    else
      ! Open dataspace
      call h5screate_simple_f(space_rank,data_dims,dspace_id,error)
    end if
     dataset_location = group//"/"//key
    ! Check if the dataset exists
     call h5lexists_f(file_id, trim(dataset_location), dataset_exists, error)

    if (.not. dataset_exists) then
      ! Create dataset if it doesn't exist already
      call h5dcreate_f(group_id,key,H5T_STD_I64LE,dspace_id,dset_id,error)
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

    ! Close dataset
    call h5dclose_f(dset_id,error)

    ! Close dataspace
    call h5sclose_f(dspace_id, error)

    ! Close the group
    call h5gclose_f(group_id, error)

    ! Close the file.
    call h5fclose_f(file_id, error)

    ! Close Fortran interface.
    call h5close_f(error)

  end subroutine hdf5_write_int
end module biocfd_hdf5_io
#endif

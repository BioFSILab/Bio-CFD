module biocfd_write_output_corner1
  use, intrinsic :: iso_fortran_env, only: dp => real64, int8,int16,int32,int64
  ! allow(use-all) - TODO: Aim to fix this in the future
  use global
#if USE_HDF5 == 1
  use hdf5
#endif
  implicit none

  private

#if USE_HDF5 == 1
  public :: writeOutput1_hdf5, body_plot, writeresult
#else
  public :: writeOutput1, body_plot, writeresult
#endif

contains

#if USE_HDF5 == 1
  subroutine hdf5_write_real(filename,scalar_input,&
                             array_input_1d,array_input_2d,array_input_3d,key,group)

    character(len=*), intent(in) :: filename
    real (dp), allocatable, optional, intent(in) :: array_input_1d(:)
    real (dp), allocatable, optional, intent(in) :: array_input_2d(:,:)
    real (dp), allocatable, optional, intent(in) :: array_input_3d(:,:,:)
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
    integer (int64), allocatable, optional, intent(in) :: array_input_1d(:)
    integer (int64), allocatable, optional, intent(in) :: array_input_2d(:,:)
    integer (int64), allocatable, optional, intent(in) :: array_input_3d(:,:,:)
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
      call h5dcreate_f(group_id,key,H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
    end if

    if (present(array_input_1d)) then
      ! Write to dataset
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,array_input_1d,data_dims,error)
    end if

    if (present(array_input_2d)) then
      ! Write to dataset
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,array_input_2d,data_dims,error)
    end if

    if (present(array_input_3d)) then
      ! Write to dataset
       call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,array_input_3d,data_dims,error)
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

      SUBROUTINE writeOutput1_hdf5
       CHARACTER(len=150)  :: filename1
       INTEGER  :: k, i, j, g
       REAL (dp), allocatable :: u1(:,:,:), v1(:,:,:), w1(:,:,:)
       character (len=11) :: dummy_1
       character (len=5) ::dummy_2

         IF((mod(ita,200_int64) ==0 .or. ita <= 2 ))then
         !$acc serial
         do g=1,nblocks
            write(dummy_1,'(A6,I5.5)') 'block_',g
            write(dummy_2,'(I5.5)') ita
            allocate(u1(2:block(g)%nz+1,2:block(g)%ny+1,2:block(g)%nx+1),&
                     v1(2:block(g)%nz+1,2:block(g)%ny+1,2:block(g)%nx+1),&
                     w1(2:block(g)%nz+1,2:block(g)%ny+1,2:block(g)%nx+1))
            DO k = 2, block(g)%nz+1
            DO j = 2, block(g)%ny+1
            DO i = 2, block(g)%nx+1
               u1(k,j,i) = 0.5_dp*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))
               v1(k,j,i) = 0.5_dp*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))
               w1(k,j,i) = 0.5_dp*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))
            END DO
            END DO
            END DO
            filename1="test_"//trim(dummy_2)//".h5"
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=u1,key='u1',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=v1,key='v1',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=w1,key='w1',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_1d=block(g)%xp,key='xp',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_1d=block(g)%yp,key='yp',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_1d=block(g)%zp,key='zp',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                scalar_input=block(g)%nx,key='zonei',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                scalar_input=block(g)%nx,key='zonej',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                scalar_input=block(g)%nz,key='zonek',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=block(g)%p,key='p',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 scalar_input=totime,key='totime',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                array_input_3d=block(g)%cell,key='cell',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                array_input_3d=block(g)%cell_n,key='cell_n',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                array_input_3d=block(g)%cell_pr,key='cell_pr',group=dummy_1)
            deallocate(u1,v1,w1)
         end do
         !$acc end serial
        !$acc wait
         ENDIF
       END SUBROUTINE writeOutput1_hdf5
#else
      SUBROUTINE writeOutput1
       CHARACTER(len=150)  :: filename1
       INTEGER  :: k, i, j, g
       REAL (dp) :: u1, v1, w1

         IF((mod(ita,200_int64) ==0 .or. ita <= 2 ))then

           Do g=1,nblocks
       WRITE(filename1,1)char_f,ita,g,re,block(2)%dx,nblocks
1     FORMAT('out/',A3,'_butter_fielddata.',i9.9,'.',i3.3,'.',f7.1,'.',f8.6,'.',i3.3,".dat")
           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables="x","y","z","u","v","w","p","totime","cellid","cell_n","cell_pr"'
            WRITE(786,*) 'zone, ', 'i = ', block(g)%nx,' j = ', block(g)%ny, ' k = ', block(g)%nz

            DO k = 2, block(g)%nz+1
            DO j = 2, block(g)%ny+1
            DO i = 2, block(g)%nx+1
               u1 = 0.5_dp*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))
               v1 = 0.5_dp*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))
               w1 = 0.5_dp*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))
               WRITE(786,*) block(g)%xp(i), block(g)%yp(j), block(g)%zp(k), u1, v1, w1, &
                            block(g)%p(i,j,k), totime, block(g)%cell(i,j,k) , &
                            block(g)%cell_n(i,j,k) , block(g)%cell_pr(i,j,k)
            END DO
            END DO
            END DO
            CLOSE(786)
        end do
        !$acc wait
         ENDIF
      END SUBROUTINE writeOutput1
#endif
      SUBROUTINE writeResult
        INTEGER::  i, j, k,g
        CHARACTER(len=70)  :: filename1
        IF(mod(ita,500_int64)==0)THEN
           Do g=1,nblocks
           WRITE(filename1,22)char_f,g,re,block(2)%dx
 22          FORMAT('out/Chkpt/',A3,'_butter_chkpt.',i3.3,'.',f6.1,'.',f8.6,".dat")
        OPEN (1,FILE=filename1,FORM='formatted')
        DO k = 1, block(g)%nz+2
        DO j = 1, block(g)%ny+2
        DO i = 1, block(g)%nx+2
          WRITE(1,*) block(g)%u(i,j,k), block(g)%v(i,j,k), block(g)%w(i,j,k), &
        block(g)%p(i,j,k), totime, ita, ita1
       END DO
       END DO
       END DO
        CLOSE(22)
        END DO
        END IF
      END SUBROUTINE writeResult

         SUBROUTINE body_plot
         INTEGER(int64) :: inode, ielem, g
         CHARACTER(len=150) :: filename1

          DO g=1,nblocks
          if (ita == 1 )then
          WRITE(filename1,108)
  108     FORMAT("out/butterfly.dat")
          OPEN(UNIT=857,FILE=filename1,STATUS='unknown')
          WRITE(857,*) 'TITLE = "FEstressplot"'
          WRITE(857,*) 'VARIABLES= "x", "y", "z","bd_n"'
          WRITE(857,*) 'ZONE NODES= ',block(g)%ibNodes,',ELEMENTS=',block(g)%ibElems,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
          DO inode = 1, block(g)%ibNodes
            WRITE(857,*)block(g)%xnode1(inode),block(g)%ynode1(inode),block(g)%znode1(inode), 99
          END DO
          DO ielem = 1, block(g)%ibElems
            WRITE(857,*) block(g)%ibElP1(ielem), block(g)%ibElP2(ielem), block(g)%ibElP3(ielem)
          END DO
          CLOSE(857)
         END IF
          END DO
        end subroutine body_plot

end module biocfd_write_output_corner1

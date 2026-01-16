module biocfd_last_conditions
  use, intrinsic :: iso_fortran_env, only : dp => real64, int64
  use biocfd_block_type, only: Block_t
#if USE_HDF5 == 1
  use biocfd_hdf5_io, only: hdf5_read_real_3d, hdf5_read_real_scalar, &
       hdf5_read_int_scalar
#endif
  implicit none
  private
#if USE_HDF5 == 1
  public :: lastConditions
#endif
contains
#if USE_HDF5 == 1
  SUBROUTINE lastConditions(blk, id, totime,ita,ita1)

    !> The Block that we want to setup
    type(Block_t), intent(inout) :: blk
    real(dp), intent(out) :: totime
    INTEGER(int64), intent(out)   ::     ita, ita1
    !> The ID number of the block, usually 1 for the coarse block and
    !> >=2 for the fine blocks
    integer(int64), intent(in) :: id
    INTEGER ::  i, j, k
    CHARACTER(len=150) :: filename
    character(len=20) :: id_as_string
    write(id_as_string, '(I0)') id

    filename="input/Checkpoint/Checkpoint_"//trim(id_as_string)//".h5"

    call hdf5_read_real_3d(filename=filename, group="/", key="u", output=blk%u)
    call hdf5_read_real_3d(filename=filename, group="/", key="v", output=blk%v)
    call hdf5_read_real_3d(filename=filename, group="/", key="w", output=blk%w)
    call hdf5_read_real_3d(filename=filename, group="/", key="p", output=blk%p)
    call hdf5_read_real_scalar(filename=filename, group="/", key="totime", output=totime)
    call hdf5_read_int_scalar(filename=filename, group="/", key="ita", output=ita)
    call hdf5_read_int_scalar(filename=filename, group="/", key="ita1", output=ita1)

    blk%ut = blk%u
    blk%vt = blk%v
    blk%wt = blk%w

  END SUBROUTINE lastConditions
#endif
end module biocfd_last_conditions

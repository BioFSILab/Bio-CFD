!> This module is used for the initialisation and finalisation of MPI.
!> It is written in such a way that the subroutines do the right thing
!> regardless of whether the code is compiled with MPI support or not.
module biocfd_mpi_helpers
  use, intrinsic :: iso_fortran_env, only: error_unit
#ifdef BIOCFD_MPI
  use mpi_f08, only: MPI_Abort, MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD, MPI_Finalize, &
                     MPI_Init_Thread, MPI_IN_PLACE, MPI_INTEGER, MPI_THREAD_SERIALIZED
#endif
#ifdef _OPENACC
  use openacc, only: acc_device_default, acc_get_num_devices
#endif

  implicit none

  private

  public :: biocfd_init, biocfd_finalize, get_block_iteration_params

contains

  !> Initialize MPI and then work out steps that we will be performing
  !> over blocks with `get_block_iteration_params`.
  subroutine biocfd_init(nblocks, start, finish, step, rank)
    !> The total number of blocks we are simulating
    integer, intent(in) :: nblocks
    !> The start, finish (both end and stop are keywords) and step size we will use to loop over blocks
    integer, intent(out) :: start, finish, step
    !> The MPI rank that we are running on (0 in the case we aren't using MPI)
    integer, intent(out) :: rank
#ifdef BIOCFD_MPI
    ! For MPI threading considerations
    integer :: required, provided
    ! An error value to check when using MPI
    integer :: ierror

    required = MPI_THREAD_SERIALIZED
    call MPI_Init_Thread(required, provided, ierror)

    if (provided < required) then
       write(error_unit, *) "MPI does not provide the required threading support. Aborting!"
       call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    call MPI_Comm_size(MPI_COMM_WORLD, step, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    if (step > nblocks) then
       if (rank == 0) then
          write(error_unit, *) "You have launched the program with more processes than blocks."
          write(error_unit, *) "This will lead to wasted resources."
          write(error_unit, *) "Please resubmit the job with the number of processors <= ", nblocks
       end if
       call MPI_Finalize()
       stop 1
    end if
#endif
    call  get_block_iteration_params(nblocks, start, finish, step, rank)
  end subroutine biocfd_init

  !> Finalize is an no-op in the case that we aren't using MPI
  subroutine biocfd_finalize()
#ifdef BIOCFD_MPI
    integer :: ierror
    call MPI_Finalize(ierror)
#endif
  end subroutine biocfd_finalize

!> Work out how work will be shared across MPI nodes
subroutine get_block_iteration_params(nblocks, start, finish, step, rank)
   !> The total number of blocks we are simulating
    integer, intent(in) :: nblocks
    !> The start, finish (both end and stop are keywords) and step size we will use to loop over blocks
    integer, intent(out) :: start, finish, step
    !> The MPI rank that we are running on (0 in the case we aren't using MPI)
    integer, intent(out) :: rank
#ifdef BIOCFD_MPI
    ! An error value to check when using MPI
    integer :: ierror
    call MPI_Comm_size(MPI_COMM_WORLD, step, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    ! MPI ranks are zero indexed, we want to start looping from 1 in
    ! fortran
    start = rank + 1
    finish = nblocks
#else
    ! If we aren't using MPI, it is straight forward. Our loop goes
    ! over every block from 1 in steps of 1. Rank is set to 0.
    start = 1
    finish = nblocks
    step = 1
    rank = 0
#endif

end subroutine get_block_iteration_params

!> Work out how work will be shared across MPI nodes when we have
!> multiple GPUs! The idea here is that we will work out how many GPUs
!> our rank has and how many GPUs there are overall. We can then split
!> the blocks proportionally across nodes. Each rank must have at
!> least one GPU for this to work.
subroutine get_block_iteration_params_gpu(nblocks, start, finish, step, rank)
   !> The total number of blocks we are simulating
    integer, intent(in) :: nblocks
    !> The start, finish (both end and stop are keywords) and step size we will use to loop over blocks
    integer, intent(out) :: start, finish, step
    !> The MPI rank that we are running on (0 in the case we aren't using MPI)
    integer, intent(out) :: rank
#if defined(BIOCFD_MPI) && defined(_OPENACC)
    ! An error value to check when using MPI
    integer :: ierror

    integer :: world_size
    integer, allocatable :: rank_gpus(:), rank_blocks(:)
    integer :: total_gpus
    !> This is blocks per GPU computed as nblocks / total_gpus (real
    !> as in a floating point number)
    real :: real_blocks_per_gpu
    ! Looping variable
    integer:: i

    call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    ! Get an array of GPUs and fill accordingly
    allocate(rank_gpus(world_size))
    rank_gpus(rank + 1) = acc_get_num_devices(acc_device_default)

    ! Gather all GPUs together such that every rank has an array of
    ! GPU numbers
    call MPI_Allgather(MPI_IN_PLACE, 0, MPI_INTEGER, rank_gpus, 1, MPI_INTEGER, MPI_COMM_WORLD, & 
                       ierror)
    total_gpus = sum(rank_gpus)

    ! This is how many blocks we need to have per GPU
    real_blocks_per_gpu = real(nblocks) / total_gpus

    ! This is the number of blocks we are going to place on each rank
    allocate(rank_blocks(world_size))
    ! This is the rounded number of blocks on each rank
    rank_blocks = nint(rank_gpus * real_blocks_per_gpu)

    start = sum(rank_blocks(1:rank))
    finish = sum(rank_blocks(1:rank+1))
    step = 1
    print *, "Testing - start = ", start, "finish = ", finish
#else
    ! If we aren't using MPI, it is straight forward. Our loop goes
    ! over every block from 1 in steps of 1. Rank is set to 0.
    start = 1
    finish = nblocks
    step = 1
    rank = 0
#endif

end subroutine get_block_iteration_params_gpu

end module biocfd_mpi_helpers

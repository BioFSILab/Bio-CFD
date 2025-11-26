module biocfd_gpu_helpers
#ifdef _OPENACC
  use openacc, only: acc_device_default, acc_get_num_devices, acc_set_device_num
#endif
#ifdef _OPENMP
  use omp_lib, only: omp_get_thread_num
#endif
  implicit none
  private

  public :: set_gpu

contains
  !> Sets the GPU to use on this particular GPU thread. Should be
  !> called inside an OpenMP parallel region to distribute work across
  !> multiple GPUs.
  subroutine set_gpu()
    ! For controlling OpenMP
    integer :: omp_thread_num
    ! For controlling OpenACC
    integer :: acc_devices

    omp_thread_num = 0
    acc_devices = 0

#ifdef _OPENMP
    omp_thread_num = omp_get_thread_num()
#endif
#ifdef _OPENACC
    ! Not checked, but apparently in nvfortran the default
    ! resolves to the same as `acc_device_nvidia` (see
    ! https://docs.nvidia.com/hpc-sdk/compilers/openacc-gs/index.html#defaults)
    ! For gfortran this can be set at runtime with an environment
    ! variable, ACC_DEVICE_TYPE. It may or may not pick up a
    ! compatible GPU if it can find it.
    acc_devices = acc_get_num_devices(acc_device_default)
    call acc_set_device_num(mod(omp_thread_num, acc_devices), acc_device_default)
#endif

  end subroutine set_gpu

end module biocfd_gpu_helpers
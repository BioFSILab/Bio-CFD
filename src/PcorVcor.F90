module biocfd_pcor_vcor
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, deltat, epsi, omega, omega1, omega2, omega3, omega4, pcitamax, &
       amgxita, couptime, dfinish, dstart, ita, mstime, nblocks, solvertime, totaltime, &
       totime
#ifdef _OPENMP
  use omp_lib, only: omp_get_max_threads, omp_get_thread_num
#endif
#ifdef _OPENACC
  use openacc, only: acc_device_default, acc_get_num_devices, acc_set_device_num
#endif
  use biocfd_fine_interp_bound, only : fineUpdate_newv_bd, fineUpdate_bd, fineUpdate_pc_bd
  use biocfd_coarse_update, only : coarseUpdate_newv, coarseUpdate_pc, coarseUpdate
  use biocfd_boundary_conditions, only : velocityBC
#ifdef BIOCFD_MPI
        use mpi_f08
#endif
  implicit none
  private

  public :: poissonSolver, updateVelocity_newv

  contains

      SUBROUTINE poissonSolver

        INTEGER(int64) :: i, j,k, n, g
        REAL (dp)    :: max_derr1, max_derr2, max_div, max_derrStdSt
        REAL (dp)    :: er_dudt, er_dvdt, er_dwdt, err_ds
        INTEGER(int64) :: max_nIterPcor
        CHARACTER(len=160) :: filename1

        ! For controlling OpenMP
        integer :: omp_threads
        integer :: omp_thread_num
        ! For controlling OpenACC
        integer :: acc_devices

        integer :: rank, num_proc, start_block
#ifdef BIOCFD_MPI
        integer :: ierror
#endif

          max_derrStdst=0._dp
          max_derr1=0._dp
          max_derr2=0._dp
          max_nIterPcor=0
          max_div=0._dp
          err_ds=0._dp
          er_dudt=0.
          er_dvdt=0.
          er_dwdt=0.

          ! Dummy values for omp/acc variables
          omp_threads = 1
          omp_thread_num = 0
          acc_devices = 0

          ! These variables are used to control MPI execution
          rank = 0
          num_proc = 1

#ifdef _OPENMP
          omp_threads = min(omp_get_max_threads(), size(block))
#endif
#ifdef _OPENACC
         ! Not checked, but apparently in nvfortran the default
         ! resolves to the same as `acc_device_nvidia` (see
         ! https://docs.nvidia.com/hpc-sdk/compilers/openacc-gs/index.html#defaults)
         ! For gfortran this can be set at runtime with an environment
         ! variable, ACC_DEVICE_TYPE. It may or may not pick up a
         ! compatible GPU if it can find it.
         acc_devices = acc_get_num_devices(acc_device_default)
#endif
#ifdef BIOCFD_MPI
         call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)
         call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
#endif
        start_block = rank + 1

        do g=start_block, size(block), num_proc
        !$acc parallel loop gang vector collapse (3) default(present)
        DO k = 1, block(g)%nz+2
        DO j = 1, block(g)%ny+2
        DO i = 1, block(g)%nx+2
           block(g)%b(i,j,k)  = 0.
           block(g)%pc(i,j,k) = 0.
           block(g)%pco(i,j,k)= 0.
        END DO
        END DO
        END DO
        !$acc end parallel loop
        block(g)%nIterPcor=0
        block(g)%derr2  = 0._dp
        block(g)%derrStdSt=0._dp
        end do


        solverTime=0.
        CALL cpu_time(dStart)
        ! This is the first place we require MPI communication
        CALL fineUpdate_newv_bd
#ifdef BIOCFD_MPI
        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        call MPI_Finalize(ierror)
        stop
#endif
        CALL coarseUpdate_newv

     CALL cpu_time(dfinish)
        coupTime=coupTime + dfinish -dstart

        !$omp parallel num_threads(omp_threads) default(none) &
        !$omp& private(dStart, dfinish, amgxita, mstime, g) &
        !$omp& shared(nblocks, acc_devices) firstprivate(omp_thread_num)

#ifdef _OPENMP
        omp_thread_num = omp_get_thread_num()
#endif

#ifdef _OPENACC
        call acc_set_device_num(mod(omp_thread_num, acc_devices), acc_device_default)
#endif

        !$omp do
        DO g=1,nblocks
           CALL computeDiv(g)    !divergence vector
           ! Do not compute Red/Black here for block 1
           if (g /= 1)  CALL REDBLACKSOR_linear(g)

        end do
        !$omp end do

        !$omp single
        CALL coarseUpdate_pc
        g=1
        CALL cpu_time(dStart)
        CALL REDBLACKSOR_linear(g)
        CALL cpu_time(dfinish)
        call fineUpdate_pc_bd
        !$omp end single
        !$omp do
        DO g=2,nblocks
         CALL cpu_time(dStart)
         amgxita=0
         CALL REDBLACKSOR_linear(g)
         CALL cpu_time(dfinish)
         msTime = msTime + dfinish-dstart
        end do
        !$omp end do
        !$omp single
          CALL coarseUpdate_pc
        !$omp end single

       !$omp do
       DO g=1,nblocks
               CALL correctPressure(g)  !pressure correction
               CALL correctVelocity(g)  !velocity correction
        END DO
        !$omp end do
        !$omp end parallel

        do g=1,nblocks
         CALL velocityBC(block(g))      !correct velocity at boundaries
        end do

         DO g=1,nblocks
         err_ds=0.
        !$acc parallel loop gang vector firstprivate (deltat)   &
        !$acc private (i, j, k, er_dudt, er_dvdt, er_dwdt)               &
        !$acc default(present) reduction (max: err_ds)
         DO n = 1, block(g)%fluidCellCount
           i = block(g)%fluidIndexPtr(n, 1)
           j = block(g)%fluidIndexPtr(n, 2)
           k = block(g)%fluidIndexPtr(n, 3)
           er_dudt = dabs((block(g)%ut(i,j,k) - block(g)%u(i,j,k)))/deltat
           er_dvdt = dabs((block(g)%vt(i,j,k) - block(g)%v(i,j,k)))/deltat
           er_dwdt = dabs((block(g)%wt(i,j,k) - block(g)%w(i,j,k)))/deltat
           err_ds = dmax1(err_ds, er_dudt,er_dvdt, er_dwdt)
         ENDDO
         !$acc end parallel loop
           block(g)%derrStdSt = err_ds


         ENDDO



        DO i=1,nblocks
               if ( block(i)%derr2 >max_derr2)then
                  max_derr2=block(i)%derr2
               end if
              if ( block(i)%derr1 >max_derr1)then
                 max_derr1=block(i)%derr1
              end if
              if ( block(i)%derrStdSt >max_derrStdSt)then
                 max_derrStdSt=block(i)%derrStdSt
              end if
              if ( block(i)%nIterPcor >max_nIterPcor)then
                 max_nIterPcor=block(i)%nIterPcor
              end if
              totalTime=totime + totalTime
          end do

            WRITE(filename1,1)
 1          FORMAT('sphere_iter.dat')
         OPEN(111,FILE=filename1,POSITION='APPEND',STATUS='unknown')
         WRITE(111,126)   ita, block(1)%nIterPcor, block(2)%nIterPcor, omega1, omega2, solverTime
         WRITE(*,16) ita, max_nIterPcor, max_derr2, max_derrStdSt, totalTime
 126      FORMAT(' ',I8, 2I10, 2F6.2,F14.9)
 16      FORMAT(' ',I8, I10, 4E15.6)
         CLOSE(111)

         DO g=1,nblocks
        !$acc parallel loop gang vector default(present) collapse (3)
         DO k = 1, block(g)%nz+2
         DO j = 1, block(g)%ny+2
         DO i = 1, block(g)%nx+2
           block(g)%u(i,j,k) = block(g)%ut(i,j,k)
           block(g)%v(i,j,k) = block(g)%vt(i,j,k)
           block(g)%w(i,j,k) = block(g)%wt(i,j,k)
         END DO
         END DO
         END DO
        !$acc end parallel loop
         END DO
        CALL fineUpdate_bd
        CALL coarseUpdate
      END SUBROUTINE poissonSolver

      SUBROUTINE computeDiv(g)

         INTEGER :: n, i, j, k,gg, counter, nx_var, ny_var
         INTEGER(int64),INTENT(IN) ::g
         gg=g
         nx_var=block(g)%nx
         ny_var=block(g)%ny
        !$acc parallel loop gang vector private (i, j, k,counter)   &
        !$acc default(present)
         DO n = 1, block(gg)%fluidCellCount
           i = block(gg)%fluidIndexPtr(n, 1)
           j = block(gg)%fluidIndexPtr(n, 2)
           k = block(gg)%fluidIndexPtr(n, 3)
           counter =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
           block(gg)%b(i,j,k) = &
             (block(gg)%ut(i,j,k) - block(gg)%ut(i-1,j,k))/block(gg)%deltax(i) +  &
             (block(gg)%vt(i,j,k) - block(gg)%vt(i,j-1,k))/block(gg)%deltay(j) +  &
             (block(gg)%wt(i,j,k) - block(gg)%wt(i,j,k-1))/block(gg)%deltaz(k)

         END DO
        !$acc end parallel loop
      END SUBROUTINE computeDiv

      SUBROUTINE correctPressure(g)

         INTEGER(int64) :: n, i, j, k,gg
         INTEGER(int64),INTENT(IN) ::g
        gg=g
        !$acc parallel loop gang vector private (i, j, k)   &
        !$acc default(present)
         DO n = 1, block(gg)%fluidCellCount
            i = block(gg)%fluidIndexPtr(n, 1)
            j = block(gg)%fluidIndexPtr(n, 2)
            k = block(gg)%fluidIndexPtr(n, 3)
            block(gg)%p(i,j,k) = block(gg)%p(i,j,k) + block(gg)%pc(i,j,k)
        END DO
        !$acc end parallel loop
      END SUBROUTINE correctPressure

      SUBROUTINE correctVelocity(g)

         INTEGER(int64) :: n, i, j, k,gg
         INTEGER(int64),INTENT(IN) ::g
         gg=g

        !$acc parallel loop gang vector private (i, j, k) firstprivate (deltat) &
        !$acc default(present)
         DO 30 n = 1, block(gg)%fluidCellCount
            i = block(gg)%fluidIndexPtr(n, 1)
            j = block(gg)%fluidIndexPtr(n, 2)
            k = block(gg)%fluidIndexPtr(n, 3)

            block(gg)%ut(i,j,k) = block(gg)%ut(i,j,k) - &
               deltat/(0.5d0*(block(gg)%deltax(i+1)+block(gg)%deltax(i))) * &
               (block(gg)%pc(i+1,j,k)-block(gg)%pc(i,j,k))
            block(gg)%vt(i,j,k) = block(gg)%vt(i,j,k) - &
               deltat/(0.5d0*(block(gg)%deltay(j+1)+block(gg)%deltay(j))) * &
               (block(gg)%pc(i,j+1,k)-block(gg)%pc(i,j,k))
            block(gg)%wt(i,j,k) = block(gg)%wt(i,j,k) - &
               deltat/(0.5d0*(block(gg)%deltaz(k+1)+block(gg)%deltaz(k))) * &
               (block(gg)%pc(i,j,k+1)-block(gg)%pc(i,j,k))
 30      CONTINUE
         !$acc end parallel loop
      END SUBROUTINE correctVelocity

      SUBROUTINE REDBLACKSOR_linear(g)

         INTEGER(int64) :: n, i, j, k, gg, nx_var, ny_var,nz_var,nxy
         REAL (dp) :: errSum,var,derr4
         INTEGER(int64),INTENT(IN) ::g

         gg=g
            if (gg == 1)then
                 omega = omega1
           elseif (gg == 2) then
                 omega=omega2
           elseif (gg == 3) then
                 omega=omega3
           else
                omega =omega4
          endif

          nx_var=block(gg)%nx
          ny_var=block(gg)%ny
          nz_var=block(gg)%nz

        derr4=11111.
         nxy= nx_var * ny_var
         block(g)%derr2 = 0._dp
         errSum = 0._dp
 3       block(g)%nIterPcor=block(g)%nIterPcor+1

        var=0.

        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k)
         DO 10 n = 1, block(gg)%redCellCount
            i = block(gg)%redCellIndexPtr(n, 1)
            j = block(gg)%redCellIndexPtr(n, 2)
            k = block(gg)%redCellIndexPtr(n, 3)

            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
              -block(gg)%Acx(i-1,1)*block(gg)%pco(i-1,j,k) &
              -block(gg)%Acx(i-1,3)*block(gg)%pco(i+1,j,k) &
              -block(gg)%Acy(j-1,1)*block(gg)%pco(i,j-1,k) &
              -block(gg)%Acy(j-1,3)*block(gg)%pco(i,j+1,k) &
              -block(gg)%Acz(k-1,1)*block(gg)%pco(i,j,k-1) &
              -block(gg)%Acz(k-1,3)*block(gg)%pco(i,j,k+1)) / &
              (block(gg)%Acx(i-1,2)+block(gg)%Acy(j-1,2)+block(gg)%Acz(k-1,2))
            block(gg)%pc(i,j,k) = (1._dp-omega)*block(gg)%pco(i,j,k) +omega*block(gg)%pc(i,j,k)

 10      CONTINUE
        !$acc end parallel loop

        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k)
         DO 20 n = 1, block(gg)%blackCellCount
             i = block(gg)%blackCellIndexPtr(n, 1)
             j = block(gg)%blackCellIndexPtr(n, 2)
             k = block(gg)%blackCellIndexPtr(n, 3)

            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
              -block(gg)%Acx(i-1,1)*block(gg)%pc(i-1,j,k) &
              -block(gg)%Acx(i-1,3)*block(gg)%pc(i+1,j,k) &
              -block(gg)%Acy(j-1,1)*block(gg)%pc(i,j-1,k) &
              -block(gg)%Acy(j-1,3)*block(gg)%pc(i,j+1,k) &
              -block(gg)%Acz(k-1,1)*block(gg)%pc(i,j,k-1) &
              -block(gg)%Acz(k-1,3)*block(gg)%pc(i,j,k+1)) / &
              (block(gg)%Acx(i-1,2)+block(gg)%Acy(j-1,2)+block(gg)%Acz(k-1,2))
            block(gg)%pc(i,j,k) = (1._dp-omega)*block(gg)%pco(i,j,k) +omega*block(gg)%pc(i,j,k)

 20      CONTINUE
        !$acc end parallel loop

       if(mod(block(gg)%nIterPcor,5_int64) ==0)then
       derr4=0.
        !$acc parallel loop gang vector reduction(max:derr4) default(present) private (i, j, k, var)
         DO 30 n = 1, block(gg)%fluidCellCount
             i = block(gg)%fluidIndexPtr(n, 1)
             j = block(gg)%fluidIndexPtr(n, 2)
             k = block(gg)%fluidIndexPtr(n, 3)
            var = abs(block(gg)%pc(i,j,k)-block(gg)%pco(i,j,k))
            derr4=dmax1(derr4,var)
 30      CONTINUE
        !$acc end parallel loop
        end if
        !$acc parallel loop gang vector collapse(3) default(present) private (i, j, k)
        DO k = 1, block(gg)%nz+2
        DO j = 1, block(gg)%ny+2
        DO i = 1, block(gg)%nx+2
            block(gg)%pco(i,j,k) = block(gg)%pc(i,j,k)

        END DO
        END DO
        END DO
        !$acc end parallel loop

         block(g)%derr2=derr4
         IF (derr4>=epsi .and. block(g)%nIterPcor <= pcItaMax) GOTO 3
      END SUBROUTINE REDBLACKSOR_linear

      SUBROUTINE updateVelocity_newv(g)

        INTEGER ::  i, j, k
         INTEGER (int64), INTENT(IN) :: g
        if (block(g)%move_check == 1) then
        !$acc parallel loop gang vector collapse(2) default(present)
        DO k = 1, block(g)%nz+2
        DO j = 1, block(g)%ny+2
        DO i = 1, block(g)%nx+2
           block(g)%ut(i,j,k) = block(g)%u(i,j,k)
           block(g)%vt(i,j,k) = block(g)%v(i,j,k)
           block(g)%wt(i,j,k) = block(g)%w(i,j,k)
       END DO
       END DO
       END DO
       !$acc end parallel loop
        endif

      END SUBROUTINE updateVelocity_newv
end module biocfd_pcor_vcor

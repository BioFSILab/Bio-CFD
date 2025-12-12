PROGRAM main
        use, intrinsic :: iso_fortran_env, only: int64, dp => real64
        USE global, only: block, coarse_flcnt_check, deltat, &
             ita, ita1, totaltime, totime, &
             pi,uc,re,intfr
        use biocfd_search, only: findDistnode, shiftSurfaceNodesInitial, computeSurfaceNorm, &
             tagging_th, tagging_th_move, block_move_check, cellcount_solid, &
             cellcount_solid_coarse, cellcount_solid_coarse_mv, change_block_coords, &
             change_block_interface, computenormdistance, computesurfacevariables, findtscells, &
             fine_block_cell, selectiveretagging_th, change_block_coords_interfaces
        use biocfd_pcor_vcor, only: poissonSolver, updateVelocity_newv
        use biocfd_boundary_conditions, only: velocityBC, solidcellbc, solidcellbc_move
        use biocfd_read_input, only: readInput, readBlockInterface, readSurfaceMeshGmsh
        use biocfd_allocate_arrays, only: allocateArrays
        use biocfd_interface_detail, only: interfaceDetail
        use biocfd_initial_conditions, only: initialConditions
        use biocfd_last_conditions, only: lastConditions
        use biocfd_coefficient_matrix, only: coefficientMatrix
        use biocfd_navier_stokes, only: non_uni_coeff, nsmomentum2order
        use biocfd_write_output_corner1, only: body_plot, writeresult
#if USE_HDF5 == 1
        use biocfd_write_output_corner1, only: write_output_hdf5
#else
        use biocfd_write_output_corner1, only: write_output_ascii
#endif
        use biocfd_forcing, only: pressureForcing1, pressureforcingfield, pressureforcingghost, &
             velocityforcing1, velocityforcingfield, velocityforcingghost
        use biocfd_mpi_helpers, only: biocfd_init, biocfd_finalize
        use biocfd_gpu_helpers, only: set_gpu
#ifdef BIOCFD_MPI
        use mpi_f08, only: MPI_Allreduce, MPI_Bcast, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, &
                           MPI_IN_PLACE, MPI_INTEGER8, MPI_Max
#endif
        IMPLICIT NONE

        INTEGER (int64) :: g
        INTEGER (int64)   :: surGeoPoints
        CHARACTER (LEN = 3)   :: char_f
        INTEGER               :: istart
        INTEGER (int64)   :: itamax, pcItaMax
        real(dp) :: aoa, aoa1, aoa2, phase_angle, piv_pt, mu_f, rho_f
        integer :: start, finish, step, rank

        CALL readInput(surGeoPoints, char_f, istart, itamax, pcItaMax, aoa, phase_angle, piv_pt, &
                       mu_f, rho_f)
        ! Need to init after we know the size of block
        call biocfd_init(size(block), start, finish, step, rank)
        CALL readBlockInterface
        CALL interfaceDetail

        ! Block 1, which is assumed to be the coarse block is
        ! allocated on every rank in MPI
        CALL allocateArrays(block(1))
        do g=start, finish, step
          if (g /= 1) then
            CALL readSurfaceMeshGmsh(block(g), surGeoPoints)
            CALL allocateArrays(block(g))
          end if
        end do

        phase_angle = phase_angle*pi/180_dp
        aoa1 = aoa*pi/180_dp
        aoa2 = -aoa1

        do g=start, finish, step
          if (g /= 1) then
            CALL findDistnode(block(g))
            CALL shiftSurfaceNodesInitial(block(g),aoa1,aoa2,piv_pt)
            CALL computeSurfaceNorm(block(g))
          end if
        end do

        totalTime=0.
        totime = 0.
        ita1 = 0

        !$omp parallel default(none) private(g) shared(block, start, finish, step)
        ! Set the device (for multi-GPU)
        call set_gpu()
        !$omp do
        do g=start, finish, step
          if (g /= 1) then
           CALL tagging_th(block(g), g)
           CALL cellCount_solid(block(g), g)
          end if
        end do
        !$omp end do
        !$omp end parallel

        ! Just run fine_block_cell on rank 1 (note that we'll need to
        ! be sure that only block(1) is being written to)
        do g=start, finish, step
          if (g == 1) then
            CALL fine_block_cell
            CALL cellCount_solid_coarse(block(1))
          end if
        end do

         ita = 0
         ita1 = 0
         totime = 0.

         do g=start, finish, step
            IF (iStart==0) call initialConditions(block(g), uc)
            IF (iStart==1) call lastConditions(block(g), g, re,totime, ita, ita1)
         end do

        do g=start, finish, step
          if (g /= 1) then
            CALL computeNormDistance(block(g))
            CALL findTScells(block(g))
          end if
          CALL coefficientMatrix(block(g), g == 1)
         end do

        do g=start, finish, step
           CALL non_uni_coeff(block(g))
        end do

        totime = totime + deltat

        do g=start, finish, step
#if USE_HDF5 == 1
        CALL write_output_hdf5(block(g),g)
#else
        CALL write_output_ascii(block(g),g,char_f)
#endif
        end do
        coarse_flcnt_check=0
        DO
        ita = ita + 1
        totime = totime + deltat

        !$omp parallel default(none) private(g) shared(block, start, finish, step, deltat, uc)
        ! Set the device (for multi-GPU)
        call set_gpu()
        !$omp do
        do g=start, finish, step
           CALL nsMomentum2order(block(g))
           if (g == 1) CALL velocityBC(block(g), deltat, uc)
           if (g /= 1) then
             CALL solidCellBC(block(g))
             CALL velocityForcing1(block(g))
           end if
           ! TODO: Not sure why velocityBC is called twice in a row for block(1)?
           if (g == 1) CALL velocityBC(block(g), deltat, uc)
        end do
        !$omp end do
        !$omp end parallel

        CALL poissonSolver(pcItaMax, start, finish, step)

        do g=start, finish, step
           if (g /= 1) CALL pressureForcing1(block(g))
        end do

        do g=start, finish, step
#if USE_HDF5 == 1
           CALL write_output_hdf5(block(g),g)
#else
           CALL write_output_ascii(block(g),g,char_f)
#endif
#ifndef BIOCFD_MPI
           ! TODO: Not yet tested on MPI but should be added!
           CALL writeResult(block(g),g,char_f)
           CALL body_plot(block(g), g)
#endif
       end do

       do g=start, finish, step
         if (g /= 1) then
           !$acc wait
           DEALLOCATE(block(g)%xcent, block(g)%ycent, block(g)%zcent, &
                      block(g)%cosAlpha, block(g)%cosBeta, block(g)%cosGamma, &
                      block(g)%element_length)
           block(g)%blk_mv_tag=0.
           CALL computeSurfaceVariables(block(g),g,phase_angle,piv_pt)
         end if
       END DO

      CALL block_move_check
#ifdef BIOCFD_MPI
      ! Here we set coarse_flcnt_check to 1 on every rank if it is 1
      ! on any rank (by computing the maximum). Since this is really
      ! a true or false flag, we should probably convert this to a
      ! logical. Then we'd use MPI_LOR instead of MPI_MAX for the
      ! reduction operation.
      call MPI_Allreduce(MPI_IN_PLACE, coarse_flcnt_check, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD)
#endif

           do g=start, finish, step
              if (g /= 1) CALL change_block_coords(block(g))
           END DO

           DO g=1,size(intfr)
              CALL change_block_coords_interfaces(intfr(g),block(intfr(g)%b_blk))
           END DO
        do g=1,size(intfr)
           CALL change_block_interface(intfr(g),block(intfr(g)%a_blk),block(intfr(g)%b_blk))
        end do

        do g=start, finish, step
          if (g == 1) then
            CALL fine_block_cell
            CALL cellCount_solid_coarse_mv
          end if
        end do
        ! cellCount_solid_coarse_mv may or may not set
        ! coarse_flcnt_check to 0. If it does then we have to
        ! broadcast it everwhere...
#ifdef BIOCFD_MPI
        call MPI_Bcast(coarse_flcnt_check, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD)
#endif
         do g=start, finish, step
            if (g /= 1) CALL computeSurfaceNorm(block(g))
         end do

#ifdef BIOCFD_MPI
          ! If we are using MPI then at this stage we need to make
          ! sure that block(1) is up-to-date on all ranks. We assume
          ! that all interfaces are from block(1) to another block
          ! This is because at the end of tagging_th_move there are
          ! fineUpdates which will need block(1) to be up to date
          call MPI_Bcast(block(1)%p, size(block(1)%p), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
          call MPI_Bcast(block(1)%u, size(block(1)%u), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
          call MPI_Bcast(block(1)%v, size(block(1)%v), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
          call MPI_Bcast(block(1)%w, size(block(1)%w), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
#endif
         ! Extra care needs to be taken with tagging_th_move because
         ! it contains interpolation functions
         do g=start, finish, step
            if (g /= 1) CALL tagging_th_move(block(g), g)
         end do

        !$omp parallel default(none) private(g) shared(block, start, finish, step)
        ! Set the device (for multi-GPU)
        call set_gpu()
        !$omp do
        DO g=start, finish, step
          if (g == 1) cycle
          CALL selectiveRetagging_th(block(g))
          block(g)%blk_mv_tag=0.
          DEALLOCATE(block(g)%index_ts,block(g)% TSIndexPtr,block(g)% interceptedIndexPtr,&
               block(g)% pNormDis,block(g)% nelp,block(g)% nelu1,block(g)% nelu2,&
               block(g)% nelv1,block(g)% nelv2,block(g)% nelw1,block(g)% nelw2,&
               block(g)% u1NormDis,block(g)% u2NormDis ,block(g)% v1NormDis,&
               block(g)% v2NormDis,block(g)% w1NormDis,block(g)% w2NormDis,&
               block(g)% solidIndexPtr)
          DEALLOCATE(block(g)%fluidIndexPtr,block(g)% redCellIndexPtr, block(g)%blackCellIndexPtr)
          DEALLOCATE(block(g)%p_ghost,block(g)% pt_ghost,block(g)% u2_ghost,block(g)% u2t_ghost,&
               block(g)% v2_ghost,block(g)% v2t_ghost,block(g)% w2_ghost, block(g)%w2t_ghost,&
               block(g)% u1_ghost,block(g)% u1t_ghost,block(g)% v1_ghost,block(g)% v1t_ghost, &
               block(g)%w1_ghost,block(g)% w1t_ghost)
          CALL cellCount_solid(block(g),g)
          call solidCellBC_move(block(g))
          call updateVelocity_newv(block(g))
          block(g)%move_check=0.
          CALL computeNormDistance(block(g))
          CALL findTScells(block(g))
          CALL velocityForcingField(block(g))
          CALL pressureForcingField(block(g))
          CALL velocityForcingGhost(block(g))
          CALL pressureForcingGhost(block(g))
        end do
        !$omp end do
        !$omp end parallel
        IF(ita>=itamax) EXIT
     END DO
     call biocfd_finalize()
      END PROGRAM main


      PROGRAM main
        use, intrinsic :: iso_fortran_env, only: int64, dp => real64, error_unit
        USE global, only: block, blk_start, coarse_flcnt_check, couptime, deltat, istart, &
             ita, ita1, ita2, itamax, nblocks, solvertime, totaltime, totime
        use biocfd_search, only: findDistnode, shiftSurfaceNodesInitial, computeSurfaceNorm, &
             tagging_th, tagging_th_move, block_move_check, cellcount_solid, &
             cellcount_solid_coarse, cellcount_solid_coarse_mv, change_block_coords, &
             change_block_interface, computenormdistance, computesurfacevariables, findtscells, &
             fine_block_cell, selectiveretagging_th
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
        use biocfd_write_output_corner1, only: write_output => write_output_hdf5
#else
        use biocfd_write_output_corner1, only: write_output => write_output_ascii
#endif
        use biocfd_forcing, only: pressureForcing1, pressureforcingfield, pressureforcingghost, &
             velocityforcing1, velocityforcingfield, velocityforcingghost

#ifdef BIOCFD_MPI
        use mpi_f08
#endif
        IMPLICIT NONE

        INTEGER (int64) :: g
        real(dp) :: dstart1, dfinish1

        integer :: rank, num_proc, start_block
#ifdef BIOCFD_MPI
        ! For MPI threading considerations
        integer :: required, provided
        integer :: ierror
#endif

        ! These variables are used to control MPI execution
        rank = 0
        num_proc = 1

        ! These two subroutines are going to be run on every node
        CALL readInput
        CALL readBlockInterface
        CALL interfaceDetail
        ! (until here)

#ifdef BIOCFD_MPI
        required = MPI_THREAD_SERIALIZED
        call MPI_Init_Thread(required, provided, ierror)

        if (provided < required) then
          write(error_unit, *) "MPI does not provide the required threading support. Aborting!"
          call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        end if

        call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

        if (num_proc > size(block)) then
          if (rank == 0) then
            write(error_unit, *) "You have launched the program with more processes than blocks."
            write(error_unit, *) "This will lead to wasted resources."
            write(error_unit, *) "Please resubmit the job with the number of processors <= ", size(block)
          end if
          call MPI_Finalize()
          stop 1
        end if
#endif

        start_block = rank + 1

        do g=start_block, size(block), num_proc
          if (g /= 1) call readSurfaceMeshGmsh(block(g))
        end do
        ! Always allocate the arrays for block(1) - the coarse block,
        ! as will be needed later for interpolation TODO: Work out
        ! whether or not we are really going to need all of these
        ! arrays
        !
        ! We need to allocate block(1) for every rank
        call allocateArrays(block(1))
        do g=start_block, size(block), num_proc
          if (g /= 1) then
            call allocateArrays(block(g))
            call findDistnode(block(g))
            call shiftSurfaceNodesInitial(block(g))
            call computeSurfaceNorm(block(g))
          end if
        end do
        totalTime=0.
        totime = 0.
        ita1 = 0
        ita2 = 0
        solverTime=0.
        coupTime=0.
         do g=start_block, size(block), num_proc
          if (g /= 1) then
            call tagging_th(block(g), g)
            print*,'11'
            call cellCount_solid(block(g), g)
          end if
        end do
        print*,'12'
        ! This should only be called once (for the coarse block)
         do g=start_block, size(block), num_proc
          if (g == 1) then
            CALL fine_block_cell
            print*,'13'
            ! This is only called for the coarse block
            call cellCount_solid_coarse(block(1))
          end if
        end do
        print*,'14'
         do g=start_block, size(block), num_proc
          IF (iStart==0) CALL initialConditions(block(g))
          IF (iStart==1) CALL lastConditions(block(g))
        end do
         do g=start_block, size(block), num_proc
          if (g /= 1) then
            call computeNormDistance(block(g))
            call findTScells(block(g))
          end if
          CALL coefficientMatrix(block(g), g)
          call non_uni_coeff(block(g))
        end do
        totime = totime + deltat
        if ((mod(ita,200_int64) ==0 .or. ita <= 2 )) then
           do g=start_block, size(block), num_proc
               ! TN: TODO: This will not work for HDF5 output!!!
            call write_output(block(g), g)
          end do
        end if
        coarse_flcnt_check=0
        print*, 'adam'

        DO
        ita = ita + 1
        ita2 = ita2 + 1
        totime = totime + deltat

        do g=start_block, size(block), num_proc
          print *, "Rank = ", rank, "g = ", g
          CALL nsMomentum2order(block(g))
          if (g == 1) call velocityBC(block(1))
          if (g /= 1) then
            call solidCellBC(block(g))
            !$acc wait
            call velocityForcing1(block(g))
          end if
        ! I'm not sure we need two calls of velocityBC one after another?
        if (g == 1) call velocityBC(block(1))
      end do

        CALL poissonSolver
        print *,7
        do g=start_block, size(block), num_proc
          if (g /= 1) CALL pressureForcing1(block(g))
        end do
        print *, "Rank = ", rank, " - got to file writing"
        if ((mod(ita,200_int64) ==0 .or. ita <= 2 )) then
          do g=start_block, size(block), num_proc
               ! TN: TODO: This will not work for HDF5 output!!!
            CALL write_output(block(g), g)
          end do
        end if
#ifndef BIOCFD_MPI
        ! TODO: Make sure these run in MPI mode
        !$acc wait
        CALL writeResult
        !$acc wait
        CALL body_plot
#endif
        do g=start_block, size(block), num_proc
          if (g == 1) cycle
           DEALLOCATE(block(g)%xcent, block(g)%ycent, block(g)%zcent,block(g)%cosAlpha, &
                block(g)%cosBeta, block(g)%cosGamma)
            block(g)%blk_mv_tag=0.
        END DO
        print *,10
        do g=start_block, size(block), num_proc
            if (g /= 1) call computeSurfaceVariables(block(g), g)
        end do
        ! TODO: block_move_check iterates over interfaces but we can
        ! probably iterate over blocks instead
        CALL block_move_check
        do g=start_block, size(block), num_proc
          if (g /= 1) call change_block_coords(block(g), g)
          ! TODO: There isn't anything stopping there being an
          ! interface where the block is g==1, but I doubt whether
          ! that is deliberate or not
          call change_block_interface(block(g), g)
          if (g == 1) call fine_block_cell
          if (g == 1) CALL cellCount_solid_coarse_mv
        end do
        ! cellCount_solid_coarse_mv may or may not set
        ! coarse_flcnt_check to 0. If it does then we have to
        ! broadcast it everwhere...
        call MPI_Bcast(coarse_flcnt_check, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD)

          print*,1
          do g=start_block, size(block), num_proc
            if (g /= 1) CALL computeSurfaceNorm(block(g))
          end do
           print*,2

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

        do g=start_block, size(block), num_proc
          if (g /= 1) then
            call tagging_th_move(block(g), g)
            call selectiveRetagging_th(block(g))
            block(g)%blk_mv_tag=0.
            DEALLOCATE(block(g)%index_ts,block(g)% TSIndexPtr,block(g)% interceptedIndexPtr,&
                 block(g)% pNormDis,block(g)% nelp,block(g)% nelu1,block(g)% nelu2,&
                 block(g)% nelv1,block(g)% nelv2,block(g)% nelw1,block(g)% nelw2,&
                 block(g)% u1NormDis,block(g)% u2NormDis ,block(g)% v1NormDis,&
                 block(g)% v2NormDis,block(g)% w1NormDis,block(g)% w2NormDis,&
                 block(g)% solidIndexPtr)
            DEALLOCATE( block(g)%fluidIndexPtr,block(g)% redCellIndexPtr, block(g)%blackCellIndexPtr)
            DEALLOCATE(block(g)%p_ghost,block(g)% pt_ghost,block(g)% u2_ghost,block(g)% u2t_ghost,&
                 block(g)% v2_ghost,block(g)% v2t_ghost,block(g)% w2_ghost, block(g)%w2t_ghost,&
                 block(g)% u1_ghost,block(g)% u1t_ghost,block(g)% v1_ghost,block(g)% v1t_ghost, &
                 block(g)%w1_ghost,block(g)% w1t_ghost)
          end if
       END DO

       do g=start_block, size(block), num_proc
         if (g /= 1) then
           call cellCount_solid(block(g), g)
           call solidCellBC_move(block(g))
           call updateVelocity_newv(block(g))
           block(g)%move_check=0.
           call computeNormDistance(block(g))
           CALL findTScells(block(g))
         end if
       end do

      do g=start_block, size(block), num_proc
        if (g /= 1) then
          call velocityForcingField(block(g))
          call pressureForcingField(block(g))
          call velocityForcingGhost(block(g))
          call pressureForcingGhost(block(g))
        end if
      end do
        IF(ita>=itamax) EXIT
        END DO

#ifdef BIOCFD_MPI
        ! call MPI_Barrier(MPI_COMM_WORLD, ierror)
        call MPI_Finalize(ierror)
#endif
      END PROGRAM main














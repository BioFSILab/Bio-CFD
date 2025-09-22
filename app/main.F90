
      PROGRAM main
        use, intrinsic :: iso_fortran_env, only: int64, dp => real64
        USE global, only: block, blk_start, coarse_flcnt_check, deltat, istart, &
             ita, ita1, ita2, itamax, nblocks, totaltime, totime
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
        IMPLICIT NONE

        INTEGER (int64) :: g
        real(dp) :: dstart1, dfinish1
        CALL readInput
        CALL readBlockInterface
        do g=blk_start, size(block)
          CALL readSurfaceMeshGmsh(block(g))
        end do
        do g=1, size(block)
          CALL allocateArrays(block(g))
        end do
        do g=blk_start, size(block)
           CALL findDistnode(block(g))
           CALL shiftSurfaceNodesInitial(block(g))
        end do
        CALL computeSurfaceNorm
        CALL interfaceDetail
        totalTime=0.
        totime = 0.
        ita1 = 0
        ita2 = 0
        CALL tagging_th
        print*,'11'
        CALL cellCount_solid
        print*,'12'
        CALL fine_block_cell
        print*,'13'
        CALL cellCount_solid_coarse
        print*,'14'
        IF (iStart==0) CALL initialConditions
        IF (iStart==1) CALL lastConditions
        CALL computeNormDistance
        CALL findTScells
        CALL coefficientMatrix
        CALL non_uni_coeff
        totime = totime + deltat
        CALL write_output
        coarse_flcnt_check=0
        print*, 'adam'
        DO
        ita = ita + 1
        ita2 = ita2 + 1
        totime = totime + deltat
        CALL nsMomentum2order
        CALL velocityBC
        CALL solidCellBC
        !$acc wait
        CALL velocityForcing1
        CALL velocityBC
        CALL poissonSolver
        print *,7
        CALL pressureForcing1
        CALL write_output
        !$acc wait
        CALL writeResult
        !$acc wait
        CALL body_plot
        DO g=blk_start, nblocks
           DEALLOCATE(block(g)%xcent, block(g)%ycent, block(g)%zcent,block(g)%cosAlpha, &
                block(g)%cosBeta, block(g)%cosGamma)
            block(g)%blk_mv_tag=0.
        END DO
        print *,10
            CALL computeSurfaceVariables
           CALL block_move_check
           CALL change_block_coords
           CALL change_block_interface
          CALL fine_block_cell
          CALL cellCount_solid_coarse_mv
           print*,1
           CALL computeSurfaceNorm
           print*,2
           CALL tagging_th_move
           CALL selectiveRetagging_th
        DO g=blk_start, nblocks
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
       END DO

            CALL cellCount_solid
        DO g=blk_start, nblocks
            call solidCellBC_move(g)
             call updateVelocity_newv(g)
            block(g)%move_check=0.
        ENDDO
           CALL computeNormDistance
        CALL findTScells
        CALL cpu_time(dStart1)
        CALL velocityForcingField
        CALL pressureForcingField
        CALL cpu_time(dFinish1)
        CALL cpu_time(dStart1)
        CALL velocityForcingGhost
        CALL pressureForcingGhost
        CALL cpu_time(dFinish1)
        IF(ita>=itamax) EXIT
        END DO
      END PROGRAM main














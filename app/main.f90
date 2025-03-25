
      PROGRAM main
        USE global
        IMPLICIT NONE

        REAL (KIND=8) :: ts, te
        INTEGER (KIND=8) :: g, i
        CALL readInput
        CALL readBlockInterface
        CALL readSurfaceMeshGmsh
        !CALL readSurfaceMeshGambit
        CALL allocateArrays
        CALL findDistnode
        CALL shiftSurfaceNodesInitial
        CALL computeSurfaceNorm
        CALL interfaceDetail
        totalTime=0.
        totime = 0.
        ita1 = 0
        ita2 = 0
        solverTime=0.
        coupTime=0.
        !CALL body_plot
        !CALL writeOutput1
        CALL tagging_th
        !CALL readTagging
        !CALL writeOutput1
        !CALL writeTagging
        print*,'11'
        CALL cellCount_solid
        print*,'12'
        CALL fine_block_cell
        print*,'13'
        CALL cellCount_solid_coarse
        print*,'14'
        !CALL cellCount
        IF (iStart==0) CALL initialConditions
        IF (iStart==1) CALL lastConditions
        CALL computeNormDistance
        CALL findTScells
        !ita = 1
        !ita1 = 1
        CALL coefficientMatrix
        !CALL amgx_matrix
        CALL non_uni_coeff
        !ita = 1
        !totime = 0.
        !CALL readComputeSumData
        !CALL readstressData
        !ita = ita + 1
        !ita1 = ita1 + 1
        totime = totime + deltat
        !CALL calculateUcMO
        !CALL calculateUc
!       CALL nsMomentum
!       print*,'1'
!       CALL velocityBC
!       print*,'2'
!       CALL solidCellBC
!       print*,'3'
!       CALL velocityForcing1
!       !CALL velocityForcing1_gh
!       print*,'after vforc'
!       CALL velocityBC
!       print*,'after vbc'
!       CALL poissonSolver
!       print*,'after psolv'
!       CALL pressureForcing1
!       !CALL pressureForcing1_gh
!       print*,'after pforc'
!       !CALL calculateUcMO
!       !CALL calculateUc
!       !CALL computeSumData
        !CALL writeOutput1
        st_flag=0
!!        CALL stressCal1
        CALL writeOutput1
        coarse_flcnt_check=0
        print*, 'adam'
        DO
                ita = ita + 1
                ita2 = ita2 + 1
                totime = totime + deltat
                !CALL nsMomentumAB
                !!$acc wait
                !print *,1
                !CALL nsMomentum
                CALL nsMomentum2order
        ! if (ita .eq. 258)then
        ! CALL writeOutput1
        ! print*,'after ns check'
        ! pause
        ! endif
                !!$acc wait
        ! print *,2
                CALL velocityBC
                !!$acc wait
                !print *,3
                CALL solidCellBC
                !$acc wait
                !print *,4
                CALL velocityForcing1
                !CALL velocityForcing1_gh
                !!$acc wait
                !print *,5
                CALL velocityBC
        ! !!$acc wait
                !print *,6
                CALL poissonSolver
        ! if (ita .eq. 258)then
        ! CALL writeOutput1
        ! print*,'after psolv check'
        ! pause
        ! endif
        ! !!$acc wait
                print *,7
                CALL pressureForcing1
                !CALL pressureForcing1_gh
                !!$acc wait
                !CALL calculateUcMO
                !CALL calculateUc
                !!$acc update host (u, v, w, ut, vt, wt, p, cell)
        ! CALL writeOutput
                CALL writeOutput1
        !       if (ita .gt.60000 )then
        !       ita1 = ita1 + 1
        !       CALL computeSumData
        !       CALL computeSumSqData
        !       CALL writeComputeSumData
        !       CALL writeComputeSumSqData
        !       !print*, ita1
        !       end if
        !       IF(ita1 .gt.0 .and. mod(ita1,10000) .eq. 0) THEN
        !           CALL computeAvgData
        !            CALL computePhaseAvgData

        !       !print*, ita1
        !       !if (st_flag .eq. 0 )then
        !       !CALL stressCal1_fl
        !       !end if
        !       !CALL stressCal2_fl
        !
        !            CALL writeAvgoutput
        !            CALL writeAvgoutput_fl
        !       END IF
                !print*,'2'
                !print *,8
                !!$acc wait
        !!        CALL stressCal2
                !print*,'3'

                !print *,9
                !CALL writeStressData_fl
        !       CALL writeStressData
        !       CALL write_fl_points
                !$acc wait
                !print *,10
                CALL writeResult
                !$acc wait
                CALL body_plot
        ! print*,solverTime, coupTime
                DO g=blk_start, nblocks
                DEALLOCATE(block(g)%xcent, block(g)%ycent, block(g)%zcent,block(g)%cosAlpha, block(g)%cosBeta, block(g)%cosGamma,block(g)%alpha3,block(g)%beta3, block(g)%gamma3)
                block(g)%blk_mv_tag=0.
                END DO
                !DEALLOCATE (cell2)
                print *,10
                CALL computeSurfaceVariables
                CALL block_move_check
                !print*,ita,block(2)%y1(1),block(2)%yp(1)
                CALL change_block_coords
                !print*,ita,block(2)%y1(1),block(2)%yp(1)
                CALL change_block_interface
                CALL fine_block_cell
                CALL cellCount_solid_coarse_mv
                !print*,ita,block(2)%y1(1),block(2)%yp(1)
                print*,1
                CALL computeSurfaceNorm
                print*,2
                CALL tagging_th_move
                ! print*,3
        !  if (block(2)%blk_mv_tag .eq. 1)then
        ! CALL writeOutput1
        ! print*,'!!!!!!!!check!!!!!!!!!!!!'
        ! pause
        ! endif
                CALL selectiveRetagging_th
                ! print*,4
                !!print*,'!!!!yp!!!!!!!!!',block(2)%yp(160)
        !print*,'After selective retagging'
        !  CALL amgx_matrix
        !print*,'After amgx_matrix'
                !CALL tagging
                DO g=blk_start, nblocks
                block(g)%blk_mv_tag=0.
                        DEALLOCATE(block(g)%index_ts,block(g)% TSIndexPtr,block(g)% interceptedIndexPtr,block(g)% pNormDis,block(g)% nelp,block(g)% nelu1,block(g)% nelu2,block(g)% nelv1,block(g)% nelv2,block(g)% nelw1,block(g)% nelw2,block(g)% u1NormDis,block(g)% u2NormDis ,block(g)% v1NormDis,block(g)% v2NormDis,block(g)% w1NormDis,block(g)% w2NormDis,block(g)% solidIndexPtr)
                DEALLOCATE( block(g)%fluidIndexPtr,block(g)% redCellIndexPtr, block(g)%blackCellIndexPtr)
                DEALLOCATE(block(g)%p_ghost,block(g)% pt_ghost,block(g)% u2_ghost,block(g)% u2t_ghost,block(g)% v2_ghost,block(g)% v2t_ghost,block(g)% w2_ghost, block(g)%w2t_ghost,block(g)% u1_ghost,block(g)% u1t_ghost,block(g)% v1_ghost,block(g)% v1t_ghost, block(g)%w1_ghost,block(g)% w1t_ghost)
        !   	DEALLOCATE(block(g)%pNormDis, block(g)%nelp, &
        !       block(g)%nelu1, block(g)%nelu2, block(g)%nelv1, &
        !       block(g)%nelv2, block(g)%nelw1, block(g)%nelw2, &
        !       block(g)%u1NormDis, block(g)%u2NormDis , block(g)%v1NormDis, &
        !       block(g)%v2NormDis , block(g)%w1NormDis, block(g)%w2NormDis)
                !  DEALLOCATE(block(g)%interceptedIndexPtr, block(g)%pNormDis,block(g)%nelp, &
                !    block(g)%nelu1, block(g)%nelu2, block(g)%nelv1, block(g)%nelv2, block(g)%u1NormDis, block(g)%u2NormDis, &
                !    block(g)%v1NormDis, block(g)%v2NormDis,block(g)%solidIndexPtr)
        !print*,'After deallocate'
                !DO i=1,nblocks
        !                DEALLOCATE(block(g)%interceptedIndexPtr,block(g)%fluidIndexPtr,block(g)%redCellIndexPtr,block(g)%blackCellIndexPtr,block(g)%solidIndexPtr)
                !END DO
        END DO
        !print*,'After deallocate'

                CALL cellCount_solid
                DO g=blk_start, nblocks
                call solidCellBC_move(g)
                call updateVelocity_newv(g)
                !call writeOutput
                !pause
        !            if (block(g)%move_check .eq. 1)then
        !       CALL writeOutput1
        !       print*,'!!!!!!!!check2!!!!!!!!!!!!'
        !               pause
        !               endif
                block(g)%move_check=0.
                ENDDO
                !CALL coarseUpdate
                !CALL cellCount
                CALL computeNormDistance
                CALL findTScells
        ! 2 CONTINUE
        !!  print*,'11'
                CALL cpu_time(dStart1)
                CALL velocityForcingField
        !!   print*,'12'
                CALL pressureForcingField
                CALL cpu_time(dFinish1)
        !!   print*,'13'
        !        PRINT*, 'time Field Forcing = ', dStart1 - dFinish1
                CALL cpu_time(dStart1)
                CALL velocityForcingGhost
        !!   print*,'14'
                CALL pressureForcingGhost
        !!   print*,'15'
                CALL cpu_time(dFinish1)
                IF(ita>=itamax) EXIT
        END DO
      END PROGRAM main














MODULE biocfd_blocks
       use, intrinsic :: iso_fortran_env, only: dp => real64, int32, int64
       IMPLICIT NONE
       private

       public :: Blocks

        type Blocks
           REAL(dp) ::  dx,dy, dz, ypth2, xpth1, xpth2, xchg,ychg, yt, ydot, yddot,&
                bfreq, yamp, xt, xdot
        REAL(dp) :: derr1,derr2,derrStdSt,a0
        REAL(dp) :: xshift, yshift, zshift, gx_shift, gy_shift,gz_shift
        INTEGER (int64):: nx, ny, nz
        INTEGER (int64):: mk, mkx1, mkx2
        INTEGER (int64):: nIterPcor, cintp, fineg
        INTEGER (int64) ::  k_startSearch, k_endSearch, &
                                     j_startSearch, j_endSearch, &
                                     i_startSearch, i_endSearch
       INTEGER (int64), ALLOCATABLE, DIMENSION (:,:,:) :: cell, cell2, cell_n, cell_pr,nodeIdTag

       REAL (dp), ALLOCATABLE, DIMENSION (:) :: deltax, deltay, deltaz, x1, y1, z1, &
            xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp

       real(dp), allocatable, dimension(:, :) :: ca_uu, ck_uu, ca_vv, ck_vv, ca_ww, ck_ww, &
                                                 ca_uv, ck_uv, ca_uw, ck_uw, ca_vu, ck_vu, &
                                                 ca_vw, ck_vw, ca_wu, ck_wu, ca_wv, ck_wv


       REAL (dp), ALLOCATABLE, DIMENSION (:, :)    :: Acx, Acy, Acz
       REAL (dp), ALLOCATABLE, DIMENSION (:, :, :) :: b, u, u_dum, ut, &
                                                      v, vt, v_dum, &
                                                      w, wt, w_dum, &
                                                      p, p_dum, pc, pco

        INTEGER (int64), ALLOCATABLE, DIMENSION (:, :) :: fluidIndexPtr, redCellIndexPtr, &
                                                          blackCellIndexPtr, nodeId
        INTEGER(int64) :: ibCellCount, solidCellCount, fluidCellCount, redCellCount, &
                          blackCellCount, TSCellCount

       INTEGER (int64), ALLOCATABLE, DIMENSION (:, :) :: TSIndexPtr, interceptedIndexPtr, &
                                                         solidIndexPtr


       INTEGER (int64), ALLOCATABLE, DIMENSION (:) :: nelp, nelu1, nelu2, nelv1, nelv2, nelw1, nelw2
       REAL (dp), ALLOCATABLE, DIMENSION (:) :: xcent, ycent, zcent, &
                                                cosAlpha, cosBeta, cosGamma, &
                                                pNormDis, u1NormDis, u2NormDis , &
                                                v1NormDis, v2NormDis, w1NormDis, w2NormDis, &
                                                p_ghost, pt_ghost, u2_ghost, u2t_ghost, &
                                                v2_ghost, v2t_ghost, w2_ghost, &
                                                w2t_ghost, u1_ghost, u1t_ghost, &
                                                v1_ghost, v1t_ghost, w1_ghost, w1t_ghost

       INTEGER (int64), ALLOCATABLE, DIMENSION (:) ::  ibSurfId,ibElP1, ibElP2, ibElP3
       INTEGER (int64), ALLOCATABLE, DIMENSION (:) :: ibNodeId,index_ts
       REAL (dp), ALLOCATABLE, DIMENSION (:) :: xnode, ynode,  znode, xnode1, ynode1, znode1
       INTEGER (int64) :: ibElems, ibNodes
       INTEGER (int64) :: move_check, move_amty, move_amtx,move_amtz,blk_mv_tag
       REAL(dp) :: ymove,ypos
       !zpos not used
       REAL (dp) :: xmove,xpos, zmove,zpos
       REAL (dp) :: inity_cent, initx_cent, nxty_cent,nxtx_cent
       REAL (dp) :: initz_cent,nxtz_cent
       INTEGER (int64) :: cpy_x_start_mv, cpy_x_end_mv, cpy_y_start_mv, cpy_y_end_mv
       INTEGER (int64) :: cpy_z_start_mv, cpy_z_end_mv
       INTEGER (int64) :: cpy_x_start, cpy_x_end, cpy_y_start, cpy_y_end
       INTEGER (int64) :: cpy_z_start, cpy_z_end
       REAL (dp) :: theta, thetaDot, thetaDDot, piv_x,piv_y, piv_z
       REAL (dp) :: alphaDot, alphaDDot, thetaDot1, thetaDDot1, thetaDot2, thetaDDot2


    end type Blocks

END MODULE biocfd_blocks


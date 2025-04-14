! allow(missing-accessibility-statement) - TODO: Attempt to remove in future
MODULE global
       use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int32, int64
       use iso_c_binding,only :c_int,c_double,c_loc,c_ptr
       IMPLICIT NONE
       CHARACTER (LEN = 128) :: line
       CHARACTER (LEN = 3)   :: char_f
       INTEGER               :: istart, id1, st_flag
       INTEGER (int64)   :: itamax, pcItaMax,amgxita,      &
                                ita, nIterPcor, nc, ita1, ital,ita2,    &
                                sumIterPc,itaSola, totIterPc, inor,blk_start, coarse_flcnt_check
       REAL (dp)        :: lx, ly, dt_order,  &
                                omega,omega1,omega2,omega3,omega4,  &
                                deltx2, delty2,  deltz2, dxmin, &
                                freq, &
                                u0, v0, w0, p0,Uavg,  &
                                eps1, epsi, epsDiv, re, rev, divmax,   &
                                fx, fy,  alpha, pct, pct1, dstart1,dfinish1, &
                                xfact,deltat, coupTime,totime, totalTime, dfinish, dstart, &
                                solverTime, pi , msTime, mindx, al, uc

       REAL (dp)       :: p_new1, p_new2, p_final, u_new, v_new, w_new
       REAL (dp)       :: alpha_m, theta_m, alpha_m1, theta_m1, mu_f, rho_f, l_c, u_tip, lwing, disp
       !!!variables for Orlanski multiple outlet
       REAL (dp)       :: y11, y12, z11, z12, y21, y22, z21, z22, uc11, uc22


        INTEGER (int64) ::nblocks, intflines

        type Blocks
           REAL(dp) ::  dx,dy, dz,ypth1, ypth2, xpth1, xpth2, xchg,ychg, yt, ydot, yddot,&
                bfreq, yamp, xt, xdot, xddot
        REAL(dp) :: derr1,derr2,derrStdSt,a0
        REAL(dp) :: xshift, yshift, zshift, gx_shift, gy_shift,gz_shift
        INTEGER (int64):: nx, ny, nz
        INTEGER (int64):: mk, mkx1, mkx2
        INTEGER (int64):: nIterPcor, cintp, fineg
        INTEGER (int64) ::  k_startSearch, k_endSearch, &
                                     j_startSearch, j_endSearch, &
                                     i_startSearch, i_endSearch
       INTEGER (int64), ALLOCATABLE, DIMENSION (:,:,:) :: cell, cell2, cell_n, cell_pr,nodeIdTag
       INTEGER (int64), ALLOCATABLE, DIMENSION (:,:,:) :: minElemcell

       REAL (dp), ALLOCATABLE, DIMENSION (:) :: deltax, deltay, deltaz, x1, y1, z1, &
            xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp, &
            ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
            ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
            ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
            ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
            ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
            ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww, &
            ca1_uv, ca2_uv, ca3_uv, ca4_uv, ca5_uv, ca6_uv, &
            ck1_uv, ck2_uv, ck3_uv, ck4_uv, ck5_uv, ck6_uv, &
            ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
            ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
            ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
            ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
            ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
            ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
            ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
            ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
            ca1_wv, ca2_wv, ca3_wv, ca4_wv, ca5_wv, ca6_wv, &
            ck1_wv, ck2_wv, ck3_wv, ck4_wv, ck5_wv, ck6_wv


       REAL (dp), ALLOCATABLE, DIMENSION (:, :)    :: A, An,Ac, Acx, Acy, Acz
       REAL (dp), ALLOCATABLE, DIMENSION (:, :, :) :: b, u, u_dum, ut, u_sum, u_avg,   &
                                                      v, vt, v_dum, v_sum, v_avg,      &
                                                      w, wt, w_dum, w_sum, w_avg,      &
                                                      p, p_sum,p_dum, p_avg, pc, pco, &
                                                      uv_sum,vw_sum,uw_sum, &
                                                      uv_avg,vw_avg,uw_avg, &
                                                      uflu_avg,vflu_avg,wflu_avg,pflu_avg, &
                                                      uflu_rms,vflu_rms,wflu_rms,pflu_rms, &
                                                      uvflu_avg,vwflu_avg,uwflu_avg, &
                                                      u2_sum,v2_sum,w2_sum,p2_sum,  &
                                                      u2_avg,v2_avg,w2_avg,p2_avg,  &
                                                      ufl,vfl,wfl,  &
                                                      resi_u, resi_v, resi_w

        REAL(sp) , ALLOCATABLE, DIMENSION (:, :, :) :: xp1, yp1, zp1
        REAL(sp) , ALLOCATABLE, DIMENSION (:, :, :) :: xpn1, ypn1, zpn1

        INTEGER (int64), ALLOCATABLE, DIMENSION (:, :) :: fluidIndexPtr, redCellIndexPtr, &
                                                          blackCellIndexPtr, nodeId
        INTEGER(int64) :: ibCellCount, solidCellCount, fluidCellCount, redCellCount, &
                          blackCellCount, TSCellCount
        REAL (dp), ALLOCATABLE, DIMENSION (:) ::dataval
        INTEGER (int32)   :: nnz,nu, nit, nit1
        INTEGER (int32), ALLOCATABLE, DIMENSION (:) :: row_ptr,col
        REAL (dp), ALLOCATABLE, DIMENSION (:) ::rhs,sol
        INTEGER (int32) :: crs_data(4)
         integer(int32) :: diag(7)

       INTEGER (int64), ALLOCATABLE, DIMENSION (:, :) :: TSIndexPtr, interceptedIndexPtr, &
                                                         solidIndexPtr, fluidInterceptedIndexPtr


       INTEGER (int64), ALLOCATABLE, DIMENSION (:) :: nelp, nelu1, nelu2, nelv1, nelv2, nelw1, nelw2
       REAL (dp), ALLOCATABLE, DIMENSION (:) :: xcent, ycent, zcent, &
                                                cosAlpha, cosBeta, cosGamma, alpha3, beta3, &
                                                gamma3, pNormDis, u1NormDis, u2NormDis , &
                                                v1NormDis, v2NormDis, w1NormDis, w2NormDis, &
                                                p_ghost, pt_ghost, u2_ghost, u2t_ghost, &
                                                v2_ghost, v2t_ghost, w2_ghost, &
                                                w2t_ghost, u1_ghost, u1t_ghost, &
                                                v1_ghost, v1t_ghost, w1_ghost, w1t_ghost
       INTEGER (int64):: ibElemCnt, ibElemCntSt
       INTEGER (int64), ALLOCATABLE, DIMENSION (:,:) :: Elemcell, ucell, vcell, wcell, pcell
       REAL (dp), ALLOCATABLE, DIMENSION (:) :: areaElem
       REAL (dp), ALLOCATABLE, DIMENSION (:,:) :: stressElem
       REAL(dp) :: modStressNode, modSIGNWSS

       REAL (dp) :: u_init, u_final, v_init, v_final, w_init, w_final, &
                    Total_Force_X, Total_Force_Y, Total_Force_Z

       INTEGER (int64), ALLOCATABLE, DIMENSION (:) :: ibNodeId, ibSurfId,ibElP1, ibElP2, ibElP3, &
                                                      index_ts
       REAL (dp), ALLOCATABLE, DIMENSION (:) :: xnode, ynode,  znode, xnode1, ynode1, znode1, bcSurf
       INTEGER (int64) :: ibElems, ibNodes
       INTEGER (int64) :: tp_pt,bt_pt,lt_pt,rt_pt
       INTEGER (int64) :: move_check, move_amty, move_amtx,move_amtz,blk_mv_tag
       REAL (dp) :: xshift_move,yshift_move,zshift_move
       REAL (dp) :: u_prev,u_curr,v_prev,v_curr,w_prev,w_curr,Total_VP_FY,Total_VP_FX
       REAL (dp) :: rhof, rhop, dp, accn_g, volp, massp,Total_FY,accnp_Y,ymove,ypos
       REAL (dp) :: Total_FX,accnp_X,xmove,xpos, zmove,zpos
       REAL (dp) :: Total_V_Fy, Total_P_Fy, Total_V_Fx,Total_P_Fx
       REAL (dp) :: inity_cent, initx_cent, nxty_cent,nxtx_cent
       REAL (dp) :: initz_cent,nxtz_cent
       INTEGER (int64) :: cpy_x_start_mv, cpy_x_end_mv, cpy_y_start_mv, cpy_y_end_mv
       INTEGER (int64) :: cpy_z_start_mv, cpy_z_end_mv
       INTEGER (int64) :: cpy_x_start, cpy_x_end, cpy_y_start, cpy_y_end
       INTEGER (int64) :: cpy_z_start, cpy_z_end
       REAL (dp) :: theta, thetaDot, thetaDDot, piv_x,piv_y, piv_z, theta_i, thetaDot_i, &
                    theta2, theta2Dot, theta2DDot, piv2_x, piv2_y, piv2_z,theta2_i, theta2Dot_i, &
                    alphaDot, alphaDDot,ac_x_al, ac_y_al, at_x_al, at_y_al, &
                    aoa, thetaDot1, thetaDDot1, thetaDot2, thetaDDot2
       REAL (dp), ALLOCATABLE, DIMENSION (:) :: xp_dum, yp_dum, zp_dum


       end type Blocks
       REAL(dp), ALLOCATABLE, DIMENSION (:) :: INSTWSS, SUMWSS, SQSUMWSS
       REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: SIGNWSS
       REAL(dp), ALLOCATABLE, DIMENSION(:):: TAWSS, OSI, RRT,WSSRMS
       REAL(dp), ALLOCATABLE, DIMENSION(:,:):: Afnode, stressNode, Anode
       REAL(dp) :: ac_x_al, ac_y_al, at_x_al, at_y_al, ac_x, ac_y, ac_z, at_x, at_y ,at_z

        type Interfaces

        INTEGER(int64), ALLOCATABLE, DIMENSION (:,:) :: px_interface_det, ux_interface_det, &
                                                        vx_interface_det, wx_interface_det, &
                                                        pz_interface_det, uz_interface_det, &
                                                        vz_interface_det, wz_interface_det, &
                                                        py_interface_det, uy_interface_det, &
                                                        vy_interface_det, wy_interface_det
        REAL(dp) :: a_blk, b_blk, xintf_start, xintf_end, yintf_start, yintf_end, &
                    zintf_start, zintf_end
        REAL(dp) :: xintf_st_new, xintf_en_new, yintf_st_new, yintf_en_new
        REAL(dp) :: zintf_st_new, zintf_en_new
        INTEGER(int64) :: a_msh, b_msh, a_intf, b_intf
        INTEGER(int64) :: counterxu,counteryu,counterxp,counteryp,counterxv,counteryv
        INTEGER(int64) :: counterxw,counteryw
        INTEGER(int64) :: counterzu,counterzp,counterzv, counterzw
        INTEGER(int64) :: cpx_start, cpy_start,cpx_end, cpy_end
        INTEGER(int64) :: cux_start, cuy_start,cux_end, cuy_end
        INTEGER(int64) :: cvx_start, cvy_start,cvx_end, cvy_end
        INTEGER(int64) :: ccellx_start, ccelly_start,ccellx_end, ccelly_end
        INTEGER(int64) :: fvx_start, fvy_start
        INTEGER(int64) :: fpx_start, fpy_start
        INTEGER(int64) :: fux_start, fuy_start
        INTEGER(int64) :: fcellx_start, fcelly_start
        end type Interfaces

        type(Blocks),allocatable ::block(:)
        type(Interfaces),allocatable ::intfr(:)


      REAL (dp) :: d_fine_x1, d_fine_y1, d_coarse_x1, d_coarse_y1, d_val_x1, d_val_y1, &
                   d_fine_frac, d_coarse_frac
      REAL (dp) :: l_intp_valx, l_intp_x1,l_intp_x2,l_intp_y1, l_intp_y2, l_intp_valy


       INTEGER(int64)  :: nfl_blk
       INTEGER(int64),ALLOCATABLE,DIMENSION(:) :: fl_blk

       ! ibm variables
       INTEGER  :: bcType
       REAL (dp)    :: aoa, piv_pt, var_surf, aoa1
       REAL (dp)    :: phase_angle, aoa2, a0y, ang_theta, alpha_t, theta_t
       INTEGER (int64)   :: surGeoPoints

        !amgx
        INTEGER (int32) :: nnz_max,nu_max
        type(c_ptr)::cptr_crs,cptr_data,cptr_col,cptr_row,cptr_rhs,cptr_sol,cptr_nit, cptr_nit1
        INTEGER(int32) ::init_stat,dest_stat,solve_stat,amgx_checker
END MODULE global

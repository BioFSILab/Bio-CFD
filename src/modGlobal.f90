       MODULE global
       use iso_c_binding,only :c_int,c_double,c_loc,c_ptr
       IMPLICIT NONE
       CHARACTER (LEN = 128) :: line
       CHARACTER (LEN = 3)   :: char_f
       INTEGER               :: istart, id1, st_flag
       INTEGER ( KIND = 8)   :: itamax, pcItaMax,amgxita,      &
                                ita, nIterPcor, nc, ita1, ital,ita2,    &
                                sumIterPc,itaSola, totIterPc, inor,blk_start, coarse_flcnt_check 
       REAL (KIND =8)        :: lx, ly, dt_order,  &
                                omega,omega1,omega2,omega3,omega4,  & 
                                deltx2, delty2,  deltz2, dxmin, &
                                freq, &
                                u0, v0, w0, p0,Uavg,  &  
                                eps1, epsi, epsDiv, re, rev, divmax,   &
                                fx, fy,  alpha, pct, pct1, dstart1,dfinish1, &
                                xfact,deltat, coupTime,totime, totalTime, dfinish, dstart, solverTime, pi , msTime, mindx, al, uc
       
       REAL (KIND = 8)       :: p_new1, p_new2, p_final, u_new, v_new, w_new 
       REAL (KIND = 8)       :: alpha_m, theta_m, alpha_m1, theta_m1, mu_f, rho_f, l_c, u_tip, lwing, disp 
       !!!variables for Orlanski multiple outlet
       REAL (KIND = 8)       :: y11, y12, z11, z12, y21, y22, z21, z22, uc11, uc22

                         
       !REAL (KIND =8)        :: u_init, u_final, v_init, v_final, w_init, w_final, xmove, ymove, &
       !                         zmove, Total_Force_X, Total_Force_Y, Total_Force_Z                               


        INTEGER (KIND=8) ::nblocks, intflines

        type Blocks
        REAL(KIND=8) :: xstart, xend, ystart, yend, zstart, zend, dx,dy, dz, ypth1, ypth2, xpth1, xpth2, xchg,ychg, yt, ydot, yddot, bfreq, yamp, xt, xdot, xddot
        REAL(KIND=8) :: derr1,derr2,derrStdSt,a0
        REAL(KIND=8) :: xshift, yshift, zshift, gx_shift, gy_shift,gz_shift
        INTEGER (KIND=8):: nx, ny, nz
        INTEGER (KIND=8):: mk, mkx1, mkx2
        INTEGER (KIND=8):: nIterPcor, cintp, fineg
        INTEGER (KIND=8) ::  k_startSearch, k_endSearch, &
                                     j_startSearch, j_endSearch, &
                                     i_startSearch, i_endSearch
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:,:,:)        :: cell, cell2, cell_n, cell_pr,nodeIdTag         
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:,:,:)    :: minElemcell

        REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)       :: deltax, deltay, deltaz, x1, y1, z1,             &
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
 

       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:, :)    :: A, An,Ac, Acx, Acy, Acz
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:, :, :) :: b, u, u_dum, ut, u_sum, u_avg,   &
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

         REAL*4 , ALLOCATABLE, DIMENSION (:, :, :) :: xp1, yp1, zp1
         REAL*4 , ALLOCATABLE, DIMENSION (:, :, :) :: xpn1, ypn1, zpn1

       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:, :) :: fluidIndexPtr,     &
                                                            redCellIndexPtr, blackCellIndexPtr, nodeId
         INTEGER(KIND=8)  ::  ibCellCount, solidCellCount, fluidCellCount, redCellCount, blackCellCount, TSCellCount
         !REAL (KIND = 8), ALLOCATABLE,TARGET, DIMENSION (:) ::dataval
        REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) ::dataval
        INTEGER (KIND=4)   :: nnz,nu, nit, nit1
        !INTEGER (KIND=4), ALLOCATABLE,TARGET, DIMENSION (:) ::row_ptr,col
        INTEGER (KIND=4), ALLOCATABLE, DIMENSION (:) :: row_ptr,col
        REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) ::rhs,sol
        !REAL (KIND = 8), ALLOCATABLE,TARGET, DIMENSION (:) ::rhs,sol
        INTEGER (KIND=4) :: crs_data(4)
         integer(KIND=4) :: diag(7)
       
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:, :) :: TSIndexPtr, interceptedIndexPtr, solidIndexPtr,     &
                                                            fluidInterceptedIndexPtr


        INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:)   :: nelp, nelu1, nelu2, nelv1, nelv2, nelw1, nelw2
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) ::  xcent, ycent, zcent, &
                                                      cosAlpha, cosBeta, cosGamma, alpha3, beta3, gamma3, &
                                                      pNormDis, u1NormDis, u2NormDis , v1NormDis, v2NormDis, w1NormDis, w2NormDis, &
                                                      p_ghost, pt_ghost, u2_ghost, u2t_ghost, v2_ghost, v2t_ghost, w2_ghost, &
                                                      w2t_ghost, u1_ghost, u1t_ghost, v1_ghost, v1t_ghost, w1_ghost, w1t_ghost
        INTEGER (KIND=8):: ibElemCnt, ibElemCntSt
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:,:)      :: Elemcell, ucell, vcell, wcell, pcell
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)           ::  areaElem                                             
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:,:)           ::   stressElem                                              
       REAL(KIND = 8):: modStressNode, modSIGNWSS  
  
       REAL (KIND =8)        :: u_init, u_final, v_init, v_final, w_init, w_final,  &
                                Total_Force_X, Total_Force_Y, Total_Force_Z                               
 
       INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:)   :: ibNodeId, ibSurfId,ibElP1, ibElP2, ibElP3, index_ts
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: xnode, ynode,  znode, xnode1, ynode1, znode1, bcSurf
       INTEGER (KIND=8)   ::  ibElems, ibNodes
       INTEGER (KIND = 8) :: tp_pt,bt_pt,lt_pt,rt_pt
       INTEGER (KIND = 8) :: move_check, move_amty, move_amtx,move_amtz,blk_mv_tag
       REAL (KIND = 8) :: xshift_move,yshift_move,zshift_move
       REAL (KIND=8) :: u_prev,u_curr,v_prev,v_curr,w_prev,w_curr,Total_VP_FY,Total_VP_FX
       REAL (KIND=8) :: rhof, rhop, dp, accn_g, volp, massp,Total_FY,accnp_Y,ymove,ypos
       REAL (KIND=8) :: Total_FX,accnp_X,xmove,xpos, zmove,zpos
       REAL (KIND=8) :: Total_V_Fy, Total_P_Fy, Total_V_Fx,Total_P_Fx
       REAL (KIND=8) :: inity_cent, initx_cent, nxty_cent,nxtx_cent
       REAL (KIND=8) :: initz_cent,nxtz_cent
       INTEGER (KIND=8) :: cpy_x_start_mv, cpy_x_end_mv, cpy_y_start_mv, cpy_y_end_mv
       INTEGER (KIND=8) :: cpy_z_start_mv, cpy_z_end_mv
       INTEGER (KIND=8) :: cpy_x_start, cpy_x_end, cpy_y_start, cpy_y_end
       INTEGER (KIND=8) :: cpy_z_start, cpy_z_end
       REAL (KIND = 8)    :: theta, thetaDot, thetaDDot, piv_x,piv_y, piv_z, theta_i, thetaDot_i, &
        theta2, theta2Dot, theta2DDot, piv2_x, piv2_y, piv2_z,theta2_i, theta2Dot_i, alphaDot, alphaDDot,ac_x_al, ac_y_al, at_x_al, at_y_al, & 
        aoa, thetaDot1, thetaDDot1, thetaDot2, thetaDDot2
       REAL (KIND=8), ALLOCATABLE, DIMENSION (:)       :: xp_dum, yp_dum, zp_dum

 
        end type Blocks
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)           :: INSTWSS,SUMWSS,  SQSUMWSS                                             
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:,:)           ::  SIGNWSS                                              
       REAL(KIND = 8), ALLOCATABLE,DIMENSION(:):: TAWSS, OSI, RRT,WSSRMS       
       REAL(KIND = 8),ALLOCATABLE, DIMENSION(:,:):: Afnode, stressNode, Anode             
       REAL (KIND = 8)    :: ac_x_al, ac_y_al, at_x_al, at_y_al, ac_x, ac_y, ac_z, at_x, at_y ,at_z

        type Interfaces

       !INTEGER (KIND=8),DIMENSION(2) :: px_coarse_intf, px_fine_intf,ux_coarse_intf, ux_fine_intf, vx_coarse_intf,vx_fine_intf,wx_coarse_intf, wx_fine_intf
       !INTEGER (KIND=8),DIMENSION(2) :: pz_coarse_intf, pz_fine_intf,uz_coarse_intf, uz_fine_intf, vz_coarse_intf,vz_fine_intf,wz_coarse_intf, wz_fine_intf
       !INTEGER (KIND=8),DIMENSION(2) :: py_coarse_intf, py_fine_intf,uy_coarse_intf, uy_fine_intf, vy_coarse_intf,vy_fine_intf,wy_coarse_intf, wy_fine_intf
        INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:,:) ::px_interface_det,ux_interface_det, vx_interface_det, wx_interface_det
        INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:,:) ::pz_interface_det,uz_interface_det, vz_interface_det, wz_interface_det
        INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:,:) ::py_interface_det,uy_interface_det, vy_interface_det, wy_interface_det
        REAL (KIND=8)  :: a_blk, b_blk,xintf_start, xintf_end,yintf_start, yintf_end, zintf_start, zintf_end
        REAL (KIND=8)  :: xintf_st_new, xintf_en_new, yintf_st_new, yintf_en_new
        REAL (KIND=8)  :: zintf_st_new, zintf_en_new
        INTEGER (KIND=8)   :: a_msh, b_msh, a_intf, b_intf
        INTEGER (KIND=8)   :: counterxu,counteryu,counterxp,counteryp,counterxv,counteryv
        INTEGER (KIND=8)   :: counterxw,counteryw
        INTEGER (KIND=8)   :: counterzu,counterzp,counterzv, counterzw
        INTEGER (KIND=8)   :: cpx_start, cpy_start,cpx_end, cpy_end 
        INTEGER (KIND=8)   :: cux_start, cuy_start,cux_end, cuy_end 
        INTEGER (KIND=8)   :: cvx_start, cvy_start,cvx_end, cvy_end 
        INTEGER (KIND=8)   :: ccellx_start, ccelly_start,ccellx_end, ccelly_end 
        INTEGER (KIND=8)   :: fvx_start, fvy_start 
        INTEGER (KIND=8)   :: fpx_start, fpy_start 
        INTEGER (KIND=8)   :: fux_start, fuy_start 
        INTEGER (KIND=8)   :: fcellx_start, fcelly_start 
        end type Interfaces

        type(Blocks),allocatable ::block(:)
        type(Interfaces),allocatable ::intfr(:)


       !REAL (KIND=8)  :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4, bl_interp_ans
   
      REAL (KIND=8) ::d_fine_x1,d_fine_y1,d_coarse_x1,d_coarse_y1,d_val_x1,d_val_y1,d_fine_frac,d_coarse_frac
      REAL (KIND=8) :: l_intp_valx, l_intp_x1,l_intp_x2,l_intp_y1, l_intp_y2, l_intp_valy
    

        INTEGER(KIND=8)  :: nfl_blk
       INTEGER (KIND=8),ALLOCATABLE,DIMENSION(:) :: fl_blk

       ! ibm variables
       INTEGER  :: bcType
       REAL (KIND = 8)    :: aoa, piv_pt, var_surf, aoa1
       REAL (KIND = 8)    :: phase_angle, aoa2, a0y, ang_theta, alpha_t, theta_t
       !INTEGER (KIND=8)   :: surGeoPoints, ibElems, ibNodes 
       INTEGER (KIND=8)   :: surGeoPoints 
      !INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:)   :: ibElP1, ibElP2, ibElP3
      !REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: xnode, ynode,  znode, xnode1, ynode1, znode1, bcSurf
                                                      
        !amgx
        INTEGER (KIND=4) :: nnz_max,nu_max
!       REAL (KIND = 8), ALLOCATABLE,TARGET, DIMENSION (:)::var_data,var_rhs,var_sol,var_row,var_col
        type(c_ptr)::cptr_crs,cptr_data,cptr_col,cptr_row,cptr_rhs,cptr_sol,cptr_nit, cptr_nit1
        INTEGER(KIND=4) ::init_stat,dest_stat,solve_stat,amgx_checker
       !!$acc commands..
      !!$acc declare create (deltax, deltay, deltaz)
      !!$acc declare create (x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp)
      !!$acc declare create (b, Acx, Acy, Acz)
      !!$acc declare create (pc, pco)
      !!$acc declare create (xnode, ynode, znode, ibElP1, ibElP2, ibElP3)
      !!$acc declare create (xnode1, ynode1, znode1)
      !!$acc declare create (xcent, ycent, zcent)
      !!$acc declare create (cosAlpha, cosBeta, cosGamma)
      !!$acc declare create (alpha3, beta3, gamma3)
      !!$acc declare create (cell)
      !!$acc declare create (minElemcell)
      !!$acc declare create (fluidIndexPtr, redCellIndexPtr, blackCellIndexPtr)
      !!$acc declare create (interceptedIndexPtr, solidIndexPtr)
      !!$acc declare create (pNormDis, u1NormDis, u2NormDis, v1NormDis, &
      !!$acc                           v2NormDis, w1NormDis, w2NormDis)
      !!$acc declare create (nelp, nelu1, nelu2, nelv1, nelv2, nelw1, nelw2)
      !!$acc declare create (u, ut, v, vt, w, wt, p)
      !!$acc declare create (resi_u, resi_v, resi_w)
      !!$acc declare create (                                &
      !!$acc ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
      !!$acc ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
      !!$acc ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
      !!$acc ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
      !!$acc ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
      !!$acc ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww, &
      !!$acc ca1_uv, ca2_uv, ca3_uv, ca4_uv, ca5_uv, ca6_uv, &
      !!$acc ck1_uv, ck2_uv, ck3_uv, ck4_uv, ck5_uv, ck6_uv, &
      !!$acc ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
      !!$acc ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
      !!$acc ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
      !!$acc ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
      !!$acc ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
      !!$acc ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
      !!$acc ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
      !!$acc ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
      !!$acc ca1_wv, ca2_wv, ca3_wv, ca4_wv, ca5_wv, ca6_wv, &
      !!$acc ck1_wv, ck2_wv, ck3_wv, ck4_wv, ck5_wv, ck6_wv)         
      !!$acc declare create (Elemcell, ucell, vcell, wcell, pcell, areaElem, stressElem)      
       END MODULE global
                
                
                
                

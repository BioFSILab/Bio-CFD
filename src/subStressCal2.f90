!csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
       SUBROUTINE stressCal2
       USE global
       IMPLICIT NONE
       INTEGER, PARAMETER :: rk = selected_real_kind(8)

       INTEGER:: i, j, k, ielem, inode, i_x1, i_y1, i_z1, g, nvib,n

       INTEGER:: i_cell, j_cell, k_cell

       REAL(KIND = 8):: diagdis, normdis, aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z

       REAL(KIND = 8):: xsurf, ysurf, zsurf, pos1_x, pos1_y, pos1_z,              &
                        psurf, p_pos1, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2,  &
                        p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,     &
                        p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e, &
                        usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,        &
                        vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,        &
                        wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,        &
                        dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e,   &
                        dvdz_e, dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, &
                        u_x1_z2, u_x2_z2, u_z1_x1, u_z2_x1, u_z1_x2, u_z2_x2,     &
                        v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1,     &
                        v_z1_x2, v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2,     &
                        w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, dudn_s, dvdn_s, dwdn_s

       !REAL(KIND = 8):: alen, area

       !REAL(KIND = 8), DIMENSION(ibNodes):: Anode, Afnode, stressNode, TAWSS, OSI

       !REAL(KIND = 8), DIMENSION(ibElems):: stressElem

       CHARACTER(LEN=150):: filename1
        REAL(KIND = 8):: modStressNode, modSIGNWSS

       REAL(KIND = 8):: alen, area, area_xz, area_yz, area_xy

       REAL(KIND = 8):: shear_x_force, shear_y_force, shear_z_force,  &
                        f_surf, f_surf_x, f_surf_y, f_surf_z

       REAL(KIND = 8):: pressureDrag, viscousDrag, viscousLift,                        &
                        pressureLift, viscousDragcoefficient, pressureDragcoefficient, &
                        viscousLiftcoefficient, PressureLiftcoefficient, area_Sx,      &
                        area_Sy, surf_area, total_v_x, total_p_x, total_v_y, total_p_y
        REAL (KIND=8)      :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z,lenEL
        REAL (KIND=8)      :: var_xcent, var_ycent, var_zcent,ylim1,ylim2,zlim1

       DO g=blk_start, nblocks
       viscousDrag = 0.
       pressureDrag = 0.
       viscousLift = 0.
       PressureLift = 0.
       surf_area = 0.
	area_Sx = 0.
       area_Sy = 0
       block(g)%stressElem = 0._rk
       block(g)%Total_V_Fx = 0._rk
       block(g)%Total_P_Fx = 0._rk

       block(g)%Total_V_Fy = 0._rk
       block(g)%Total_P_Fy = 0._rk
       total_v_x=0._rk
       total_p_y=0._rk
       total_v_y=0._rk
       total_p_x=0._rk
       !!!$acc parallel loop gang vector reduction(+: pressureDrag, viscousDrag, viscousLift, pressureLift, area_Sx, area_Sy, surf_area)  &
       !!!$acc private (i_x1, i_y1, i_z1, i_cell, j_cell, k_cell, diagdis, normdis,                    &
       !!!$acc          aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z,                       &
       !!!$acc          xsurf, ysurf, zsurf, pos1_x, pos1_y, pos1_z,                                   &
       !!!$acc          psurf, p_pos1, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2,                       &
       !!!$acc          p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,                          &
       !!!$acc          p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e,                      &
       !!!$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,                             &
       !!!$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,                             &
       !!!$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,                             &
       !!!$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e,                        &
       !!!$acc          dvdz_e, dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1,                      &
       !!!$acc          u_x1_z2, u_x2_z2, u_z1_x1, u_z2_x1, u_z1_x2, u_z2_x2,                          &
       !!!$acc          v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1,                          &
       !!!$acc          v_z1_x2, v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2,                          &
       !!!$acc          w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, dudn_s, dvdn_s, dwdn_s,                    &
       !!!$acc          alen, area, area_xz, area_yz, area_xy, shear_x_force, shear_y_force,           &
       !!!$acc          shear_z_force, f_surf, f_surf_x, f_surf_y, f_surf_z, viscousDragcoefficient,   &
       !!!$acc          pressureDragcoefficient,viscousLiftcoefficient, PressureLiftcoefficient)       &
       !!!$acc present xcent, ycent, zcent, x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp, &
       !!!$acc          ut, u, vt, v, wt, w, p, deltax, deltay, deltaz, cosAlpha, cosBeta, cosGamma,     &
       !!!$acc          block(g)%alpha3, block(g)%beta3, block(g)%gamma3, xnode1, ynode1, znode1, Elemcell, ucell, vcell, wcell, pcell, areaElem, stressElem)    &
       !!!$acc firstprivate (nx, ny, nz, deltat)
        !!!$acc parallel loop gang vector   &
         !$acc parallel loop gang vector reduction(+: pressureDrag, viscousDrag, viscousLift, pressureLift, area_Sx, area_Sy, surf_area,total_v_x,total_p_x,total_v_y,total_p_y)  &
         !$acc private (i_x1, i_y1, i_z1, i_cell, j_cell, k_cell, diagdis, normdis,                    &
         !$acc          aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z,                       &
         !$acc          xsurf, ysurf, zsurf, pos1_x, pos1_y, pos1_z,                                   &
         !$acc          psurf, p_pos1, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2,                       &
         !$acc          p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,                          &
         !$acc          p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e,                      &
         !$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,                             &
         !$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,                             &
         !$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,                             &
         !$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e,                        &
         !$acc          dvdz_e, dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1,                      &
         !$acc          u_x1_z2, u_x2_z2, u_z1_x1, u_z2_x1, u_z1_x2, u_z2_x2,                          &
         !$acc          v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1,                          &
         !$acc          v_z1_x2, v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2,                          &
         !$acc          w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, dudn_s, dvdn_s, dwdn_s,                    &
         !$acc          alen, area, area_xz, area_yz, area_xy, shear_x_force, shear_y_force,           &
         !$acc          shear_z_force, f_surf, f_surf_x, f_surf_y, f_surf_z)   &
         !$acc default(present) &
         !$acc firstprivate (block(g)%nx, block(g)%ny, block(g)%nz, deltat)
        DO ielem = 1, block(g)%ibElems
       !IF((block(g)%ycent(ielem).ge.-3.30).and.(block(g)%ycent(ielem).le.18.03)) THEN
!***********************interpolation points****************************
       i_cell = block(g)%Elemcell(ielem,1)
       j_cell = block(g)%Elemcell(ielem,2)
       k_cell = block(g)%Elemcell(ielem,3)

       xsurf = block(g)%xcent(ielem)
       ysurf = block(g)%ycent(ielem)
       zsurf = block(g)%zcent(ielem)

	del_X = block(g)%x1(i_cell+1)-block(g)%x1(i_cell)
       del_Y = block(g)%y1(j_cell+1)-block(g)%y1(j_cell)
       del_Z = block(g)%z1(k_cell+1)-block(g)%z1(k_cell)

       diagdis = dsqrt(del_X**2 + del_Y**2 + del_Z**2)
       normdis = diagdis
	pos1_x = xsurf + normdis*block(g)%cosAlpha(ielem)
	pos1_y = ysurf + normdis*block(g)%cosBeta(ielem)
	pos1_z = zsurf + normdis*block(g)%cosGamma(ielem)

!**************************velocity at the surface**********************
       IF (block(g)%ibSurfID(ielem)==50) THEN
        block(g)% thetaDot  =  0.
        block(g)% thetaDDot =  0.
         usurf     =  0.
	      vsurf     =  0.
	      wsurf     =  0.
	      ac_z      =  0.  !-thetaDot**2*(zcent(nelp(index_ts(n))) - piv_z)
         ac_y      =  0.  !-thetaDot**2*(ycent(nelp(index_ts(n))) - piv_y)
         at_z      =  0.  ! thetaDDot*(ycent(nelp(index_ts(n))) - piv_y)
         at_y      =  0.  !
       ELSEIF (block(g)%ibSurfID(ielem)==51) THEN
        block(g)% thetaDot  = block(g)% thetaDot1
	     block(g)% thetaDDot = block(g)% thetaDDot1
         !usurf     =  0.
         !vsurf     = -thetaDot*(zcent(ielem) - piv_z)
	      !wsurf     =  thetaDot*(ycent(ielem) - piv_y)
	      ac_z      = -block(g)%thetaDot**2*(block(g)%zcent(ielem) -block(g)% piv_z)
         ac_y      = -block(g)%thetaDot**2*(block(g)%ycent(ielem) -block(g)% piv_y)
         at_z      = block(g)% thetaDDot*(block(g)%ycent(ielem)-block(g)% piv_y)
         at_y      = -block(g)%thetaDDot*(block(g)%zcent(ielem)-block(g)% piv_z)
       ELSEIF (block(g)%ibSurfId(ielem)==52) THEN
	     block(g)% thetaDot  = block(g)% thetaDot2
	     block(g)% thetaDDot = block(g)% thetaDDot2
         !usurf     =  0.
	      !wsurf     =  thetaDot*(ycent(ielem) - piv_y)
         !vsurf     = -thetaDot*(zcent(ielem) - piv_z)
	      ac_z      = -block(g)%thetaDot**2*(block(g)%zcent(ielem) -block(g)% piv_z)
         ac_y      = -block(g)%thetaDot**2*(block(g)%ycent(ielem) - block(g)%piv_y)
         at_z      =  block(g)%thetaDDot*(block(g)%ycent(ielem) -block(g)% piv_y)
         at_y      = -block(g)%thetaDDot*(block(g)%zcent(ielem) - block(g)%piv_z)
       ENDIF
       usurf     =  0.
       wsurf     = block(g)% thetaDot*(block(g)%ycent(ielem) -block(g)% piv_y)
       vsurf     = -block(g)%thetaDot*(block(g)%zcent(ielem) -block(g)% piv_z) + block(g)%ydot
       dpdn      = -((ac_z + at_z)*block(g)%cosAlpha(ielem)  + (ac_y + at_y)*block(g)%cosBeta(ielem)) -block(g)%yddot*block(g)%cosBeta(ielem)
!*******************velocity interpolation at point 2******************

!******************u velocity interpolation at point 2******************
        i_x1 = block(g)%ucell(ielem,1)
        i_y1 = block(g)%ucell(ielem,2)
        i_z1 = block(g)%ucell(ielem,3)

       !interpolation along x @ z1 plane
         u_x1_z1 = block(g)%ut(i_x1-1, i_y1, i_z1)   + (block(g)%ut(i_x1, i_y1, i_z1)   - block(g)%ut(i_x1-1, i_y1, i_z1))  *(pos1_x - block(g)%xu(i_x1))/(block(g)%xu(i_x1+1) - block(g)%xu(i_x1))
         u_x2_z1 = block(g)%ut(i_x1-1, i_y1+1, i_z1) + (block(g)%ut(i_x1, i_y1+1, i_z1) - block(g)%ut(i_x1-1, i_y1+1, i_z1))*(pos1_x - block(g)%xu(i_x1))/(block(g)%xu(i_x1+1) - block(g)%xu(i_x1))

         !interpolation along x @ z2 plane
         u_x1_z2 = block(g)%ut(i_x1-1, i_y1, i_z1+1)   + (block(g)%ut(i_x1, i_y1, i_z1+1)   - block(g)%ut(i_x1-1, i_y1, i_z1+1))  *(pos1_x - block(g)%xu(i_x1))/(block(g)%xu(i_x1+1) - block(g)%xu(i_x1))
         u_x2_z2 = block(g)%ut(i_x1-1, i_y1+1, i_z1+1) + (block(g)%ut(i_x1, i_y1+1, i_z1+1) - block(g)%ut(i_x1-1, i_y1+1, i_z1+1))*(pos1_x - block(g)%xu(i_x1))/(block(g)%xu(i_x1+1) - block(g)%xu(i_x1))

         !This will be used for dudz calculation point 1
         u_z1 = u_x1_z1 + (u_x2_z1 - u_x1_z1)*(pos1_y - block(g)%yu(i_y1))/(block(g)%yu(i_y1+1)-block(g)%yu(i_y1))
         u_z2 = u_x1_z2 + (u_x2_z2 - u_x1_z2)*(pos1_y - block(g)%yu(i_y1))/(block(g)%yu(i_y1+1)-block(g)%yu(i_y1))

         !This will be used for dudy calculation at point 1
         u_y1 = u_x1_z1 + (u_x1_z2 - u_x1_z1)*(pos1_z - block(g)%zu(i_z1))/(block(g)%zu(i_z1+1)-block(g)%zu(i_z1))
         u_y2 = u_x2_z1 + (u_x2_z2 - u_x2_z1)*(pos1_z - block(g)%zu(i_z1))/(block(g)%zu(i_z1+1)-block(g)%zu(i_z1))

         !interpolation along z @ x1 plane
         u_z1_x1 = block(g)%ut(i_x1-1, i_y1, i_z1)   + (block(g)%ut(i_x1-1, i_y1, i_z1+1) - block(g)%ut(i_x1-1, i_y1, i_z1))  *(pos1_z - block(g)%zu(i_z1))/(block(g)%zu(i_z1+1) - block(g)%zu(i_z1))
         u_z2_x1 = block(g)%ut(i_x1-1, i_y1+1, i_z1)   + (block(g)%ut(i_x1-1, i_y1+1, i_z1+1) - block(g)%ut(i_x1-1, i_y1+1, i_z1))  *(pos1_z - block(g)%zu(i_z1))/(block(g)%zu(i_z1+1) - block(g)%zu(i_z1))

         !interpolation along z @ x2 plane
         u_z1_x2 = block(g)%ut(i_x1, i_y1, i_z1)   + (block(g)%ut(i_x1, i_y1, i_z1+1) - block(g)%ut(i_x1, i_y1, i_z1))  *(pos1_z - block(g)%zu(i_z1))/(block(g)%zu(i_z1+1) - block(g)%zu(i_z1))
         u_z2_x2 = block(g)%ut(i_x1, i_y1+1, i_z1)   + (block(g)%ut(i_x1, i_y1+1, i_z1+1) - block(g)%ut(i_x1, i_y1+1, i_z1))  *(pos1_z - block(g)%zu(i_z1))/(block(g)%zu(i_z1+1) - block(g)%zu(i_z1))

         !This will be used for dudx calculation at point 1
         u_x1 = u_z1_x1 + (u_z2_x1 - u_z1_x1)*(pos1_y - block(g)%yu(i_y1))/(block(g)%yu(i_y1+1)-block(g)%yu(i_y1))
         u_x2 = u_z1_x2 + (u_z2_x2 - u_z1_x2)*(pos1_y - block(g)%yu(i_y1))/(block(g)%yu(i_y1+1)-block(g)%yu(i_y1))

         u_pos1 = u_x1 + (u_x2 - u_x1)*(pos1_x - block(g)%xu(i_x1))/(block(g)%xu(i_x1+1)-block(g)%xu(i_x1))

         h2 = dabs(block(g)%xu(i_x1+1) - pos1_x)
         h1 = dabs(block(g)%xu(i_x1)   - pos1_x)
         dudx_e = (h1**2*u_x2 - h2**2*u_x1 + (h2**2- h1**2)*u_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%yu(i_y1+1) - pos1_y)
         h1 = dabs(block(g)%yu(i_y1)   - pos1_y)
         dudy_e = (h1**2*u_y2 - h2**2*u_y1 + (h2**2- h1**2)*u_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%zu(i_z1+1) - pos1_z)
         h1 = dabs(block(g)%zu(i_z1)   - pos1_z)
         dudz_e = (h1**2*u_z2 - h2**2*u_z1 + (h2**2- h1**2)*u_pos1)/(h1*h2*(h1+h2)+1e-16)

         dudn_e = dudx_e*block(g)%cosAlpha(ielem) + dudy_e*block(g)%cosBeta(ielem) + dudz_e*block(g)%cosGamma(ielem)

         dudn_s = (2./normdis)*(u_pos1 - usurf) - dudn_e

!******************v velocity interpolation in point 2******************
        i_x1 = block(g)%vcell(ielem,1)
        i_y1 = block(g)%vcell(ielem,2)
        i_z1 = block(g)%vcell(ielem,3)

       !interpolation along x @ z1 plane
         v_x1_z1 = block(g)%vt(i_x1, i_y1-1, i_z1)   + (block(g)%vt(i_x1+1, i_y1-1, i_z1) - block(g)%vt(i_x1, i_y1-1, i_z1))*(pos1_x - block(g)%xv(i_x1))/(block(g)%xv(i_x1+1) - block(g)%xv(i_x1))
         v_x2_z1 = block(g)%vt(i_x1, i_y1, i_z1)     + (block(g)%vt(i_x1+1, i_y1, i_z1)   - block(g)%vt(i_x1, i_y1, i_z1))*(pos1_x - block(g)%xv(i_x1))/(block(g)%xv(i_x1+1) - block(g)%xv(i_x1))

         !interpolation along x @ z2 plane
         v_x1_z2 = block(g)%vt(i_x1, i_y1-1, i_z1+1)   + (block(g)%vt(i_x1+1, i_y1-1, i_z1+1) - block(g)%vt(i_x1, i_y1-1, i_z1+1))*(pos1_x - block(g)%xv(i_x1))/(block(g)%xv(i_x1+1) - block(g)%xv(i_x1))
         v_x2_z2 = block(g)%vt(i_x1, i_y1, i_z1+1)     + (block(g)%vt(i_x1+1, i_y1, i_z1+1)   - block(g)%vt(i_x1, i_y1, i_z1+1))*(pos1_x - block(g)%xv(i_x1))/(block(g)%xv(i_x1+1) - block(g)%xv(i_x1))

         !This will be used for dvdz calculation point 1
         v_z1 = v_x1_z1 + (v_x2_z1 - v_x1_z1)*(pos1_y - block(g)%yv(i_y1))/(block(g)%yv(i_y1+1)-block(g)%yv(i_y1))
         v_z2 = v_x1_z2 + (v_x2_z2 - v_x1_z2)*(pos1_y - block(g)%yv(i_y1))/(block(g)%yv(i_y1+1)-block(g)%yv(i_y1))

         !This will be used for dvdy calculation at point 1
         v_y1 = v_x1_z1 + (v_x1_z2 - v_x1_z1)*(pos1_z - block(g)%zv(i_z1))/(block(g)%zv(i_z1+1)-block(g)%zv(i_z1))
         v_y2 = v_x2_z1 + (v_x2_z2 - v_x2_z1)*(pos1_z - block(g)%zv(i_z1))/(block(g)%zv(i_z1+1)-block(g)%zv(i_z1))

         !interpolation along z @ x1 plane
         v_z1_x1 = block(g)%vt(i_x1, i_y1-1, i_z1)   + (block(g)%vt(i_x1, i_y1-1, i_z1+1) - block(g)%vt(i_x1, i_y1-1, i_z1))  *(pos1_z - block(g)%zv(i_z1))/(block(g)%zv(i_z1+1) - block(g)%zv(i_z1))
         v_z2_x1 = block(g)%vt(i_x1, i_y1, i_z1)   + (block(g)%vt(i_x1, i_y1, i_z1+1) - block(g)%vt(i_x1, i_y1, i_z1))  *(pos1_z - block(g)%zv(i_z1))/(block(g)%zv(i_z1+1) - block(g)%zv(i_z1))

         !interpolation along z @ x2 plane
         v_z1_x2 = block(g)%vt(i_x1+1, i_y1-1, i_z1)   + (block(g)%vt(i_x1+1, i_y1-1, i_z1+1) - block(g)%vt(i_x1+1, i_y1-1, i_z1))  *(pos1_z - block(g)%zv(i_z1))/(block(g)%zv(i_z1+1) - block(g)%zv(i_z1))
         v_z2_x2 = block(g)%vt(i_x1+1, i_y1, i_z1)   + (block(g)%vt(i_x1+1, i_y1, i_z1+1) - block(g)%vt(i_x1+1, i_y1, i_z1))  *(pos1_z - block(g)%zv(i_z1))/(block(g)%zv(i_z1+1) - block(g)%zv(i_z1))

         !This will be used for dvdx calculation at point 1
         v_x1 = v_z1_x1 + (v_z2_x1 - v_z1_x1)*(pos1_y - block(g)%yv(i_y1))/(block(g)%yv(i_y1+1)-block(g)%yv(i_y1))
         v_x2 = v_z1_x2 + (v_z2_x2 - v_z1_x2)*(pos1_y - block(g)%yv(i_y1))/(block(g)%yv(i_y1+1)-block(g)%yv(i_y1))

         v_pos1 = v_x1 + (v_x2 - v_x1)*(pos1_x - block(g)%xv(i_x1))/(block(g)%xv(i_x1+1)-block(g)%xv(i_x1))

         h2 = dabs(block(g)%xv(i_x1+1) - pos1_x)
         h1 = dabs(block(g)%xv(i_x1)   - pos1_x)
         dvdx_e = (h1**2*v_x2 - h2**2*v_x1 + (h2**2- h1**2)*v_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%yv(i_y1+1) - pos1_y)
         h1 = dabs(block(g)%yv(i_y1)   - pos1_y)
         dvdy_e = (h1**2*v_y2 - h2**2*v_y1 + (h2**2- h1**2)*v_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%zv(i_z1+1) - pos1_z)
         h1 = dabs(block(g)%zv(i_z1)   - pos1_z)
         dvdz_e = (h1**2*v_z2 - h2**2*v_z1 + (h2**2- h1**2)*v_pos1)/(h1*h2*(h1+h2)+1e-16)


         dvdn_e = dvdx_e*block(g)%cosAlpha(ielem) + dvdy_e*block(g)%cosBeta(ielem) +  dvdz_e*block(g)%cosGamma(ielem)

	  dvdn_s = (2./normdis)*(v_pos1 - vsurf) - dvdn_e

!******************w velocity interpolation in point 2******************
        i_x1 = block(g)%wcell(ielem,1)
        i_y1 = block(g)%wcell(ielem,2)
        i_z1 = block(g)%wcell(ielem,3)

       !interpolation along x @ z1 plane
         w_x1_z1 = block(g)%wt(i_x1, i_y1, i_z1-1)   + (block(g)%wt(i_x1+1, i_y1, i_z1-1) - block(g)%wt(i_x1, i_y1, i_z1-1))*(pos1_x - block(g)%xw(i_x1))/(block(g)%xw(i_x1+1) - block(g)%xw(i_x1))
         w_x2_z1 = block(g)%wt(i_x1, i_y1+1, i_z1-1)   + (block(g)%wt(i_x1+1, i_y1+1, i_z1-1) - block(g)%wt(i_x1, i_y1+1, i_z1-1))*(pos1_x - block(g)%xw(i_x1))/(block(g)%xw(i_x1+1) - block(g)%xw(i_x1))
         !interpolation along x @ z2 plane
         w_x1_z2 = block(g)%wt(i_x1, i_y1, i_z1)   + (block(g)%wt(i_x1+1, i_y1, i_z1) - block(g)%wt(i_x1, i_y1, i_z1))*(pos1_x - block(g)%xw(i_x1))/(block(g)%xw(i_x1+1) - block(g)%xw(i_x1))
         w_x2_z2 = block(g)%wt(i_x1, i_y1+1, i_z1)   + (block(g)%wt(i_x1+1, i_y1+1, i_z1) - block(g)%wt(i_x1, i_y1+1, i_z1))*(pos1_x - block(g)%xw(i_x1))/(block(g)%xw(i_x1+1) - block(g)%xw(i_x1))

         !This will be used for dwdz calculation point 1
         w_z1 = w_x1_z1 + (w_x2_z1 - w_x1_z1)*(pos1_y - block(g)%yw(i_y1))/(block(g)%yw(i_y1+1)-block(g)%yw(i_y1))
         w_z2 = w_x1_z2 + (w_x2_z2 - w_x1_z2)*(pos1_y - block(g)%yw(i_y1))/(block(g)%yw(i_y1+1)-block(g)%yw(i_y1))

         !This will be used for dwdy calculation at point 1
         w_y1 = w_x1_z1 + (w_x1_z2 - w_x1_z1)*(pos1_z - block(g)%zw(i_z1))/(block(g)%zw(i_z1+1)-block(g)%zw(i_z1))
         w_y2 = w_x2_z1 + (w_x2_z2 - w_x2_z1)*(pos1_z - block(g)%zw(i_z1))/(block(g)%zw(i_z1+1)-block(g)%zw(i_z1))

         !interpolation along z @ x1 plane
         w_z1_x1 = block(g)%wt(i_x1, i_y1, i_z1-1)   + (block(g)%wt(i_x1, i_y1, i_z1) - block(g)%wt(i_x1, i_y1, i_z1-1))  *(pos1_z - block(g)%zw(i_z1))/(block(g)%zw(i_z1+1) - block(g)%zw(i_z1))
         w_z2_x1 = block(g)%wt(i_x1, i_y1+1, i_z1-1)   + (block(g)%wt(i_x1, i_y1+1, i_z1) - block(g)%wt(i_x1, i_y1+1, i_z1-1))  *(pos1_z - block(g)%zw(i_z1))/(block(g)%zw(i_z1+1) - block(g)%zw(i_z1))

         !interpolation along z @ x2 plane
         w_z1_x2 = block(g)%wt(i_x1+1, i_y1, i_z1-1)   + (block(g)%wt(i_x1+1, i_y1, i_z1) - block(g)%wt(i_x1+1, i_y1, i_z1-1))  *(pos1_z - block(g)%zw(i_z1))/(block(g)%zw(i_z1+1) - block(g)%zw(i_z1))
         w_z2_x2 = block(g)%wt(i_x1+1, i_y1+1, i_z1-1)   + (block(g)%wt(i_x1+1, i_y1+1, i_z1) - block(g)%wt(i_x1+1, i_y1+1, i_z1-1))  *(pos1_z - block(g)%zw(i_z1))/(block(g)%zw(i_z1+1) - block(g)%zw(i_z1))

         !This will be used for dwdx calculation at point 1
         w_x1 = w_z1_x1 + (w_z2_x1 - w_z1_x1)*(pos1_y - block(g)%yw(i_y1))/(block(g)%yw(i_y1+1)-block(g)%yw(i_y1))
         w_x2 = w_z1_x2 + (w_z2_x2 - w_z1_x2)*(pos1_y - block(g)%yw(i_y1))/(block(g)%yw(i_y1+1)-block(g)%yw(i_y1))

         w_pos1 = w_x1 + (w_x2 - w_x1)*(pos1_x - block(g)%xw(i_x1))/(block(g)%xw(i_x1+1)-block(g)%xw(i_x1))

         h2 = dabs(block(g)%xw(i_x1+1) - pos1_x)
         h1 = dabs(block(g)%xw(i_x1)   - pos1_x)
         dwdx_e = (h1**2*w_x2 - h2**2*w_x1 + (h2**2- h1**2)*w_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%yw(i_y1+1) - pos1_y)
         h1 = dabs(block(g)%yw(i_y1)   - pos1_y)
         dwdy_e = (h1**2*w_y2 - h2**2*w_y1 + (h2**2- h1**2)*w_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%zw(i_z1+1) - pos1_z)
         h1 = dabs(block(g)%zw(i_z1)   - pos1_z)
         dwdz_e = (h1**2*w_z2 - h2**2*w_z1 + (h2**2- h1**2)*w_pos1)/(h1*h2*(h1+h2)+1e-16)


         dwdn_e = dwdx_e*block(g)%cosAlpha(ielem) + dwdy_e*block(g)%cosBeta(ielem) +  dwdz_e*block(g)%cosGamma(ielem)

         dwdn_s = (2./normdis)*(w_pos1 - wsurf) - dwdn_e

!***********************calculate area of the elements******************
       alen = sqrt(block(g)%alpha3(ielem)**2.+block(g)%beta3(ielem)**2.+ block(g)%gamma3(ielem)**2.)
       area = alen/2.
       area_yz = 0.5*abs(block(g)%alpha3(ielem))
       area_xz = 0.5*abs(block(g)%beta3(ielem))
       area_xy = 0.5*abs(block(g)%gamma3(ielem))
       !areaElem(ielem) = area

!*********non-dimensional viscous stress & force calculation************
       stx1 = dudn_s - (dudn_s*block(g)%cosAlpha(ielem) + dvdn_s*block(g)%cosBeta(ielem) + dwdn_s*block(g)%cosGamma(ielem))*block(g)%cosAlpha(ielem)

       sty1 = dvdn_s - (dudn_s*block(g)%cosAlpha(ielem) + dvdn_s*block(g)%cosBeta(ielem) + dwdn_s*block(g)%cosGamma(ielem))*block(g)%cosBeta(ielem)

       stz1 = dwdn_s - (dudn_s*block(g)%cosAlpha(ielem) + dvdn_s*block(g)%cosBeta(ielem) + dwdn_s*block(g)%cosGamma(ielem))*block(g)%cosGamma(ielem)

       shear_x_force = (1.0/re)*stx1*area*block(g)%rhof

       shear_y_force = (1.0/re)*sty1*area*block(g)%rhof

       shear_z_force = (1.0/re)*stz1*area*block(g)%rhof

      !shear_x_force = (1.0/re)*stx1*area

      !shear_y_force = (1.0/re)*sty1*area
      !
      !shear_z_force = (1.0/re)*stz1*area
       !stressElem(ielem) = 0.0035*stx1
      !block(g)%stressElem(ielem, 1) = (1.0/re)*stx1*block(g)%rhof
      !block(g)%stressElem(ielem, 2) = (1.0/re)*sty1*block(g)%rhof
      !block(g)%stressElem(ielem, 3) = (1.0/re)*stz1*block(g)%rhof
       block(g)%stressElem(ielem, 1) = (1.0/re)*stx1
       block(g)%stressElem(ielem, 2) = (1.0/re)*sty1
       block(g)%stressElem(ielem, 3) = (1.0/re)*stz1
!************************presssure interpolation************************

!*******************pressure interpolation at point 2*******************
       !$acc loop seq
       DO i = 2, block(g)%nx+1
       if(pos1_x>=block(g)%xp(i).and.pos1_x<block(g)%xp(i+1)) i_x1 = i
       END DO
	!$acc loop seq
       DO j = 2, block(g)%ny+1
       if(pos1_y>=block(g)%yp(j).and.pos1_y<block(g)%yp(j+1)) i_y1 = j
       END DO
	!$acc loop seq
       DO k = 1, block(g)%nz+2
       if(pos1_z>=block(g)%zp(k).and.pos1_z<block(g)%zp(k+1)) i_z1 = k
       END DO

       block(g)%pcell(ielem,1) = i_x1
       block(g)%pcell(ielem,2) = i_y1
       block(g)%pcell(ielem,3) = i_z1

        !interpolation along x  @ z1 plane
         p_x1_z1 = block(g)%p(i_x1, i_y1, i_z1)   + (block(g)%p(i_x1+1, i_y1, i_z1)   - block(g)%p(i_x1, i_y1, i_z1))  *(pos1_x - block(g)%xp(i_x1))/(block(g)%xp(i_x1+1) - block(g)%xp(i_x1))
         p_x2_z1 = block(g)%p(i_x1, i_y1+1, i_z1) + (block(g)%p(i_x1+1, i_y1+1, i_z1) - block(g)%p(i_x1, i_y1+1, i_z1))*(pos1_x - block(g)%xp(i_x1))/(block(g)%xp(i_x1+1) - block(g)%xp(i_x1))

         !interpolation along x  @ z2 plane
         p_x1_z2 = block(g)%p(i_x1, i_y1, i_z1+1)   + (block(g)%p(i_x1+1, i_y1, i_z1+1)   - block(g)%p(i_x1, i_y1, i_z1+1))  *(pos1_x - block(g)%xp(i_x1))/(block(g)%xp(i_x1+1) - block(g)%xp(i_x1))
         p_x2_z2 = block(g)%p(i_x1, i_y1+1, i_z1+1) + (block(g)%p(i_x1+1, i_y1+1, i_z1+1) - block(g)%p(i_x1, i_y1+1, i_z1+1))*(pos1_x - block(g)%xp(i_x1))/(block(g)%xp(i_x1+1) - block(g)%xp(i_x1))

         !This will be used for dpdz calculation point 1
         p_z1 = p_x1_z1 + (p_x2_z1 - p_x1_z1)*(pos1_y - block(g)%yp(i_y1))/(block(g)%yp(i_y1+1)-block(g)%yp(i_y1))
         p_z2 = p_x1_z2 + (p_x2_z2 - p_x1_z2)*(pos1_y - block(g)%yp(i_y1))/(block(g)%yp(i_y1+1)-block(g)%yp(i_y1))

         !This will be used for dpdy calculation at point 1
         p_y1 = p_x1_z1 + (p_x1_z2 - p_x1_z1)*(pos1_z - block(g)%zp(i_z1))/(block(g)%zp(i_z1+1)-block(g)%zp(i_z1))
         p_y2 = p_x2_z1 + (p_x2_z2 - p_x2_z1)*(pos1_z - block(g)%zp(i_z1))/(block(g)%zp(i_z1+1)-block(g)%zp(i_z1))

         !interpolation along z @ x1 plane
         p_z1_x1 = block(g)%p(i_x1, i_y1, i_z1)   + (block(g)%p(i_x1, i_y1, i_z1+1) - block(g)%p(i_x1, i_y1, i_z1))  *(pos1_z - block(g)%zp(i_z1))/(block(g)%zp(i_z1+1) - block(g)%zp(i_z1))
         p_z2_x1 = block(g)%p(i_x1, i_y1+1, i_z1)   + (block(g)%p(i_x1, i_y1+1, i_z1+1) - block(g)%p(i_x1, i_y1+1, i_z1))  *(pos1_z - block(g)%zp(i_z1))/(block(g)%zp(i_z1+1) - block(g)%zp(i_z1))

         !interpolation along z @ x2 plane
         p_z1_x2 = block(g)%p(i_x1+1, i_y1, i_z1)   + (block(g)%p(i_x1+1, i_y1, i_z1+1) - block(g)%p(i_x1+1, i_y1, i_z1))  *(pos1_z - block(g)%zp(i_z1))/(block(g)%zp(i_z1+1) - block(g)%zp(i_z1))
         p_z2_x2 = block(g)%p(i_x1+1, i_y1+1, i_z1)   + (block(g)%p(i_x1+1, i_y1+1, i_z1+1) - block(g)%p(i_x1+1, i_y1+1, i_z1))  *(pos1_z - block(g)%zp(i_z1))/(block(g)%zp(i_z1+1) - block(g)%zp(i_z1))

         !This will be used for dpdx calculation at point 1
         p_x1 = p_z1_x1 + (p_z2_x1 - p_z1_x1)*(pos1_y - block(g)%yp(i_y1))/(block(g)%yp(i_y1+1)-block(g)%yp(i_y1))
         p_x2 = p_z1_x2 + (p_z2_x2 - p_z1_x2)*(pos1_y - block(g)%yp(i_y1))/(block(g)%yp(i_y1+1)-block(g)%yp(i_y1))

         p_pos1 = p_x1 + (p_x2 - p_x1)*(pos1_x - block(g)%xp(i_x1))/(block(g)%xp(i_x1+1)-block(g)%xp(i_x1))

         h2 = dabs(block(g)%xp(i_x1+1) - pos1_x)
         h1 = dabs(block(g)%xp(i_x1)   - pos1_x)
         dpdx_e = (h1**2*p_x2 - h2**2*p_x1 + (h2**2- h1**2)*p_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%yp(i_y1+1) - pos1_y)
         h1 = dabs(block(g)%yp(i_y1)   - pos1_y)
         dpdy_e = (h1**2*p_y2 - h2**2*p_y1 + (h2**2- h1**2)*p_pos1)/(h1*h2*(h1+h2)+1e-16)

         h2 = dabs(block(g)%zp(i_z1+1) - pos1_z)
         h1 = dabs(block(g)%zp(i_z1)   - pos1_z)
         dpdz_e = (h1**2*p_z2 - h2**2*p_z1 + (h2**2- h1**2)*p_pos1)/(h1*h2*(h1+h2)+1e-16)

         dpdn_e = dpdx_e*block(g)%cosAlpha(ielem) + dpdy_e*block(g)%cosBeta(ielem) + dpdz_e*block(g)%cosGamma(ielem)

         bval = dpdn  !dpdn=-dudt
         aval = (dpdn_e - dpdn)/(2*diagdis)
         cval = (p_pos1 - (dpdn_e + dpdn)*diagdis*.5)

         psurf = cval

         f_surf = cval*area

          f_surf_x = -1._rk*f_surf*block(g)%cosAlpha(ielem)*block(g)%rhof

          f_surf_y = -1._rk*f_surf*block(g)%cosBeta(ielem)*block(g)%rhof

	   f_surf_z = -1._rk* f_surf*block(g)%cosGamma(ielem)*block(g)%rhof

        !f_surf_x = -f_surf*block(g)%cosAlpha(ielem)

        !f_surf_y = -f_surf*block(g)%cosBeta(ielem)

	 !f_surf_z = -f_surf*block(g)%cosGamma(ielem)
!        write(12,*)'st',ielem,dudn_s,dvdn_s,dwdn_s,cval
!12        format(A2,I6,4F12.6)
!***********************drag calculation********************************
	  viscousDrag = viscousDrag + shear_x_force
         pressureDrag = pressureDrag + f_surf_x
         viscousLift = viscousLift + shear_y_force
         PressureLift = PressureLift + f_surf_y
         surf_area = surf_area + area
	  area_Sx = area_Sx + area_xz
         area_Sy = area_Sy + area_yz
!***********************************************************************
      ! block(g)%Total_V_Fx = block(g)%Total_V_Fx + shear_x_force
      ! block(g)%Total_P_Fx = block(g)%Total_P_Fx + f_surf_x

      ! block(g)%Total_V_Fy = block(g)%Total_V_Fy + shear_y_force
      ! block(g)%Total_P_Fy = block(g)%Total_P_Fy + f_surf_y
        total_v_x = total_v_x + shear_x_force
        total_p_x = total_p_x + f_surf_x

        total_v_y = total_v_y + shear_y_force
        total_p_y = total_p_y + f_surf_y
        !print*,ielem,block(g)%Total_V_Fy,block(g)%Total_P_Fy

	!END IF
	END DO
	!$acc end parallel
!***********************************************************************
          !!!$acc update host (stressElem, areaElem)
       ! !!$acc wait
        !print*,ielem,block(g)Total_V_Fy,block(g)%Total_P_Fy
        block(g)%Total_V_Fx=total_v_x
        block(g)%Total_V_Fy=total_v_y
        block(g)%Total_P_Fx=total_p_x
        block(g)%Total_P_Fy=total_p_y
        block(g)%Total_VP_FX = block(g)%Total_V_Fx + block(g)%Total_P_Fx
        block(g)%Total_VP_FY = block(g)%Total_V_Fy + block(g)%Total_P_Fy
        !area_Sx= 0.5 * area_Sx
        !area_Sy= 0.5 * area_Sy
        area_Sx= 1.  !0.5 * area_Sx
        area_Sy= 1.  !0.5 * area_Sy
        viscousDragcoefficient= 2* (viscousDrag/area_Sx)
        pressureDragcoefficient= 2* (pressureDrag/area_Sy)
        viscousLiftcoefficient= 2* (viscousLift/area_Sx)
        pressureLiftcoefficient= 2* (pressureLift/area_Sy)

        !*********************drag file writing*********************************

        WRITE(filename1,19)block(g)%dx
 19        FORMAT('dragcoff_',F6.4,'.dat')
       OPEN(899,file=filename1,Access='Append',status='unknown')
       WRITE(899,*) viscousDragcoefficient, pressureDragcoefficient, totime
       CLOSE(899)
        WRITE(filename1,29)block(g)%dx
 29       FORMAT('liftcoff_',F6.4,'.dat')
       OPEN(999,file=filename1,Access='Append',status='unknown')
       WRITE(999,*) viscousLiftcoefficient,  pressureLiftcoefficient, totime
       CLOSE(999)
        END DO

!*******************calculate stresses at node*******************
!      Anode = 0._rk
!      Afnode = 0._rk
!      stressNode = 0._rk
!
!        zlim1=19.2
!       ylim2=43.689
!        DO g=1,nblocks
!        block(g)%ibElemCntSt=0.
!        END DO
!       DO ielem = 1, block(g)%ibElems
!          p1x = block(g)%xnode1(block(g)%ibElP1(ielem))                       !x coordinate element node 1
!          p1y = block(g)%ynode1(block(g)%ibElP1(ielem))                       !y coordinate element node 1
!          p1z = block(g)%znode1(block(g)%ibElP1(ielem))                       !z coordinate element node 1
!
!          p2x = block(g)%xnode1(block(g)%ibElP2(ielem))                       !x coordinate element node 2
!          p2y = block(g)%ynode1(block(g)%ibElP2(ielem))                       !y coordinate element node 2
!          p2z = block(g)%znode1(block(g)%ibElP2(ielem))                       !z coordinate element node 2
!
!          p3x = block(g)%xnode1(block(g)%ibElP3(ielem))                       !x coordinate element node 3
!          p3y = block(g)%ynode1(block(g)%ibElP3(ielem))                       !y coordinate element node 3
!          p3z = block(g)%znode1(block(g)%ibElP3(ielem))                       !z coordinate element node 3
!
!
!          var_xcent =  (p2x+p1x+p3x)/3._rk                  !centroid x coordinate element
!          var_ycent =  (p2y+p1y+p3y)/3._rk                  !centroid y coordinate element
!          var_zcent =  (p2z+p1z+p3z)/3._rk                  !centroid z coordinate element

!         if ( var_ycent.le. ylim2 )then
!               if (var_zcent .le. zlim1)then
!                       g=1
!               else
!                       g=3
!               end if
!          else
!                       g=2
!          end if

!          block(g)%ibElemCntSt= block(g)%ibElemCntSt+1
!          nvib = block(g)%ibElemCntSt
!      DO k = 1, 3
!      Anode(block(g)%ibElP1(ielem), k) = Anode(block(g)%ibElP1(ielem), k) + block(g)%areaElem(nvib)
!      Afnode(block(g)%ibElP1(ielem), k) = Afnode(block(g)%ibElP1(ielem), k) + block(g)%stressElem(nvib, k)*block(g)%areaElem(nvib)
!
!      Anode(block(g)%ibElP2(ielem), k) = Anode(block(g)%ibElP2(ielem), k) + block(g)%areaElem(nvib)
!      Afnode(block(g)%ibElP2(ielem), k) = Afnode(block(g)%ibElP2(ielem), k) + block(g)%stressElem(nvib, k)*block(g)%areaElem(nvib)
!
!      Anode(block(g)%ibElP3(ielem), k) = Anode(block(g)%ibElP3(ielem), k) + block(g)%areaElem(nvib)
!      Afnode(block(g)%ibElP3(ielem), k) = Afnode(block(g)%ibElP3(ielem), k) + block(g)%stressElem(nvib, k)*block(g)%areaElem(nvib)
!      END DO

!       ENDDO
!
!
!      DO inode = 1, block(g)%ibNodes
!        DO k = 1, 3
!        stressNode(inode, k) = Afnode(inode, k)/(Anode(inode, k) + 1e-14)
!        END DO

!        modStressNode = dsqrt(stressNode(inode, 1)**2 +  stressNode(inode, 2)**2 + stressNode(inode, 3)**2 )
!        INSTWSS(inode) = modStressNode
!        SUMWSS(inode) = SUMWSS(inode) + modStressNode
!        SQSUMWSS(inode) = SQSUMWSS(inode) + modStressNode**2
!
!        DO k = 1, 3
!        SIGNWSS(inode, k) = SIGNWSS(inode, k) + stressNode(inode, k)
!        END DO
!      END DO
!

!      TAWSS = 0._rk
!      OSI = 0._rk
!      RRT = 0._rk
!wwwwwwwwwwwwwwwwwwwww Tec360 FE FILE WRITING wwwwwwwwwwwwwwwwwwwwwwwwwW
!      IF(ita.le.itamax .and.  mod(ita,500).eq.0) THEN

!       !print*,'in tawss'
!        DO inode = 1, block(g)%ibNodes
!       !print *,ita1
!          !TAWSS(inode) = SUMWSS(inode)/ita1
!          TAWSS(inode) = SUMWSS(inode)/ita2
!         ! modSIGNWSS = dsqrt(SIGNWSS(inode, 1)**2 + SIGNWSS(inode, 2)**2 + SIGNWSS(inode, 3)**2)
!         !OSI(inode) = 0.5*(1._rk-(modSIGNWSS/(SUMWSS(inode) + 1e-14)))

!         ! RRT(inode) = 1.0/(TAWSS(inode)*(1.0 - 2.0*OSI(inode)) + 1e-14)
!          WSSRMS(inode) = sqrt(SQSUMWSS(inode)/ita2)
!          !WSSRMS(inode) = sqrt(SQSUMWSS(inode)/ita1)
!        ! TavgStressNode1(inode) = SIGNWSS(inode,1)/ita1
!        ! TavgStressNode2(inode) = SIGNWSS(inode,2)/ita1
!        ! TavgStressNode3(inode) = SIGNWSS(inode,3)/ita1
!        END DO

!        WRITE(filename1,108) ita, re
!108     FORMAT('out/hv_aorta_stressdata.',i9.9,'.',f6.1,".dat")
!        OPEN(UNIT=857,FILE=filename1,STATUS='unknown')
!        WRITE(857,*) 'TITLE = "FEstressplot"'
!       ! WRITE(857,*) 'VARIABLES= "x", "y", "z", "St1", "St2", "St3", "TAWSS", "OSI", "RRT"'
!       ! WRITE(857,*) 'ZONE NODES= ',block(g)%ibNodes,',ELEMENTS= ',block(g)%ibElems,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
!        WRITE(857,*) 'VARIABLES= "x", "y", "z",  "TAWSS", "WSSRMS","INSTWSS"'
!        WRITE(857,*) 'ZONE NODES= ',block(g)%ibNodes,',ELEMENTS= ',block(g)%ibElems,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
!        DO inode = 1, block(g)%ibNodes
!          WRITE(857,*) block(g)%xnode1(inode),block(g)%ynode1(inode),block(g)%znode1(inode), TAWSS(inode),WSSRMS(inode), INSTWSS(inode)
!        END DO
!        DO ielem = 1, block(g)%ibElems
!          WRITE(857,*) block(g)%ibElP1(ielem), block(g)%ibElP2(ielem), block(g)%ibElP3(ielem)
!        END DO
!        CLOSE(857)
!      END IF
!wwwwwwwwwwwwwwwwwwwwwww END OF FE360 FILE WRITING wwwwwwwwwwwwwwwwwwwww

!*****************end of WSS, TAWSS, OSI calculation*******************

	END SUBROUTINE stressCal2


      SUBROUTINE writeStressData
       USE global
       IMPLICIT NONE
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER::  i, j, k, inode, g
        DO g=1,nblocks
       IF(mod(ita,5000)==0.OR.ita==itamax)THEN
       OPEN (1,FILE='stressdata',FORM='formatted')
        DO inode = 1, block(g)%ibNodes
         WRITE(1,*) SQSUMWSS(inode),SUMWSS(inode), SIGNWSS(inode,1), SIGNWSS(inode,2), SIGNWSS(inode,3), ita1, ita
        END DO
       CLOSE(1)
       END IF
        END DO
      END SUBROUTINE  writeStressData









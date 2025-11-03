module biocfd_stress_calculation
  use, intrinsic :: iso_fortran_env, only: dp => real64
  USE global, only: block, blk_start, mu_f, rho_f, totime
  IMPLICIT NONE

  private

  public :: stressCal1

  contains

       SUBROUTINE stressCal1

       INTEGER:: i, j, k, ielem, i_x1, i_y1, i_z1

       INTEGER:: i_cell, j_cell, k_cell,g

       REAL(dp):: diagdis, normdis, aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z

       REAL(dp):: xsurf, ysurf, zsurf, pos1_x, pos1_y, pos1_z,              &
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

       REAL(dp):: alen, area, area_xz, area_yz, area_xy

       REAL(dp):: shear_x_force, shear_y_force, shear_z_force,  &
                        f_surf, f_surf_x, f_surf_y, f_surf_z

       REAL(dp):: pressureDrag, viscousDrag, viscousLift,                        &
                        pressureLift, viscousDragcoefficient, pressureDragcoefficient, &
                        viscousLiftcoefficient, PressureLiftcoefficient, area_Sx,      &
                        area_Sy, surf_area

        real(dp) :: ac_y, ac_z, at_y, at_z

        CHARACTER(len=150) :: filename1

       DO g=blk_start, size(block)

       viscousDrag = 0.
       pressureDrag = 0.
       viscousLift = 0.
       PressureLift = 0.
       surf_area = 0.
       area_Sx = 0.
       area_Sy = 0
        !$acc parallel loop gang vector reduction(+: pressureDrag, viscousDrag, viscousLift, pressureLift, area_Sx, area_Sy, surf_area)  &
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
        !$acc          shear_z_force, f_surf, f_surf_x, f_surf_y, f_surf_z,ac_z,ac_y,ac_x,at_y,at_z)  &
        !$acc default(present)    &
        !$acc firstprivate (block(g)%nx, block(g)%ny, block(g)%nz, deltat, rho_f, re, mu_f)
        !DO ielem = 1, block(g)%ibElemCnt
        DO ielem = 1, block(g)%ibElems
     !  IF((block(g)%zcent(ielem).ge.0.0).and.(block(g)%zcent(ielem).le.2.0)) THEN
!***********************interpolation points****************************
       !$acc loop seq
       do i = 2, block(g)%nx+1
       if((block(g)%xcent(ielem)>=block(g)%x1(i)).and.(block(g)%xcent(ielem)<block(g)%x1(i+1)))then
       i_cell = i
       end if
       end do
       !$acc loop seq
       do j = 2, block(g)%ny+1
       if((block(g)%ycent(ielem)>=block(g)%y1(j)).and.(block(g)%ycent(ielem)<block(g)%y1(j+1)))then
       j_cell = j
       end if
       end do
       !$acc loop seq
       do k = 2, block(g)%nz+1
       if((block(g)%zcent(ielem)>=block(g)%z1(k)).and.(block(g)%zcent(ielem)<block(g)%z1(k+1)))then
       k_cell = k
       end if
       end do

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

!**************************velocity and pressure at the surface**********************
       IF (block(g)%ibSurfID(ielem)==50) THEN
        block(g)% thetaDot  =  0.
        block(g)% thetaDDot =  0.
             usurf = 0. + block(g)%xdot
             vsurf = 0. + block(g)%ydot
             wsurf = 0.
         ac_z      =  0.  !-thetaDot**2*(zcent(nelp(index_ts(n))) - piv_z)
         ac_y      =  0.  !-thetaDot**2*(ycent(nelp(index_ts(n))) - piv_y)
         at_z      =  0.  ! thetaDDot*(ycent(nelp(index_ts(n))) - piv_y)
         at_y      =  0.  !
       ELSEIF (block(g)%ibSurfID(ielem)==51) THEN
             block(g)% thetaDot  = block(g)% thetaDot1
            block(g)% thetaDDot = block(g)% thetaDDot1
             usurf    = 0. + block(g)%xdot
             vsurf    = -block(g)%thetaDot*(block(g)%zcent(ielem) - block(g)%piv_z) + block(g)%ydot
             wsurf    = block(g)%thetaDot*(block(g)%ycent(ielem) - block(g)%piv_y)
         ac_z      = -block(g)%thetaDot**2*(block(g)%zcent(ielem) -block(g)% piv_z)
         ac_y      = -block(g)%thetaDot**2*(block(g)%ycent(ielem) -block(g)% piv_y)
         at_z      = block(g)% thetaDDot*(block(g)%ycent(ielem)-block(g)% piv_y)
         at_y      = -block(g)%thetaDDot*(block(g)%zcent(ielem)-block(g)% piv_z)
       ELSEIF (block(g)%ibSurfId(ielem)==52) THEN
            block(g)% thetaDot  = block(g)% thetaDot2
            block(g)% thetaDDot = block(g)% thetaDDot2
             usurf = 0. +block(g)%xdot
             vsurf    = -block(g)%thetaDot*(block(g)%zcent(ielem) - block(g)%piv_z)+ block(g)%ydot  ! + ydot
             wsurf    = block(g)%thetaDot*(block(g)%ycent(ielem) - block(g)%piv_y)  ! + ydot
         ac_z      = -block(g)%thetaDot**2*(block(g)%zcent(ielem) -block(g)% piv_z)
         ac_y      = -block(g)%thetaDot**2*(block(g)%ycent(ielem) - block(g)%piv_y)
         at_z      =  block(g)%thetaDDot*(block(g)%ycent(ielem) -block(g)% piv_y)
         at_y      = -block(g)%thetaDDot*(block(g)%zcent(ielem) - block(g)%piv_z)
       ENDIF
!       usurf     =  0.
!       wsurf     = block(g)% thetaDot*(block(g)%ycent(ielem) -block(g)% piv_y)
!       vsurf     = -block(g)%thetaDot*(block(g)%zcent(ielem) -block(g)% piv_z)
!       dpdn      = -((ac_z + at_z)*block(g)%cosAlpha(ielem)  + (ac_y + at_y)*block(g)%cosBeta(ielem))
         dpdn = -((ac_z + at_z)*block(g)%cosGamma(ielem) + (ac_y + at_y)*block(g)%cosBeta(ielem))-block(g)%yddot*block(g)%cosBeta(ielem)

!*******************velocity interpolation at point 2******************

!******************u velocity interpolation at point 2******************
       !$acc loop seq
       DO i = 2, block(g)%nx+1
       if(pos1_x>=block(g)%xu(i).and.pos1_x<block(g)%xu(i+1)) i_x1 = i
       END DO
       !$acc loop seq
       DO j = 2, block(g)%ny+1
       if(pos1_y>=block(g)%yu(j).and.pos1_y<block(g)%yu(j+1)) i_y1 = j
       END DO
       !$acc loop seq
       DO k = 2, block(g)%nz+1
       if(pos1_z>=block(g)%zu(k).and.pos1_z<block(g)%zu(k+1)) i_z1 = k
       END DO

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
       !$acc loop seq
       DO i = 2, block(g)%nx+1
       if(pos1_x>=block(g)%xv(i).and.pos1_x<block(g)%xv(i+1)) i_x1 = i
       END DO
       !$acc loop seq
       DO j = 2, block(g)%ny+1
       if(pos1_y>=block(g)%yv(j).and.pos1_y<block(g)%yv(j+1)) i_y1 = j
       END DO
       !$acc loop seq
       DO k = 2, block(g)%nz+1
       if(pos1_z>=block(g)%zv(k).and.pos1_z<block(g)%zv(k+1)) i_z1 = k
       END DO

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
       !$acc loop seq
       DO i = 2, block(g)%nx+1
       if(pos1_x>=block(g)%xw(i).and.pos1_x<block(g)%xw(i+1)) i_x1 = i
       END DO
       !$acc loop seq
       DO j = 2, block(g)%ny+1
       if(pos1_y>=block(g)%yw(j).and.pos1_y<block(g)%yw(j+1)) i_y1 = j
       END DO
       !$acc loop seq
       DO k = 2, block(g)%nz+1
       if(pos1_z>=block(g)%zw(k).and.pos1_z<block(g)%zw(k+1)) i_z1 = k
       END DO

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
       ! alen = sqrt(block(g)%alpha3(ielem)**2.+block(g)%beta3(ielem)**2.+ block(g)%gamma3(ielem)**2.)
       alen = block(g)%element_length(ielem)
       area = alen/2.
       area_yz = 0.5*abs(alen * block(g)%cosAlpha(ielem))
       area_xz = 0.5*abs(alen * block(g)%cosBeta(ielem))
       area_xy = 0.5*abs(alen * block(g)%cosGamma(ielem))

!*********non-dimensional viscous stress & force calculation************
       stx1 = dudn_s - (dudn_s*block(g)%cosAlpha(ielem) + dvdn_s*block(g)%cosBeta(ielem) + dwdn_s*block(g)%cosGamma(ielem))*block(g)%cosAlpha(ielem)

       sty1 = dvdn_s - (dudn_s*block(g)%cosAlpha(ielem) + dvdn_s*block(g)%cosBeta(ielem) + dwdn_s*block(g)%cosGamma(ielem))*block(g)%cosBeta(ielem)

       stz1 = dwdn_s - (dudn_s*block(g)%cosAlpha(ielem) + dvdn_s*block(g)%cosBeta(ielem) + dwdn_s*block(g)%cosGamma(ielem))*block(g)%cosGamma(ielem)

       shear_x_force = mu_f*stx1*area
       shear_y_force = mu_f*sty1*area
       shear_z_force = mu_f*stz1*area

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
         cval = p_pos1 - (dpdn_e + dpdn)*diagdis*.5

         psurf = cval

         f_surf = cval*area

        f_surf_x = -f_surf*block(g)%cosAlpha(ielem)*rho_f

        f_surf_y = -f_surf*block(g)%cosBeta(ielem)*rho_f

       f_surf_z = -f_surf*block(g)%cosGamma(ielem)*rho_f

!***********************drag calculation********************************
         viscousDrag = viscousDrag + shear_x_force
         pressureDrag = pressureDrag + f_surf_x
         viscousLift = viscousLift + shear_y_force
         PressureLift = PressureLift + f_surf_y
         surf_area = surf_area + area
         area_Sx = area_Sx + area_xz
         area_Sy = area_Sy + area_yz

       END DO
        !$acc end parallel

        area_Sx= 0.5 * area_Sx
        area_Sy= 0.5 * area_Sy
       ! area_Sx= 1.!0.5 * area_Sx
       ! area_Sy= 1.!0.5 * area_Sy
        viscousDragcoefficient= 2* (viscousDrag/area_Sx)
        pressureDragcoefficient= 2* (pressureDrag/area_Sy)
        viscousLiftcoefficient= 2* (viscousLift/area_Sx)
        pressureLiftcoefficient= 2* (pressureLift/area_Sy)

        !*********************drag file writing*********************************

        WRITE(filename1,19)g
 19        FORMAT('dragcoff_',I4.4,'.dat')
       OPEN(899,file=filename1,Access='Append',status='unknown')
       WRITE(899,*) viscousDragcoefficient, pressureDragcoefficient, totime
       CLOSE(899)
        WRITE(filename1,29)g
 29       FORMAT('liftcoff_',I4.4,'.dat')
       OPEN(999,file=filename1,Access='Append',status='unknown')
       WRITE(999,*) viscousLiftcoefficient,  pressureLiftcoefficient, totime
       CLOSE(999)

        END DO

       END SUBROUTINE stressCal1
end module biocfd_stress_calculation

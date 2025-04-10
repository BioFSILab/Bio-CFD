module biocfd_forcing
  use global
  implicit none

  contains
SUBROUTINE pressureForcing1
      USE global
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1, g
      REAL (KIND = 8) :: n1, pos1_x, pos1_y, pos1_z, pt1, &
                         aval, bval, cval, p_pos1, sur2nodeDis, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, &
                         p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1, p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e

      dpdn = 0._rk
        !g=2
        DO g=blk_start,nblocks


 !$acc parallel loop gang vector                                                                                          &
 !$acc private (n,diagCell, n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc           p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,                &
 !$acc           p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc           i_z1,ac_z,ac_y,ac_x,at_y,at_z,ac_x_al,ac_y_al,at_x_al,at_y_al)         &
 !$acc default(present)  &
 !$acc firstprivate (block(g)%nx, block(g)%ny,block(g)%nz)

!!$omp parallel do private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, &
!!$omp cval, p_pos1, sur2nodeDis, dpdn,p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1, &
!!$omp  p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) firstprivate (g) num_threads(40)
      DO n = 1, block(g)%ibCellCount

         !dpdn = (-(v_curr-v_prev)/deltat)*block(g)%cosBeta(block(g)%nelp(n))
         IF (block(g)%ibSurfId(block(g)%nelp(n))==50) THEN
            !dpdn = 0.
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !
            ac_x_al = 0.
            ac_y_al = 0.
            at_x_al = 0.
            at_y_al = 0.
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==51) THEN
            block(g)%thetaDot  = block(g)%thetaDot1
            block(g)%thetaDDot = block(g)%thetaDDot1
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
          !!ac_x_al = -block(g)%alphaDot**2*(block(g)%xcent(block(g)%nelp(n)) - block(g)%piv_x )
          !!ac_y_al = -block(g)%alphaDot**2*(block(g)%ycent(block(g)%nelp(n)) -block(g)%piv_y)
          !!at_x_al =  block(g)%alphaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
          !!at_y_al = -block(g)%alphaDDot*(block(g)%xcent(block(g)%nelp(n)) - block(g)%piv_x)
            !dpdn = -((ac_z + at_z)*block(g)%cosAlpha(block(g)%nelp(n)) + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n)) !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==52) THEN
            block(g)%thetaDot  = block(g)%thetaDot2
            block(g)%thetaDDot = block(g)%thetaDDot2
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
          !!ac_x_al = -block(g)%alphaDot**2*(block(g)%xcent(block(g)%nelp(n)) - block(g)%piv_x )
          !!ac_y_al = -block(g)%alphaDot**2*(block(g)%ycent(block(g)%nelp(n)) -block(g)%piv_y)
          !!at_x_al =  block(g)%alphaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
          !!at_y_al = -block(g)%alphaDDot*(block(g)%xcent(block(g)%nelp(n)) - block(g)%piv_x)
            !dpdn = -((ac_z + at_z)*block(g)%cosAlpha(block(g)%nelp(n)) + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n)) !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))
        ENDIF
         dpdn = -((ac_z + at_z)*block(g)%cosGamma(block(g)%nelp(n)) + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n)))-block(g)%yddot*block(g)%cosBeta(block(g)%nelp(n))
         !dpdn = -((ac_z + at_z)*-block(g)%cosGamma(block(g)%nelp(n)) + (ac_y + at_y + ac_y_al + at_y_al)*-block(g)%cosBeta(block(g)%nelp(n)) + (ac_x_al+at_x_al)*-block(g)%cosAlpha(block(g)%nelp(n)))  !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))

         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

         sur2nodeDis = block(g)%pNormDis(n)

         pt1 = 1.5_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xp(i) + pt1*block(g)%cosAlpha(block(g)%nelp(n))
         pos1_y = block(g)%yp(j) + pt1*block(g)%cosBeta(block(g)%nelp(n))
         pos1_z = block(g)%zp(k) + pt1*block(g)%cosGamma(block(g)%nelp(n))

          !$acc loop seq
         DO il = 1, block(g)%nx+1
            if(pos1_x>=block(g)%xp(il).and.pos1_x<block(g)%xp(il+1)) i_x1 = il
         END DO
          !$acc loop seq
         DO jl = 1, block(g)%ny+1
            if(pos1_y>=block(g)%yp(jl).and.pos1_y<block(g)%yp(jl+1)) i_y1 = jl
         END DO
          !$acc loop seq
         DO kl = 1, block(g)%nz+1
            if(pos1_z>=block(g)%zp(kl).and.pos1_z<block(g)%zp(kl+1)) i_z1 = kl
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

         dpdn_e = dpdx_e*block(g)%cosAlpha(block(g)%nelp(n)) + dpdy_e*block(g)%cosBeta(block(g)%nelp(n)) + dpdz_e*block(g)%cosGamma(block(g)%nelp(n))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*.5

         block(g)%p(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
! !$omp end parallel do
 !$acc end parallel
      ENDDO

     !print*, "Leaving PressureForcing"
END SUBROUTINE pressureForcing1
!***********************************************************************

!***********************************************************************
SUBROUTINE velocityForcing1
      USE global
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1, g
      REAL (KIND = 8) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, h1, h2, &
                         usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2, &
                         vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2, &
                         wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2, &
                         dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e, &
                         dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1, &
                         u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
                         v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2

!!$acc parallel loop gang vector         &
!!$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1,           &
!!$acc          aval, bval, cval, sur2nodeDis, h1, h2,               &
!!$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,   &
!!$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,   &
!!$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,   &
!!$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e,  &
!!$acc          dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1,  &
!!$acc          u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
!!$acc          v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
!!$acc present (block(g)%interceptedIndexPtr, block(g)%deltax, block(g)%deltay, block(g)%deltaz, block(g)%cosAlpha, block(g)%cosBeta, block(g)%cosGamma,    &
!!$acc           ut, block(g)%u2NormDis, block(g)%u1NormDis, block(g)%nelu2, block(g)%nelu1, xu, yu, zu,      &
!!$acc           vt, block(g)%v2NormDis, block(g)%v1NormDis, block(g)%nelv2, block(g)%nelv1, xv, yv, zv,      &
!!$acc           wt, block(g)%w2NormDis, block(g)%w1NormDis, block(g)%nelw2, block(g)%nelw1, xw, yw, zw)      &
!!$acc firstprivate (block(g)%nx, block(g)%ny,block(g)%nz)


        !g=2
        DO g=blk_start, nblocks
 !$acc parallel loop gang vector         &
 !$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis, h1, h2,               &
 !$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,   &
 !$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,   &
 !$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,   &
 !$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e,  &
 !$acc          dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1,  &
 !$acc          u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
 !$acc          v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc default(present)   &
 !$acc firstprivate (block(g)%nx, block(g)%ny,block(g)%nz)


!!$omp parallel do private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1,           &
!!$omp         aval, bval, cval, sur2nodeDis, h1, h2,               &
!!$omp         usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,   &
!!$omp         vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,   &
!!$omp         wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,   &
!!!$omp         dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e,  &
!!!$omp         dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1,  &
!!!$omp         u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
!!!$omp         v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
!!!$omp         firstprivate (g)  num_threads(40)

      DO n = 1, block(g)%ibCellCount

         !usurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelu2(n))==50) THEN
             usurf = 0. + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu2(n))==51) THEN
             usurf = 0. + block(g)%xdot
           !!usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu2(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelu2(n))==52) THEN
             usurf = 0. +block(g)%xdot
           !!usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu2(n)) - block(g)%piv_y)
         ENDIF
         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

!***********************U(i,j,k)****************************************
         !usurf = u_curr
         sur2nodeDis = block(g)%u2NormDis(n)

         pt1 = 1.5_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i+1) + pt1*block(g)%cosAlpha(block(g)%nelu2(n))
         pos1_y = block(g)%yu(j)   + pt1*block(g)%cosBeta(block(g)%nelu2(n))
         pos1_z = block(g)%zu(k)   + pt1*block(g)%cosGamma(block(g)%nelu2(n))

         !$acc loop seq
         DO il = 1, block(g)%nx+2
            if(pos1_x>=block(g)%xu(il).and.pos1_x<block(g)%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, block(g)%ny+1
            if(pos1_y>=block(g)%yu(jl).and.pos1_y<block(g)%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, block(g)%nz+1
            if(pos1_z>=block(g)%zu(kl).and.pos1_z<block(g)%zu(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         IF(i_x1==block(g)%nx+2) i_x1 = block(g)%nx+1
         !IF(i_y1.EQ.block(g)%ny+2) i_y1 = block(g)%ny+1
         !IF(i_z1.EQ.block(g)%nz+2) i_z1 = block(g)%nz+1

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

         dudn_e = dudx_e*block(g)%cosAlpha(block(g)%nelu2(n)) + dudy_e*block(g)%cosBeta(block(g)%nelu2(n)) + dudz_e*block(g)%cosGamma(block(g)%nelu2(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%ut(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!******************************U(i-1,j,k)*******************************
         !usurf = u_curr
         !usurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelu1(n))==50) THEN
             usurf = 0. + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu1(n))==51) THEN
             usurf = 0. + block(g)%xdot
           !!usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu1(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelu1(n))==52) THEN
             usurf = 0. + block(g)%xdot
           !!usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu1(n)) - block(g)%piv_y)
         ENDIF
         sur2nodeDis = block(g)%u1NormDis(n)

         pt1 = 1.5_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i) + pt1*block(g)%cosAlpha(block(g)%nelu1(n))
         pos1_y = block(g)%yu(j) + pt1*block(g)%cosBeta(block(g)%nelu1(n))
         pos1_z = block(g)%zu(k) + pt1*block(g)%cosGamma(block(g)%nelu1(n))

         !$acc loop seq
         DO il = 1, block(g)%nx+1
            if(pos1_x>=block(g)%xu(il).and.pos1_x<block(g)%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, block(g)%ny+1
            if(pos1_y>=block(g)%yu(jl).and.pos1_y<block(g)%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, block(g)%nz+1
            if(pos1_z>=block(g)%zu(kl).and.pos1_z<block(g)%zu(kl+1)) i_z1 = kl
         END DO

         IF(i_x1==1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.block(g)%nx+2) i_x1 = block(g)%nx+1
         !IF(i_y1.EQ.block(g)%ny+2) i_y1 = block(g)%ny+1
         !IF(i_z1.EQ.block(g)%nz+2) i_z1 = block(g)%nz+1

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

         dudn_e = dudx_e*block(g)%cosAlpha(block(g)%nelu1(n)) + dudy_e*block(g)%cosBeta(block(g)%nelu1(n)) + dudz_e*block(g)%cosGamma(block(g)%nelu1(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%ut(i-1,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************V(i,j,k)*************************************
         !vsurf = v_curr
         !vsurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelv2(n))==50) THEN
                 vsurf = 0. + block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv2(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z) + block(g)%ydot
           !!vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z) + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(n))-block(g)%piv_x)
         ELSEIF (block(g)%ibSurfId(block(g)%nelv2(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z)+ block(g)%ydot  ! + ydot
           !!vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z) + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(n))-block(g)%piv_x)
         ENDIF

         sur2nodeDis = block(g)%v2NormDis(n)

         pt1 = 1.5_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) + pt1*block(g)%cosAlpha(block(g)%nelv2(n))
         pos1_y = block(g)%yv(j+1) + pt1*block(g)%cosBeta(block(g)%nelv2(n))
         pos1_z = block(g)%zv(k) + pt1*block(g)%cosGamma(block(g)%nelv2(n))

         !$acc loop seq
         DO il = 1, block(g)%nx+1
            if(pos1_x>=block(g)%xv(il).and.pos1_x<block(g)%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, block(g)%ny+2
            if(pos1_y>=block(g)%yv(jl).and.pos1_y<block(g)%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, block(g)%nz+1
            if(pos1_z>=block(g)%zv(kl).and.pos1_z<block(g)%zv(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.block(g)%nx+2) i_x1 = block(g)%nx+1
         IF(i_y1==block(g)%ny+2) i_y1 = block(g)%ny+1
         !IF(i_z1.EQ.block(g)%nz+2) i_z1 = block(g)%nz+1

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


         dvdn_e = dvdx_e*block(g)%cosAlpha(block(g)%nelv2(n)) + dvdy_e*block(g)%cosBeta(block(g)%nelv2(n)) +  dvdz_e*block(g)%cosGamma(block(g)%nelv2(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%vt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************V(i,j-1,k)*************************************
         !vsurf = v_curr
         !vsurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelv1(n))==50) THEN
           vsurf = 0.+ block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv1(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z) + block(g)%ydot
           !!vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)  + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(n))-block(g)%piv_x)
         ELSEIF (block(g)%ibSurfId(block(g)%nelv1(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)+ block(g)%ydot  ! + ydot
           !!vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)  + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(n))-block(g)%piv_x)
         ENDIF
         sur2nodeDis = block(g)%v1NormDis(n)

         pt1 = 1.5_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) + pt1*block(g)%cosAlpha(block(g)%nelv1(n))
         pos1_y = block(g)%yv(j) + pt1*block(g)%cosBeta(block(g)%nelv1(n))
         pos1_z = block(g)%zv(k) + pt1*block(g)%cosGamma(block(g)%nelv1(n))

         !$acc loop seq
         DO il = 1, block(g)%nx+1
            if(pos1_x>=block(g)%xv(il).and.pos1_x<block(g)%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, block(g)%ny+1
            if(pos1_y>=block(g)%yv(jl).and.pos1_y<block(g)%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, block(g)%nz+1
            if(pos1_z>=block(g)%zv(kl).and.pos1_z<block(g)%zv(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         IF(i_y1==1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.block(g)%nx+2) i_x1 = block(g)%nx+1
         !IF(i_y1.EQ.block(g)%ny+2) i_y1 = block(g)%ny+1
         !IF(i_z1.EQ.block(g)%nz+2) i_z1 = block(g)%nz+1

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


         dvdn_e = dvdx_e*block(g)%cosAlpha(block(g)%nelv1(n)) + dvdy_e*block(g)%cosBeta(block(g)%nelv1(n)) +  dvdz_e*block(g)%cosGamma(block(g)%nelv1(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%vt(i,j-1,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************W(i,j,k)*************************************
         !wsurf = w_curr
         !wsurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelw2(n))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw2(n))==51) THEN
           block(g)%thetaDot = block(g)%thetaDot1
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw2(n))==52) THEN
           block(g)%thetaDot = block(g)%thetaDot2
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(n)) - block(g)%piv_y)  ! + ydot
         ENDIF
         sur2nodeDis = block(g)%w2NormDis(n)

         pt1 = 1.5_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) + pt1*block(g)%cosAlpha(block(g)%nelw2(n))
         pos1_y = block(g)%yw(j) + pt1*block(g)%cosBeta(block(g)%nelw2(n))
         pos1_z = block(g)%zw(k+1) + pt1*block(g)%cosGamma(block(g)%nelw2(n))

         !$acc loop seq
         DO il = 1, block(g)%nx+1
            if(pos1_x>=block(g)%xw(il).and.pos1_x<block(g)%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, block(g)%ny+1
            if(pos1_y>=block(g)%yw(jl).and.pos1_y<block(g)%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, block(g)%nz+2
            if(pos1_z>=block(g)%zw(kl).and.pos1_z<block(g)%zw(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.block(g)%nx+2) i_x1 = block(g)%nx+1
         !IF(i_y1.EQ.block(g)%ny+2) i_y1 = block(g)%ny+1
         IF(i_z1==block(g)%nz+2) i_z1 = block(g)%nz+1

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


         dwdn_e = dwdx_e*block(g)%cosAlpha(block(g)%nelw2(n)) + dwdy_e*block(g)%cosBeta(block(g)%nelw2(n)) +  dwdz_e*block(g)%cosGamma(block(g)%nelw2(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%wt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************W(i,j,k-1)*************************************
         !wsurf = w_curr
         wsurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelw1(n))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw1(n))==51) THEN
           block(g)%thetaDot = block(g)%thetaDot1
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw1(n))==52) THEN
           block(g)%thetaDot = block(g)%thetaDot2
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(n)) - block(g)%piv_y)  ! + ydot
         ENDIF
         sur2nodeDis = block(g)%w1NormDis(n)

         pt1 = 1.5_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) + pt1*block(g)%cosAlpha(block(g)%nelw1(n))
         pos1_y = block(g)%yw(j) + pt1*block(g)%cosBeta(block(g)%nelw1(n))
         pos1_z = block(g)%zw(k) + pt1*block(g)%cosGamma(block(g)%nelw1(n))

         !$acc loop seq
         DO il = 1, block(g)%nx+1
            if(pos1_x>=block(g)%xw(il).and.pos1_x<block(g)%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, block(g)%ny+1
            if(pos1_y>=block(g)%yw(jl).and.pos1_y<block(g)%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, block(g)%nz+1
            if(pos1_z>=block(g)%zw(kl).and.pos1_z<block(g)%zw(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         IF(i_z1==1) i_z1 = 2
         !IF(i_x1.EQ.block(g)%nx+2) i_x1 = block(g)%nx+1
         !IF(i_y1.EQ.block(g)%ny+2) i_y1 = block(g)%ny+1
         !IF(i_z1.EQ.block(g)%nz+2) i_z1 = block(g)%nz+1

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


         dwdn_e = dwdx_e*block(g)%cosAlpha(block(g)%nelw1(n)) + dwdy_e*block(g)%cosBeta(block(g)%nelw1(n)) +  dwdz_e*block(g)%cosGamma(block(g)%nelw1(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%wt(i,j,k-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************
      ENDDO
!!!$omp end parallel do
!$acc end parallel
      ENDDO

     !print*, "Leaving VelocityForcing"
END SUBROUTINE velocityForcing1



!***********************************************************************
SUBROUTINE pressureForcingGhost
      USE global
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: g,n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (KIND = 8) :: n1, pos1_x, pos1_y, pos1_z, pt1, &
                         aval, bval, cval, p_pos1, sur2nodeDis, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, &
                         p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1, p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e



        !g=2
        DO g=blk_start, nblocks
      dpdn = 0._rk
        !print*,'inside pressureghost'
      !!$acc parallel loop present(block(g)%ibSurfId, block(g)%xcent, block(g)%ycent, block(g)%zcent, cell2, block(g)%TSIndexPtr, pNormDis, nelp, deltax, deltay, deltaz, block(g)%cosAlpha, block(g)%cosBeta, block(g)%cosGamma, xp, yp, zp, p, block(g)%p_ghost, block(g)%pt_ghost, cell, block(g)%index_ts)
 !$acc parallel loop gang vector                                                                    &
 !$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc           p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,                &
 !$acc           p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc           i_z1,ac_z,ac_y,ac_x,at_y,at_z,ac_x_al,ac_y_al,at_x_al,at_y_al)         &
 !$acc default(present)  &
 !$acc firstprivate (block(g)%nx, block(g)%ny,block(g)%nz)
      DO n = 1, block(g)%TSCellCount

         !dpdn = (-(v_curr-v_prev)/deltat)*block(g)%cosBeta(block(g)%nelp(n))
        !print *,n
        i = block(g)%TSIndexPtr(n, 1)
        j = block(g)%TSIndexPtr(n, 2)
        k = block(g)%TSIndexPtr(n, 3)
        !PRINT*,I,J,K
       ! IF (block(g)%cell2(i,j,k).EQ.2) THEN
        IF (block(g)%ibSurfId(block(g)%nelp(block(g)%index_ts(n)))==50) THEN
            !dpdn = 0.
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !-block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_x_al = 0.
            ac_y_al = 0.
            at_x_al = 0.
            at_y_al = 0.
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(block(g)%index_ts(n)))==51) THEN
            block(g)%thetaDot  = block(g)%thetaDot1
            block(g)%thetaDDot = block(g)%thetaDDot1
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)

          !!ac_x_al = -block(g)%alphaDot**2*(block(g)%xcent(block(g)%index_ts(n)) - block(g)%piv_x )
          !!ac_y_al = -block(g)%alphaDot**2*(block(g)%ycent(block(g)%index_ts(n)) -block(g)%piv_y)
          !!at_x_al =  block(g)%alphaDDot*(block(g)%ycent(block(g)%index_ts(n)) - block(g)%piv_y)
          !!at_y_al = -block(g)%alphaDDot*(block(g)%xcent(block(g)%index_ts(n)) - block(g)%piv_x)

            !dpdn = -((ac_z + at_z)*-block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(n))) + (ac_y + at_y)*-block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(n))))  !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(block(g)%index_ts(n)))==52) THEN
            block(g)%thetaDot  = block(g)%thetaDot2
            block(g)%thetaDDot = block(g)%thetaDDot2
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n)))  - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
          !!ac_x_al = -block(g)%alphaDot**2*(block(g)%xcent(block(g)%index_ts(n)) - block(g)%piv_x )
          !!ac_y_al = -block(g)%alphaDot**2*(block(g)%ycent(block(g)%index_ts(n)) -block(g)%piv_y)
          !!at_x_al =  block(g)%alphaDDot*(block(g)%ycent(block(g)%index_ts(n)) - block(g)%piv_y)
          !!at_y_al = -block(g)%alphaDDot*(block(g)%xcent(block(g)%index_ts(n)) - block(g)%piv_x)
        ENDIF
            dpdn = -((ac_z + at_z)* -block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(n))) + (ac_y + at_y)* -block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(n)))) -(block(g)%yddot*(-block(g)%cosBeta(block(g)%index_ts(n))))
         !dpdn = -((ac_z + at_z)*-block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(n))) + (ac_y + at_y + ac_y_al + at_y_al)*-block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(n))) + (ac_x_al+at_x_al)*-block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(n))))  !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))

         sur2nodeDis = -block(g)%pNormDis(block(g)%index_ts(n))

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xp(i) + pt1*-block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(n)))
         pos1_y = block(g)%yp(j) + pt1*-block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(n)))
         pos1_z = block(g)%zp(k) + pt1*-block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(n)))
         !ac_x = -block(g)%thetaDot**2*(block(g)%xcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_x)
         !ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
         !at_x =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
         !at_y = -block(g)%thetaDDot*(block(g)%xcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_x)
         !dpdn = -((ac_x + at_x)*-block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(n))) + (ac_y + at_y)*-block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(n))) + yddot*-block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(n))))
         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xp(il).and.pos1_x<block(g)%xp(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yp(jl).and.pos1_y<block(g)%yp(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zp(kl).and.pos1_z<block(g)%zp(kl+1)) i_z1 = kl
         END DO
         !IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
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

         dpdn_e = dpdx_e*-block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(n))) + dpdy_e*-block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(n))) + dpdz_e*-block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(n)))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*.5

         block(g)%p_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%p_ghost(n) =  p(i,j,k)
         !ENDIF
         block(g)%pt_ghost(n) = block(g)%p(i,j,k)
         !ENDIF
      ENDDO
      !$acc end parallel
      ENDDO
END SUBROUTINE pressureForcingGhost


!***********************************************************************
SUBROUTINE velocityForcingGhost
      USE global
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: g,n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (KIND = 8) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, h1, h2, &
                         usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2, &
                         vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2, &
                         wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2, &
                         dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e, &
                         dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1, &
                         u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
                         v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2

        !g=2
        !g=2
        DO g=blk_start,nblocks
      !!$acc parallel loop present(block(g)%ibSurfId, block(g)%index_ts, block(g)%xcent, block(g)%ycent, block(g)%zcent, u, v, w, block(g)%TSIndexPtr, cell2, block(g)%u2_ghost, block(g)%u2t_ghost, block(g)%v2_ghost, block(g)%v2t_ghost, block(g)%w2_ghost, block(g)%w2t_ghost, block(g)%u1_ghost, block(g)%u1t_ghost, block(g)%v1_ghost, block(g)%v1t_ghost, block(g)%w1_ghost, block(g)%w1t_ghost, block(g)%u2NormDis, block(g)%u1NormDis, block(g)%v2NormDis, block(g)%v1NormDis, block(g)%w2NormDis, block(g)%w1NormDis, block(g)%nelu2, block(g)%nelu1, block(g)%nelv2, block(g)%nelv1, block(g)%nelw2, block(g)%nelw1, deltax, deltay, deltaz, block(g)%cosAlpha, block(g)%cosBeta, cosgamma, xu, yu, zu, xv, yv, zv, xw, yw, zw, ut, vt, wt, cell)
 !$acc parallel loop gang vector         &
 !$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis, h1, h2,               &
 !$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,   &
 !$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,   &
 !$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,   &
 !$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e,  &
 !$acc          dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1,  &
 !$acc          u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
 !$acc          v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc default(present)   &
 !$acc firstprivate (block(g)%nx, block(g)%ny,block(g)%nz)
      DO n = 1, block(g)%TSCellCount

        i = block(g)%TSIndexPtr(n, 1)
        j = block(g)%TSIndexPtr(n, 2)
        k = block(g)%TSIndexPtr(n, 3)
        !IF (block(g)%cell2(i,j,k).EQ.2) THEN
!***********************U(i,j,k)****************************************
         !usurf = u_curr
         IF (block(g)%ibSurfID(block(g)%nelu2(block(g)%index_ts(n)))==50) THEN
             usurf = 0. + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu2(block(g)%index_ts(n)))==51) THEN
             usurf = 0. + block(g)%xdot
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu2(block(g)%index_ts(n))) -block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelu2(block(g)%index_ts(n)))==52) THEN
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu2(block(g)%index_ts(n)))-block(g)%piv_y)
           usurf = 0. + block(g)%xdot
         ENDIF

         !usurf = block(g)%thetaDot*(block(g)%ycent(block(g)%nelu2(block(g)%index_ts(n)))-block(g)%piv_y)
         sur2nodeDis = -block(g)%u2NormDis(block(g)%index_ts(n))

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i+1) + pt1*-block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(n)))
         pos1_y = block(g)%yu(j)   + pt1*-block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(n)))
         pos1_z = block(g)%zu(k)   + pt1*-block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xu(il).and.pos1_x<block(g)%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yu(jl).and.pos1_y<block(g)%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zu(kl).and.pos1_z<block(g)%zu(kl+1)) i_z1 = kl
         END DO

         !IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
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

         dudn_e = dudx_e*-block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(n))) + dudy_e*-block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(n))) + dudz_e*-block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(n)))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%u2_ghost(n) = u(i,j,k)
         !ENDIF
         block(g)%u2t_ghost(n) = block(g)%u(i,j,k)
!***********************************************************************

!******************************U(i-1,j,k)*******************************
         !usurf = u_curr
         IF (block(g)%ibSurfID(block(g)%nelu1(block(g)%index_ts(n)))==50) THEN
             usurf = 0. + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu1(block(g)%index_ts(n)))==51) THEN
             usurf = 0. + block(g)%xdot
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu1(block(g)%index_ts(n)))-block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelu1(block(g)%index_ts(n)))==52) THEN
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu1(block(g)%index_ts(n)))-block(g)%piv_y)
           usurf = 0. + block(g)%xdot
         ENDIF


         !usurf = block(g)%thetaDot*(block(g)%ycent(block(g)%nelu1(block(g)%index_ts(n)))-block(g)%piv_y)
         sur2nodeDis = -block(g)%u1NormDis(block(g)%index_ts(n))

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i) + pt1*-block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(n)))
         pos1_y = block(g)%yu(j) + pt1*-block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(n)))
         pos1_z = block(g)%zu(k) + pt1*-block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xu(il).and.pos1_x<block(g)%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yu(jl).and.pos1_y<block(g)%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zu(kl).and.pos1_z<block(g)%zu(kl+1)) i_z1 = kl
         END DO

         !IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
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

         dudn_e = dudx_e*-block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(n))) + dudy_e*-block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(n))) + dudz_e*-block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(n)))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%u1_ghost(n) = u(i-1,j,k)
         !ENDIF
         block(g)%u1t_ghost(n) = block(g)%u(i-1,j,k)
!***********************************************************************

!**************************V(i,j,k)*************************************
         !vsurf = v_curr
         IF (block(g)%ibSurfID(block(g)%nelv2(block(g)%index_ts(n)))==50) THEN
           vsurf = 0.+ block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv2(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(block(g)%index_ts(n))) - block(g)%piv_z) + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(block(g)%index_ts(n))) - block(g)%piv_x)+ block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv2(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(block(g)%index_ts(n))) - block(g)%piv_z) + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(block(g)%index_ts(n))) - block(g)%piv_x) + block(g)%ydot  ! + ydot
         ENDIF

         !vsurf =  -block(g)%thetaDot*(block(g)%xcent(block(g)%nelv2(block(g)%index_ts(n))) - block(g)%piv_x) + ydot
         sur2nodeDis = -block(g)%v2NormDis(block(g)%index_ts(n))

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) + pt1*-block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(n)))
         pos1_y = block(g)%yv(j+1) + pt1*-block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(n)))
         pos1_z = block(g)%zv(k) + pt1*-block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xv(il).and.pos1_x<block(g)%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yv(jl).and.pos1_y<block(g)%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zv(kl).and.pos1_z<block(g)%zv(kl+1)) i_z1 = kl
         END DO

         ! IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
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


         dvdn_e = dvdx_e*-block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(n))) + dvdy_e*-block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(n))) +  dvdz_e*-block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(n)))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%v2_ghost(n) = v(i,j,k)
         !ENDIF
         block(g)%v2t_ghost(n) = block(g)%v(i,j,k)
!***********************************************************************

!**************************V(i,j-1,k)*************************************
         !vsurf = v_curr
         IF (block(g)%ibSurfID(block(g)%nelv1(block(g)%index_ts(n)))==50) THEN
           vsurf = 0.+ block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv1(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(block(g)%index_ts(n))) - block(g)%piv_z)  + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(block(g)%index_ts(n))) - block(g)%piv_x) + block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv1(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(block(g)%index_ts(n))) - block(g)%piv_z)  + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(block(g)%index_ts(n))) - block(g)%piv_x)+ block(g)%ydot  ! + ydot
         ENDIF


         !vsurf =  -block(g)%thetaDot*(block(g)%xcent(block(g)%nelv1(block(g)%index_ts(n))) - block(g)%piv_x) + ydot
         sur2nodeDis = -block(g)%v1NormDis(block(g)%index_ts(n))

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) + pt1*-block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(n)))
         pos1_y = block(g)%yv(j) + pt1*-block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(n)))
         pos1_z = block(g)%zv(k) + pt1*-block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xv(il).and.pos1_x<block(g)%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yv(jl).and.pos1_y<block(g)%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zv(kl).and.pos1_z<block(g)%zv(kl+1)) i_z1 = kl
         END DO

         !IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
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


         dvdn_e = dvdx_e*-block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(n))) + dvdy_e*-block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(n))) +  dvdz_e*-block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(n)))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%v1_ghost(n) = v(i,j-1,k)
         !ENDIF
         block(g)%v1t_ghost(n) = block(g)%v(i,j-1,k)
!***********************************************************************

!**************************W(i,j,k)*************************************
         !wsurf = w_curr
         wsurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelw2(block(g)%index_ts(n)))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw2(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           wsurf    =  block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(block(g)%index_ts(n))) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw2(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           wsurf    =  block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(block(g)%index_ts(n))) - block(g)%piv_y)  ! + ydot
         ENDIF


         sur2nodeDis = -block(g)%w2NormDis(block(g)%index_ts(n))

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) + pt1*-block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(n)))
         pos1_y = block(g)%yw(j) + pt1*-block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(n)))
         pos1_z = block(g)%zw(k+1) + pt1*-block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xw(il).and.pos1_x<block(g)%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yw(jl).and.pos1_y<block(g)%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zw(kl).and.pos1_z<block(g)%zw(kl+1)) i_z1 = kl
         END DO
         !IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
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


         dwdn_e = dwdx_e*-block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(n))) + dwdy_e*-block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(n))) +  dwdz_e*-block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(n)))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%w2_ghost(n) = w(i,j,k)
         !ENDIF
         block(g)%w2t_ghost(n) = block(g)%w(i,j,k)

!***********************************************************************

!**************************W(i,j,k-1)*************************************
         !wsurf = w_curr
         wsurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelw1(block(g)%index_ts(n)))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw1(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           wsurf    =  block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(block(g)%index_ts(n))) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw1(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           wsurf    =  block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(block(g)%index_ts(n))) - block(g)%piv_y)  ! + ydot
         ENDIF


         sur2nodeDis = -block(g)%w1NormDis(block(g)%index_ts(n))

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) + pt1*-block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(n)))
         pos1_y = block(g)%yw(j) + pt1*-block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(n)))
         pos1_z = block(g)%zw(k) + pt1*-block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xw(il).and.pos1_x<block(g)%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yw(jl).and.pos1_y<block(g)%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zw(kl).and.pos1_z<block(g)%zw(kl+1)) i_z1 = kl
         END DO
         !IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
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


         dwdn_e = dwdx_e*-block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(n))) + dwdy_e*-block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(n))) +  dwdz_e*-block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(n)))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%w1_ghost(n) = w(i,j,k-1)
         !ENDIF
         block(g)%w1t_ghost(n) = block(g)%w(i,j,k-1)
!***********************************************************************
        !ENDIF
      ENDDO
      !$acc end parallel
      ENDDO

END SUBROUTINE velocityForcingGhost
!***********************************************************************

SUBROUTINE pressureForcingField
      USE global
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: g,n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (KIND = 8) :: n1, pos1_x, pos1_y, pos1_z, pt1, &
                         aval, bval, cval, p_pos1, sur2nodeDis, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, &
                         p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1, p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e

        !g=2
       DO g=blk_start,nblocks
      dpdn = 0._rk
      !!$acc parallel loop present(block(g)%ibSurfId, block(g)%xcent, block(g)%ycent, block(g)%zcent, block(g)%interceptedIndexPtr, block(g)%pNormDis, nelp, deltax, deltay, deltaz, block(g)%cosAlpha, block(g)%cosBeta, block(g)%cosGamma, xp, yp, zp, p)
 !$acc parallel loop gang vector                                                                                          &
 !$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc           p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,                &
 !$acc           p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc           i_z1,ac_z,ac_y,ac_x,at_y,at_z,ac_x_al,ac_y_al,at_x_al,at_y_al)         &
 !$acc default(present)  &
 !$acc firstprivate (block(g)%nx, block(g)%ny,block(g)%nz)
      DO n = 1, block(g)%ibCellCount
         IF (block(g)%ibSurfId(block(g)%nelp(n))==50) THEN
            !dpdn = 0.
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !
            ac_x_al = 0.
            ac_y_al = 0.
            at_x_al = 0.
            at_y_al = 0.
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==51) THEN
            block(g)%thetaDot  = block(g)%thetaDot1
            block(g)%thetaDDot = block(g)%thetaDDot1
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
          !!ac_x_al = -block(g)%alphaDot**2*(block(g)%xcent(block(g)%nelp(n)) -block(g)%piv_x )
          !!ac_y_al = -block(g)%alphaDot**2*(block(g)%ycent(block(g)%nelp(n)) -block(g)%piv_y )
          !!at_x_al =  block(g)%alphaDDot*(block(g)%ycent(block(g)%nelp(n)) -block(g)%piv_y )
          !!at_y_al = -block(g)%alphaDDot*(block(g)%xcent(block(g)%nelp(n)) -block(g)%piv_x )
            !dpdn = -((ac_z + at_z)*block(g)%cosAlpha(block(g)%nelp(n)) + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n)) !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==52) THEN
            block(g)%thetaDot  = block(g)%thetaDot2
            block(g)%thetaDDot = block(g)%thetaDDot2
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
          !!ac_x_al = -block(g)%alphaDot**2*(block(g)%xcent(block(g)%nelp(n)) -block(g)%piv_x )
          !!ac_y_al = -block(g)%alphaDot**2*(block(g)%ycent(block(g)%nelp(n)) -block(g)%piv_y )
          !!at_x_al =  block(g)%alphaDDot*(block(g)%ycent(block(g)%nelp(n)) -block(g)%piv_y )
          !!at_y_al = -block(g)%alphaDDot*(block(g)%xcent(block(g)%nelp(n)) -block(g)%piv_x )
            !dpdn = -((ac_z + at_z)*block(g)%cosAlpha(block(g)%nelp(n)) + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n)) !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))
        ENDIF
         dpdn = -((ac_z + at_z)*block(g)%cosGamma(block(g)%nelp(n)) + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n)))-block(g)%yddot*block(g)%cosBeta(block(g)%nelp(n))
         !dpdn = -((ac_z + at_z)*-block(g)%cosGamma(block(g)%nelp(n)) + (ac_y + at_y + ac_y_al + at_y_al)*-block(g)%cosBeta(block(g)%nelp(n)) + (ac_x_al+at_x_al)*-block(g)%cosAlpha(block(g)%nelp(n)))  !+ yddot*block(g)%cosBeta(block(g)%nelp(n)))

         !ac_x = -block(g)%thetaDot**2*(block(g)%xcent(block(g)%nelp(n)) - block(g)%piv_x)
         !ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
         !at_x =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
         !at_y = -block(g)%thetaDDot*(block(g)%xcent(block(g)%nelp(n)) - block(g)%piv_x)
         !dpdn = -((ac_x + at_x)*block(g)%cosAlpha(block(g)%nelp(n)) + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n)) + yddot*block(g)%cosBeta(block(g)%nelp(n)))

         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

         sur2nodeDis = block(g)%pNormDis(n)

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xp(i) + pt1*block(g)%cosAlpha(block(g)%nelp(n))
         pos1_y = block(g)%yp(j) + pt1*block(g)%cosBeta(block(g)%nelp(n))
         pos1_z = block(g)%zp(k) + pt1*block(g)%cosGamma(block(g)%nelp(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xp(il).and.pos1_x<block(g)%xp(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yp(jl).and.pos1_y<block(g)%yp(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zp(kl).and.pos1_z<block(g)%zp(kl+1)) i_z1 = kl
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

         dpdn_e = dpdx_e*block(g)%cosAlpha(block(g)%nelp(n)) + dpdy_e*block(g)%cosBeta(block(g)%nelp(n)) + dpdz_e*block(g)%cosGamma(block(g)%nelp(n))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*.5

         block(g)%p(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
      !$acc end parallel
      ENDDO
END SUBROUTINE pressureForcingField

!*****************************************************
!*****************************************************************************
SUBROUTINE velocityForcingField
      USE global
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: g,n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (KIND = 8) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, h1, h2, &
                         usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2, &
                         vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2, &
                         wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2, &
                         dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e, &
                         dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1, &
                         u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
                         v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2

        !g=2
        DO g=blk_start,nblocks
      !!$acc parallel loop present(block(g)%ibSurfId, block(g)%xcent, block(g)%ycent, block(g)%zcent, u, v, w, block(g)%interceptedIndexPtr, block(g)%u2NormDis, block(g)%u1NormDis, block(g)%v2NormDis, block(g)%v1NormDis, block(g)%w2NormDis, block(g)%w1NormDis, block(g)%nelu2, block(g)%nelu1, block(g)%nelv2, block(g)%nelv1, block(g)%nelw2, block(g)%nelw1, deltax, block(g)%deltay, deltaz, block(g)%cosAlpha, block(g)%cosBeta, cosgamma, xu, yu, zu, xv, yv, zv, xw, yw, zw, ut, vt, wt)
 !$acc parallel loop gang vector         &
 !$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis, h1, h2,               &
 !$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,   &
 !$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,   &
 !$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,   &
 !$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e,  &
 !$acc          dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1,  &
 !$acc          u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
 !$acc          v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc default(present)   &
 !$acc firstprivate (block(g)%nx, block(g)%ny,block(g)%nz)
      DO n = 1, block(g)%ibCellCount

         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

!***********************U(i,j,k)****************************************
         !usurf = u_curr
         IF (block(g)%ibSurfID(block(g)%nelu2(n))==50) THEN
             usurf = 0.+ block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu2(n))==51) THEN
             usurf = 0. + block(g)%xdot
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu2(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelu2(n))==52) THEN
             usurf = 0. + block(g)%xdot
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu2(n)) - block(g)%piv_y)
         ENDIF

         !usurf = block(g)%thetaDot*(block(g)%ycent(block(g)%nelu2(n)) - block(g)%piv_y)
         sur2nodeDis = block(g)%u2NormDis(n)

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i+1) + pt1*block(g)%cosAlpha(block(g)%nelu2(n))
         pos1_y = block(g)%yu(j)   + pt1*block(g)%cosBeta(block(g)%nelu2(n))
         pos1_z = block(g)%zu(k)   + pt1*block(g)%cosGamma(block(g)%nelu2(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xu(il).and.pos1_x<block(g)%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yu(jl).and.pos1_y<block(g)%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zu(kl).and.pos1_z<block(g)%zu(kl+1)) i_z1 = kl
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

         dudn_e = dudx_e*block(g)%cosAlpha(block(g)%nelu2(n)) + dudy_e*block(g)%cosBeta(block(g)%nelu2(n)) + dudz_e*block(g)%cosGamma(block(g)%nelu2(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!******************************U(i-1,j,k)*******************************
         !usurf = u_curr
         IF (block(g)%ibSurfID(block(g)%nelu1(n))==50) THEN
             usurf = 0. + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu1(n))==51) THEN
             usurf = 0. + block(g)%xdot
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu1(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelu1(n))==52) THEN
             usurf = 0. + block(g)%xdot
           !usurf = block(g)%alphaDot * (block(g)%ycent(block(g)%nelu1(n)) - block(g)%piv_y)
         ENDIF

         !usurf = block(g)%thetaDot*(block(g)%ycent(block(g)%nelu1(n)) - block(g)%piv_y)
         sur2nodeDis = block(g)%u1NormDis(n)

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i) + pt1*block(g)%cosAlpha(block(g)%nelu1(n))
         pos1_y = block(g)%yu(j) + pt1*block(g)%cosBeta(block(g)%nelu1(n))
         pos1_z = block(g)%zu(k) + pt1*block(g)%cosGamma(block(g)%nelu1(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xu(il).and.pos1_x<block(g)%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yu(jl).and.pos1_y<block(g)%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zu(kl).and.pos1_z<block(g)%zu(kl+1)) i_z1 = kl
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

         dudn_e = dudx_e*block(g)%cosAlpha(block(g)%nelu1(n)) + dudy_e*block(g)%cosBeta(block(g)%nelu1(n)) + dudz_e*block(g)%cosGamma(block(g)%nelu1(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u(i-1,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************V(i,j,k)*************************************
         !vsurf = v_curr
         IF (block(g)%ibSurfID(block(g)%nelv2(n))==50) THEN
           vsurf = 0.+ block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv2(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z)
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z) + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(n)) - block(g)%piv_x)+ block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv2(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z)! + ydot
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z) + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(n)) - block(g)%piv_x)+ block(g)%ydot
         ENDIF

         !vsurf =  -block(g)%thetaDot*(block(g)%xcent(block(g)%nelv2(n)) - block(g)%piv_x) + ydot
         sur2nodeDis = block(g)%v2NormDis(n)

         pt1 = 1.21_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) + pt1*block(g)%cosAlpha(block(g)%nelv2(n))
         pos1_y = block(g)%yv(j+1) + pt1*block(g)%cosBeta(block(g)%nelv2(n))
         pos1_z = block(g)%zv(k) + pt1*block(g)%cosGamma(block(g)%nelv2(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xv(il).and.pos1_x<block(g)%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yv(jl).and.pos1_y<block(g)%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zv(kl).and.pos1_z<block(g)%zv(kl+1)) i_z1 = kl
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


         dvdn_e = dvdx_e*block(g)%cosAlpha(block(g)%nelv2(n)) + dvdy_e*block(g)%cosBeta(block(g)%nelv2(n)) +  dvdz_e*block(g)%cosGamma(block(g)%nelv2(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************V(i,j-1,k)*************************************
         !vsurf = v_curr
         IF (block(g)%ibSurfID(block(g)%nelv1(n))==50) THEN
           vsurf = 0.+ block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv1(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)  + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(n)) - block(g)%piv_x) + block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv1(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)! + ydot
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)  + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(n)) - block(g)%piv_x) + block(g)%ydot
         ENDIF

         !vsurf =  -block(g)%thetaDot*(block(g)%xcent(block(g)%nelv1(n)) - block(g)%piv_x) + ydot
         sur2nodeDis = block(g)%v1NormDis(n)

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) + pt1*block(g)%cosAlpha(block(g)%nelv1(n))
         pos1_y = block(g)%yv(j) + pt1*block(g)%cosBeta(block(g)%nelv1(n))
         pos1_z = block(g)%zv(k) + pt1*block(g)%cosGamma(block(g)%nelv1(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xv(il).and.pos1_x<block(g)%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yv(jl).and.pos1_y<block(g)%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zv(kl).and.pos1_z<block(g)%zv(kl+1)) i_z1 = kl
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


         dvdn_e = dvdx_e*block(g)%cosAlpha(block(g)%nelv1(n)) + dvdy_e*block(g)%cosBeta(block(g)%nelv1(n)) +  dvdz_e*block(g)%cosGamma(block(g)%nelv1(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v(i,j-1,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************W(i,j,k)*************************************
         !wsurf = w_curr
         IF (block(g)%ibSurfID(block(g)%nelw2(n))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw2(n))==51) THEN
           block(g)%thetaDot = block(g)%thetaDot1
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw2(n))==52) THEN
           block(g)%thetaDot = block(g)%thetaDot2
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(n)) - block(g)%piv_y)  ! + ydot
         ENDIF

         !wsurf = 0._rk
         sur2nodeDis = block(g)%w2NormDis(n)
         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) + pt1*block(g)%cosAlpha(block(g)%nelw2(n))
         pos1_y = block(g)%yw(j) + pt1*block(g)%cosBeta(block(g)%nelw2(n))
         pos1_z = block(g)%zw(k+1) + pt1*block(g)%cosGamma(block(g)%nelw2(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xw(il).and.pos1_x<block(g)%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yw(jl).and.pos1_y<block(g)%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zw(kl).and.pos1_z<block(g)%zw(kl+1)) i_z1 = kl
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


         dwdn_e = dwdx_e*block(g)%cosAlpha(block(g)%nelw2(n)) + dwdy_e*block(g)%cosBeta(block(g)%nelw2(n)) +  dwdz_e*block(g)%cosGamma(block(g)%nelw2(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************

!**************************W(i,j,k-1)*************************************
         !wsurf = w_curr
         wsurf = 0._rk
         IF (block(g)%ibSurfID(block(g)%nelw1(n))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw1(n))==51) THEN
           block(g)%thetaDot = block(g)%thetaDot1
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(n)) - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw1(n))==52) THEN
           block(g)%thetaDot = block(g)%thetaDot2
           wsurf    = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(n)) - block(g)%piv_y)  ! + ydot
         ENDIF

         sur2nodeDis = block(g)%w1NormDis(n)

         pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) + pt1*block(g)%cosAlpha(block(g)%nelw1(n))
         pos1_y = block(g)%yw(j) + pt1*block(g)%cosBeta(block(g)%nelw1(n))
         pos1_z = block(g)%zw(k) + pt1*block(g)%cosGamma(block(g)%nelw1(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=block(g)%xw(il).and.pos1_x<block(g)%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=block(g)%yw(jl).and.pos1_y<block(g)%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=block(g)%zw(kl).and.pos1_z<block(g)%zw(kl+1)) i_z1 = kl
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


         dwdn_e = dwdx_e*block(g)%cosAlpha(block(g)%nelw1(n)) + dwdy_e*block(g)%cosBeta(block(g)%nelw1(n)) +  dwdz_e*block(g)%cosGamma(block(g)%nelw1(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w(i,j,k-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!***********************************************************************
      ENDDO
      !$acc end parallel
      ENDDO

END SUBROUTINE velocityForcingField
end module biocfd_forcing

module biocfd_forcing
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use biocfd_interpolation, only: linear_interpolation, bilinear_interpolation
  use biocfd_block_type, only: Block_t
  implicit none

  private

  public :: pressureForcing1, pressureforcingfield, pressureforcingghost
  public :: velocityforcing1, velocityforcingfield, velocityforcingghost
  public :: compute_value_and_derivatives

  contains
SUBROUTINE pressureForcing1(blk)
      type(Block_t), intent(inout) :: blk
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn, &
                   dpdn_e, ac_y, ac_z, at_y, at_z
      real(dp) :: derivatives(3)
      dpdn = 0._dp

 !$acc parallel loop gang vector                                                                                          &
 !$acc private (n, n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc           dpdn_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc           i_z1,ac_z,ac_y,at_y,at_z)         &
 !$acc default(present)  &
 !$acc private(derivatives)
      DO n = 1, blk%ibCellCount
         IF (blk%ibSurfId(blk%nelp(n))==50) THEN
            !dpdn = 0.
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !
        ELSEIF (blk%ibSurfId(blk%nelp(n))==51) THEN
            blk%thetaDot  = blk%thetaDot1
            blk%thetaDDot = blk%thetaDDot1
            ac_z = -blk%thetaDot**2*(blk%zcent(blk%nelp(n)) - blk%piv_z)
            ac_y = -blk%thetaDot**2*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_z =  blk%thetaDDot*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_y = -blk%thetaDDot*(blk%zcent(blk%nelp(n)) - blk%piv_z)
        ELSEIF (blk%ibSurfId(blk%nelp(n))==52) THEN
            blk%thetaDot  = blk%thetaDot2
            blk%thetaDDot = blk%thetaDDot2
            ac_z = -blk%thetaDot**2*(blk%zcent(blk%nelp(n)) - blk%piv_z)
            ac_y = -blk%thetaDot**2*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_z =  blk%thetaDDot*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_y = -blk%thetaDDot*(blk%zcent(blk%nelp(n)) - blk%piv_z)
        ENDIF
         dpdn = -((ac_z + at_z)*blk%cosGamma(blk%nelp(n)) + (ac_y + at_y) * &
                  blk%cosBeta(blk%nelp(n)))-blk%yddot * &
                     blk%cosBeta(blk%nelp(n))

         i = blk%interceptedIndexPtr(n, 1)
         j = blk%interceptedIndexPtr(n, 2)
         k = blk%interceptedIndexPtr(n, 3)

         sur2nodeDis = blk%pNormDis(n)

         pt1 = 1.5_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xp(i) + pt1*blk%cosAlpha(blk%nelp(n))
         pos1_y = blk%yp(j) + pt1*blk%cosBeta(blk%nelp(n))
         pos1_z = blk%zp(k) + pt1*blk%cosGamma(blk%nelp(n))

          !$acc loop seq
         DO il = 1, blk%nx+1
            if(pos1_x>=blk%xp(il).and.pos1_x<blk%xp(il+1)) i_x1 = il
         END DO
          !$acc loop seq
         DO jl = 1, blk%ny+1
            if(pos1_y>=blk%yp(jl).and.pos1_y<blk%yp(jl+1)) i_y1 = jl
         END DO
          !$acc loop seq
         DO kl = 1, blk%nz+1
            if(pos1_z>=blk%zp(kl).and.pos1_z<blk%zp(kl+1)) i_z1 = kl
         END DO

        call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                           blk%xp, blk%yp, blk%zp, 0, &
                                           blk%p, p_pos1, derivatives)

        dpdn_e = derivatives(1) * blk%cosAlpha(blk%nelp(n)) &
               + derivatives(2) * blk%cosBeta(blk%nelp(n)) &
               + derivatives(3) * blk%cosGamma(blk%nelp(n))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*0.5_dp

         blk%p(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
 !$acc end parallel loop

END SUBROUTINE pressureForcing1

SUBROUTINE velocityForcing1(blk)
      type(Block_t), intent(inout) :: blk
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, &
                         usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
                         dudn_e, dvdn_e, dwdn_e
       real(dp) :: derivatives(3)

 !$acc parallel loop gang vector         &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis,                       &
 !$acc          usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
 !$acc          dudn_e, dvdn_e, dwdn_e, &
 !$acc          k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc private(derivatives) &
 !$acc default(present)

      DO n = 1, blk%ibCellCount

         IF (blk%ibSurfID(blk%nelu2(n))==50) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfID(blk%nelu2(n))==51) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfId(blk%nelu2(n))==52) THEN
             usurf = blk%xdot
         ENDIF
         i = blk%interceptedIndexPtr(n, 1)
         j = blk%interceptedIndexPtr(n, 2)
         k = blk%interceptedIndexPtr(n, 3)

!***********************U(i,j,k)****************************************
         sur2nodeDis = blk%u2NormDis(n)

         pt1 = 1.5_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xu(i+1) + pt1*blk%cosAlpha(blk%nelu2(n))
         pos1_y = blk%yu(j)   + pt1*blk%cosBeta(blk%nelu2(n))
         pos1_z = blk%zu(k)   + pt1*blk%cosGamma(blk%nelu2(n))

         !$acc loop seq
         DO il = 1, blk%nx+2
            if(pos1_x>=blk%xu(il).and.pos1_x<blk%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, blk%ny+1
            if(pos1_y>=blk%yu(jl).and.pos1_y<blk%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, blk%nz+1
            if(pos1_z>=blk%zu(kl).and.pos1_z<blk%zu(kl+1)) i_z1 = kl
         END DO
         IF(i_x1==blk%nx+2) i_x1 = blk%nx+1

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xu, blk%yu, blk%zu, 1, &
                                            blk%ut, u_pos1, derivatives)

         dudn_e =  derivatives(1) * blk%cosAlpha(blk%nelu2(n)) &
                 + derivatives(2) * blk%cosBeta(blk%nelu2(n)) &
                 + derivatives(3) * blk%cosGamma(blk%nelu2(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         blk%ut(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval

!******************************U(i-1,j,k)*******************************
         IF (blk%ibSurfID(blk%nelu1(n))==50) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfID(blk%nelu1(n))==51) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfId(blk%nelu1(n))==52) THEN
             usurf = blk%xdot
         ENDIF
         sur2nodeDis = blk%u1NormDis(n)

         pt1 = 1.5_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xu(i) + pt1*blk%cosAlpha(blk%nelu1(n))
         pos1_y = blk%yu(j) + pt1*blk%cosBeta(blk%nelu1(n))
         pos1_z = blk%zu(k) + pt1*blk%cosGamma(blk%nelu1(n))

         !$acc loop seq
         DO il = 1, blk%nx+1
            if(pos1_x>=blk%xu(il).and.pos1_x<blk%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, blk%ny+1
            if(pos1_y>=blk%yu(jl).and.pos1_y<blk%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, blk%nz+1
            if(pos1_z>=blk%zu(kl).and.pos1_z<blk%zu(kl+1)) i_z1 = kl
         END DO

         IF(i_x1==1) i_x1 = 2

        call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xu, blk%yu, blk%zu, 1, &
                                            blk%ut, u_pos1, derivatives)

        dudn_e = derivatives(1) * blk%cosAlpha(blk%nelu1(n)) &
               + derivatives(2) * blk%cosBeta(blk%nelu1(n)) &
               + derivatives(3) * blk%cosGamma(blk%nelu1(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         blk%ut(i-1,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j,k)*************************************
         IF (blk%ibSurfID(blk%nelv2(n))==50) THEN
                 vsurf = blk%ydot
         ELSEIF (blk%ibSurfID(blk%nelv2(n))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv2(n)) &
                   - blk%piv_z) + blk%ydot
         ELSEIF (blk%ibSurfId(blk%nelv2(n))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv2(n)) &
                   - blk%piv_z)+ blk%ydot  ! + ydot
         ENDIF

         sur2nodeDis = blk%v2NormDis(n)

         pt1 = 1.5_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 &
               + blk%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xv(i) + pt1*blk%cosAlpha(blk%nelv2(n))
         pos1_y = blk%yv(j+1) + pt1*blk%cosBeta(blk%nelv2(n))
         pos1_z = blk%zv(k) + pt1*blk%cosGamma(blk%nelv2(n))

         !$acc loop seq
         DO il = 1, blk%nx+1
            if(pos1_x>=blk%xv(il).and.pos1_x<blk%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, blk%ny+2
            if(pos1_y>=blk%yv(jl).and.pos1_y<blk%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, blk%nz+1
            if(pos1_z>=blk%zv(kl).and.pos1_z<blk%zv(kl+1)) i_z1 = kl
         END DO

         IF(i_y1==blk%ny+2) i_y1 = blk%ny+1

        call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xv, blk%yv, blk%zv, 2, &
                                            blk%vt, v_pos1, derivatives)
        dvdn_e = derivatives(1) * blk%cosAlpha(blk%nelv2(n)) &
                  + derivatives(2) * blk%cosBeta(blk%nelv2(n)) &
                  + derivatives(3) * blk%cosGamma(blk%nelv2(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         blk%vt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j-1,k)*************************************
         IF (blk%ibSurfID(blk%nelv1(n))==50) THEN
           vsurf = 0._dp+ blk%ydot
         ELSEIF (blk%ibSurfID(blk%nelv1(n))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           vsurf    = -blk%thetaDot*(blk%zcent(blk%nelv1(n)) - blk%piv_z) &
                      + blk%ydot
         ELSEIF (blk%ibSurfId(blk%nelv1(n))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           vsurf    = -blk%thetaDot*(blk%zcent(blk%nelv1(n)) - blk%piv_z) &
                      + blk%ydot  ! + ydot
         ENDIF
         sur2nodeDis = blk%v1NormDis(n)

         pt1 = 1.5_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xv(i) + pt1*blk%cosAlpha(blk%nelv1(n))
         pos1_y = blk%yv(j) + pt1*blk%cosBeta(blk%nelv1(n))
         pos1_z = blk%zv(k) + pt1*blk%cosGamma(blk%nelv1(n))

         !$acc loop seq
         DO il = 1, blk%nx+1
            if(pos1_x>=blk%xv(il).and.pos1_x<blk%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, blk%ny+1
            if(pos1_y>=blk%yv(jl).and.pos1_y<blk%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, blk%nz+1
            if(pos1_z>=blk%zv(kl).and.pos1_z<blk%zv(kl+1)) i_z1 = kl
         END DO

         IF(i_y1==1) i_y1 = 2

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xv, blk%yv, blk%zv, 2, &
                                            blk%vt, v_pos1, derivatives)

           dvdn_e = derivatives(1) * blk%cosAlpha(blk%nelv1(n)) &
                  + derivatives(2) * blk%cosBeta(blk%nelv1(n)) &
                  + derivatives(3) * blk%cosGamma(blk%nelv1(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         blk%vt(i,j-1,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k)*************************************
         IF (blk%ibSurfID(blk%nelw2(n))==50) THEN
           wsurf = 0.
         ELSEIF (blk%ibSurfID(blk%nelw2(n))==51) THEN
           blk%thetaDot = blk%thetaDot1
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw2(n)) - blk%piv_y)
         ELSEIF (blk%ibSurfId(blk%nelw2(n))==52) THEN
           blk%thetaDot = blk%thetaDot2
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw2(n)) - blk%piv_y)  ! + ydot
         ENDIF
         sur2nodeDis = blk%w2NormDis(n)

         pt1 = 1.5_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xw(i) + pt1*blk%cosAlpha(blk%nelw2(n))
         pos1_y = blk%yw(j) + pt1*blk%cosBeta(blk%nelw2(n))
         pos1_z = blk%zw(k+1) + pt1*blk%cosGamma(blk%nelw2(n))

         !$acc loop seq
         DO il = 1, blk%nx+1
            if(pos1_x>=blk%xw(il).and.pos1_x<blk%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, blk%ny+1
            if(pos1_y>=blk%yw(jl).and.pos1_y<blk%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, blk%nz+2
            if(pos1_z>=blk%zw(kl).and.pos1_z<blk%zw(kl+1)) i_z1 = kl
         END DO

         IF(i_z1==blk%nz+2) i_z1 = blk%nz+1

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xw, blk%yw, blk%zw, 3, &
                                            blk%wt, w_pos1, derivatives)
          dwdn_e = derivatives(1) *blk%cosAlpha(blk%nelw2(n)) &
                 + derivatives(2) *blk%cosBeta(blk%nelw2(n)) &
                 + derivatives(3) *blk%cosGamma(blk%nelw2(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         blk%wt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k-1)*************************************
         wsurf = 0._dp
         IF (blk%ibSurfID(blk%nelw1(n))==50) THEN
           wsurf = 0.
         ELSEIF (blk%ibSurfID(blk%nelw1(n))==51) THEN
           blk%thetaDot = blk%thetaDot1
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw1(n)) - blk%piv_y)
         ELSEIF (blk%ibSurfId(blk%nelw1(n))==52) THEN
           blk%thetaDot = blk%thetaDot2
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw1(n)) - blk%piv_y)  ! + ydot
         ENDIF
         sur2nodeDis = blk%w1NormDis(n)

         pt1 = 1.5_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xw(i) + pt1*blk%cosAlpha(blk%nelw1(n))
         pos1_y = blk%yw(j) + pt1*blk%cosBeta(blk%nelw1(n))
         pos1_z = blk%zw(k) + pt1*blk%cosGamma(blk%nelw1(n))

         !$acc loop seq
         DO il = 1, blk%nx+1
            if(pos1_x>=blk%xw(il).and.pos1_x<blk%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = 1, blk%ny+1
            if(pos1_y>=blk%yw(jl).and.pos1_y<blk%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, blk%nz+1
            if(pos1_z>=blk%zw(kl).and.pos1_z<blk%zw(kl+1)) i_z1 = kl
         END DO

         IF(i_z1==1) i_z1 = 2

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xw, blk%yw, blk%zw, 3, &
                                            blk%wt, w_pos1, derivatives)

         dwdn_e =    derivatives(1) * blk%cosAlpha(blk%nelw1(n)) &
                   + derivatives(2) * blk%cosBeta(blk%nelw1(n)) &
                   + derivatives(3) * blk%cosGamma(blk%nelw1(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         blk%wt(i,j,k-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
!$acc end parallel loop

END SUBROUTINE velocityForcing1

SUBROUTINE pressureForcingGhost(blk)
      type(Block_t), intent(inout) :: blk
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1, &
                         aval, bval, cval, p_pos1, sur2nodeDis, dpdn, &
                         dpdn_e, ac_y, ac_z, at_y, at_z
       real(dp) :: derivatives(3)

      dpdn = 0._dp
 !$acc parallel loop gang vector                                                                    &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc           dpdn_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc           i_z1,ac_z,ac_y,at_y,at_z)         &
 !$acc default(present) private(derivatives)
      DO n = 1, blk%TSCellCount

        i = blk%TSIndexPtr(n, 1)
        j = blk%TSIndexPtr(n, 2)
        k = blk%TSIndexPtr(n, 3)
        IF (blk%ibSurfId(blk%nelp(blk%index_ts(n)))==50) THEN
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !-block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
        ELSEIF (blk%ibSurfId(blk%nelp(blk%index_ts(n)))==51) THEN
            blk%thetaDot  = blk%thetaDot1
            blk%thetaDDot = blk%thetaDDot1
            ac_z = -blk%thetaDot**2*(blk%zcent(blk%nelp(blk%index_ts(n))) &
                   - blk%piv_z)
            ac_y = -blk%thetaDot**2*(blk%ycent(blk%nelp(blk%index_ts(n))) &
                   - blk%piv_y)
            at_z =  blk%thetaDDot*(blk%ycent(blk%nelp(blk%index_ts(n))) &
                   - blk%piv_y)
            at_y = -blk%thetaDDot*(blk%zcent(blk%nelp(blk%index_ts(n))) &
                   - blk%piv_z)

        ELSEIF (blk%ibSurfId(blk%nelp(blk%index_ts(n)))==52) THEN
            blk%thetaDot  = blk%thetaDot2
            blk%thetaDDot = blk%thetaDDot2
            ac_z = -blk%thetaDot**2*(blk%zcent(blk%nelp(blk%index_ts(n))) &
                   - blk%piv_z)
            ac_y = -blk%thetaDot**2*(blk%ycent(blk%nelp(blk%index_ts(n))) &
                   - blk%piv_y)
            at_z =  blk%thetaDDot*(blk%ycent(blk%nelp(blk%index_ts(n)))  &
                   - blk%piv_y)
            at_y = -blk%thetaDDot*(blk%zcent(blk%nelp(blk%index_ts(n))) &
                   - blk%piv_z)
        ENDIF
            dpdn = ((ac_z + at_z)* blk%cosGamma(blk%nelp(blk%index_ts(n))) &
                  + (ac_y + at_y)* blk%cosBeta(blk%nelp(blk%index_ts(n)))) &
                   + (blk%yddot*(blk%cosBeta(blk%index_ts(n))))

         sur2nodeDis = -blk%pNormDis(blk%index_ts(n))

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 &
                             + blk%deltay(j)**2 &
                             + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xp(i) - pt1*blk%cosAlpha(blk%nelp(blk%index_ts(n)))
         pos1_y = blk%yp(j) - pt1*blk%cosBeta(blk%nelp(blk%index_ts(n)))
         pos1_z = blk%zp(k) - pt1*blk%cosGamma(blk%nelp(blk%index_ts(n)))
         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xp(il).and.pos1_x<blk%xp(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yp(jl).and.pos1_y<blk%yp(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zp(kl).and.pos1_z<blk%zp(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                           blk%xp, blk%yp, blk%zp, 0, &
                                           blk%p, p_pos1, derivatives)

         dpdn_e =  -1 * (derivatives(1) * blk%cosAlpha(blk%nelp(blk%index_ts(n))) &
                        + derivatives(2) * blk%cosBeta(blk%nelp(blk%index_ts(n))) &
                        + derivatives(3) * blk%cosGamma(blk%nelp(blk%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*0.5_dp

         blk%p_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         blk%pt_ghost(n) = blk%p(i,j,k)
      ENDDO
      !$acc end parallel loop
END SUBROUTINE pressureForcingGhost

SUBROUTINE velocityForcingGhost(blk)
      type(Block_t), intent(inout) :: blk
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, &
                         usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
                         dudn_e, dvdn_e, dwdn_e
      real(dp) :: derivatives(3)

 !$acc parallel loop gang vector         &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis, &
 !$acc          usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
 !$acc          dudn_e, dvdn_e, dwdn_e, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc default(present) private(derivatives)
      DO n = 1, blk%TSCellCount

        i = blk%TSIndexPtr(n, 1)
        j = blk%TSIndexPtr(n, 2)
        k = blk%TSIndexPtr(n, 3)
        !IF (block(g)%cell2(i,j,k).EQ.2) THEN
!***********************U(i,j,k)****************************************
         IF (blk%ibSurfID(blk%nelu2(blk%index_ts(n)))==50) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfID(blk%nelu2(blk%index_ts(n)))==51) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfId(blk%nelu2(blk%index_ts(n)))==52) THEN
           usurf = blk%xdot
         ENDIF

         sur2nodeDis = -blk%u2NormDis(blk%index_ts(n))

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 &
               + blk%deltay(j)**2 &
               + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xu(i+1) - pt1*blk%cosAlpha(blk%nelu2(blk%index_ts(n)))
         pos1_y = blk%yu(j) - pt1*blk%cosBeta(blk%nelu2(blk%index_ts(n)))
         pos1_z = blk%zu(k) - pt1*blk%cosGamma(blk%nelu2(blk%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xu(il).and.pos1_x<blk%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yu(jl).and.pos1_y<blk%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zu(kl).and.pos1_z<blk%zu(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xu, blk%yu, blk%zu, &
                                            1, blk%ut, u_pos1, derivatives)

         dudn_e = -1 * (derivatives(1)*blk%cosAlpha(blk%nelu2(blk%index_ts(n))) &
                  + derivatives(2)*blk%cosBeta(blk%nelu2(blk%index_ts(n))) &
                  + derivatives(3)*blk%cosGamma(blk%nelu2(blk%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         blk%u2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         blk%u2t_ghost(n) = blk%u(i,j,k)
!******************************U(i-1,j,k)*******************************
         IF (blk%ibSurfID(blk%nelu1(blk%index_ts(n)))==50) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfID(blk%nelu1(blk%index_ts(n)))==51) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfId(blk%nelu1(blk%index_ts(n)))==52) THEN
           usurf = blk%xdot
         ENDIF

         sur2nodeDis = -blk%u1NormDis(blk%index_ts(n))

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xu(i) - pt1*blk%cosAlpha(blk%nelu1(blk%index_ts(n)))
         pos1_y = blk%yu(j) - pt1*blk%cosBeta(blk%nelu1(blk%index_ts(n)))
         pos1_z = blk%zu(k) - pt1*blk%cosGamma(blk%nelu1(blk%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xu(il).and.pos1_x<blk%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yu(jl).and.pos1_y<blk%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zu(kl).and.pos1_z<blk%zu(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xu, blk%yu, blk%zu, &
                                            1, blk%ut, u_pos1, derivatives)

         dudn_e = -1 * (derivatives(1) * blk%cosAlpha(blk%nelu1(blk%index_ts(n))) &
                      + derivatives(2) * blk%cosBeta(blk%nelu1(blk%index_ts(n))) &
                      + derivatives(3) * blk%cosGamma(blk%nelu1(blk%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         blk%u1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         blk%u1t_ghost(n) = blk%u(i-1,j,k)
!**************************V(i,j,k)*************************************
         IF (blk%ibSurfID(blk%nelv2(blk%index_ts(n)))==50) THEN
           vsurf = 0._dp+ blk%ydot
         ELSEIF (blk%ibSurfID(blk%nelv2(blk%index_ts(n)))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv2(blk%index_ts(n))) &
                   - blk%piv_z) &
                   + (-blk%alphaDot)*(blk%xcent(blk%nelv2(blk%index_ts(n))) &
                   - blk%piv_x)+ blk%ydot
         ELSEIF (blk%ibSurfId(blk%nelv2(blk%index_ts(n)))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv2(blk%index_ts(n))) &
                   - blk%piv_z) &
                   + (-blk%alphaDot)*(blk%xcent(blk%nelv2(blk%index_ts(n))) &
                   - blk%piv_x) + blk%ydot  ! + ydot
         ENDIF

         sur2nodeDis = -blk%v2NormDis(blk%index_ts(n))

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xv(i) - pt1*blk%cosAlpha(blk%nelv2(blk%index_ts(n)))
         pos1_y = blk%yv(j+1) - pt1*blk%cosBeta(blk%nelv2(blk%index_ts(n)))
         pos1_z = blk%zv(k) - pt1*blk%cosGamma(blk%nelv2(blk%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xv(il).and.pos1_x<blk%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yv(jl).and.pos1_y<blk%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zv(kl).and.pos1_z<blk%zv(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xv, blk%yv, blk%zv, &
                                            2, blk%vt, v_pos1, derivatives)
          dvdn_e = -1 * (derivatives(1) * blk%cosAlpha(blk%nelv2(blk%index_ts(n))) &
                       + derivatives(2) * blk%cosBeta(blk%nelv2(blk%index_ts(n))) &
                       + derivatives(3) * blk%cosGamma(blk%nelv2(blk%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         blk%v2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%v2_ghost(n) = v(i,j,k)
         !ENDIF
         blk%v2t_ghost(n) = blk%v(i,j,k)
!**************************V(i,j-1,k)*************************************
         IF (blk%ibSurfID(blk%nelv1(blk%index_ts(n)))==50) THEN
           vsurf = blk%ydot
         ELSEIF (blk%ibSurfID(blk%nelv1(blk%index_ts(n)))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           vsurf    = -blk%thetaDot*(blk%zcent(blk%nelv1(blk%index_ts(n))) &
                      - blk%piv_z)  &
                      + (-blk%alphaDot)*(blk%xcent(blk%nelv1(blk%index_ts(n)))&
                      - blk%piv_x) + blk%ydot
         ELSEIF (blk%ibSurfId(blk%nelv1(blk%index_ts(n)))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           vsurf    = -blk%thetaDot*(blk%zcent(blk%nelv1(blk%index_ts(n))) &
                      - blk%piv_z)  &
                      + (-blk%alphaDot)*(blk%xcent(blk%nelv1(blk%index_ts(n)))&
                      - blk%piv_x)+ blk%ydot  ! + ydot
         ENDIF

         sur2nodeDis = -blk%v1NormDis(blk%index_ts(n))

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 &
               + blk%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xv(i) - pt1*blk%cosAlpha(blk%nelv1(blk%index_ts(n)))
         pos1_y = blk%yv(j) - pt1*blk%cosBeta(blk%nelv1(blk%index_ts(n)))
         pos1_z = blk%zv(k) - pt1*blk%cosGamma(blk%nelv1(blk%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xv(il).and.pos1_x<blk%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yv(jl).and.pos1_y<blk%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zv(kl).and.pos1_z<blk%zv(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xv, blk%yv, blk%zv, &
                                            2, blk%vt, v_pos1, derivatives)
          dvdn_e = -1 * (derivatives(1) * blk%cosAlpha(blk%nelv1(blk%index_ts(n))) &
                       + derivatives(2) * blk%cosBeta(blk%nelv1(blk%index_ts(n))) &
                       + derivatives(3) * blk%cosGamma(blk%nelv1(blk%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         blk%v1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         blk%v1t_ghost(n) = blk%v(i,j-1,k)
!**************************W(i,j,k)*************************************
         wsurf = 0._dp
         IF (blk%ibSurfID(blk%nelw2(blk%index_ts(n)))==50) THEN
           wsurf = 0.
         ELSEIF (blk%ibSurfID(blk%nelw2(blk%index_ts(n)))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           wsurf = blk%thetaDot*(blk%ycent(blk%nelw2(blk%index_ts(n))) &
                   - blk%piv_y)
         ELSEIF (blk%ibSurfId(blk%nelw2(blk%index_ts(n)))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           wsurf = blk%thetaDot*(blk%ycent(blk%nelw2(blk%index_ts(n))) &
                   - blk%piv_y)  ! + ydot
         ENDIF


         sur2nodeDis = -blk%w2NormDis(blk%index_ts(n))

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xw(i) - pt1*blk%cosAlpha(blk%nelw2(blk%index_ts(n)))
         pos1_y = blk%yw(j) - pt1*blk%cosBeta(blk%nelw2(blk%index_ts(n)))
         pos1_z = blk%zw(k+1) - pt1*blk%cosGamma(blk%nelw2(blk%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xw(il).and.pos1_x<blk%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yw(jl).and.pos1_y<blk%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zw(kl).and.pos1_z<blk%zw(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xw, blk%yw, blk%zw, &
                                            3, blk%wt, w_pos1, derivatives)
         dwdn_e = -1 * (derivatives(1) * blk%cosAlpha(blk%nelw2(blk%index_ts(n))) &
                      + derivatives(2) * blk%cosBeta(blk%nelw2(blk%index_ts(n))) &
                      + derivatives(3) * blk%cosGamma(blk%nelw2(blk%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         blk%w2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%w2_ghost(n) = w(i,j,k)
         !ENDIF
         blk%w2t_ghost(n) = blk%w(i,j,k)
!**************************W(i,j,k-1)*************************************
         wsurf = 0._dp
         IF (blk%ibSurfID(blk%nelw1(blk%index_ts(n)))==50) THEN
           wsurf = 0.
         ELSEIF (blk%ibSurfID(blk%nelw1(blk%index_ts(n)))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           wsurf = blk%thetaDot*(blk%ycent(blk%nelw1(blk%index_ts(n))) &
                   - blk%piv_y)
         ELSEIF (blk%ibSurfId(blk%nelw1(blk%index_ts(n)))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           wsurf =  blk%thetaDot*(blk%ycent(blk%nelw1(blk%index_ts(n))) &
                   - blk%piv_y)  ! + ydot
         ENDIF


         sur2nodeDis = -blk%w1NormDis(blk%index_ts(n))

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xw(i) - pt1*blk%cosAlpha(blk%nelw1(blk%index_ts(n)))
         pos1_y = blk%yw(j) - pt1*blk%cosBeta(blk%nelw1(blk%index_ts(n)))
         pos1_z = blk%zw(k) - pt1*blk%cosGamma(blk%nelw1(blk%index_ts(n)))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xw(il).and.pos1_x<blk%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yw(jl).and.pos1_y<blk%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zw(kl).and.pos1_z<blk%zw(kl+1)) i_z1 = kl
         END DO
         !IF (cell(i_x1, i_y1, i_z1).EQ.0) THEN
         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xw, blk%yw, blk%zw, &
                                            3, blk%wt, w_pos1, derivatives)
         dwdn_e = -1 * (derivatives(1) * blk%cosAlpha(blk%nelw1(blk%index_ts(n))) &
                      + derivatives(2) * blk%cosBeta(blk%nelw1(blk%index_ts(n))) &
                      + derivatives(3) * blk%cosGamma(blk%nelw1(blk%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         blk%w1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         blk%w1t_ghost(n) = blk%w(i,j,k-1)
      ENDDO
      !$acc end parallel loop

END SUBROUTINE velocityForcingGhost
!***********************************************************************

SUBROUTINE pressureForcingField(blk)
      type(Block_t), intent(inout) :: blk
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1, &
                   aval, bval, cval, p_pos1, sur2nodeDis, dpdn, &
                   dpdn_e, ac_y, ac_z, at_y, at_z
       real(dp) :: derivatives(3)

      dpdn = 0._dp
 !$acc parallel loop gang vector                                                                                          &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc          dpdn_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc          i_z1,ac_z,ac_y,at_y,at_z)         &
 !$acc default(present) private(derivatives)
      DO n = 1, blk%ibCellCount
         IF (blk%ibSurfId(blk%nelp(n))==50) THEN
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !
        ELSEIF (blk%ibSurfId(blk%nelp(n))==51) THEN
            blk%thetaDot  = blk%thetaDot1
            blk%thetaDDot = blk%thetaDDot1
            ac_z = -blk%thetaDot**2*(blk%zcent(blk%nelp(n)) - blk%piv_z)
            ac_y = -blk%thetaDot**2*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_z =  blk%thetaDDot*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_y = -blk%thetaDDot*(blk%zcent(blk%nelp(n)) - blk%piv_z)
        ELSEIF (blk%ibSurfId(blk%nelp(n))==52) THEN
            blk%thetaDot  = blk%thetaDot2
            blk%thetaDDot = blk%thetaDDot2
            ac_z = -blk%thetaDot**2*(blk%zcent(blk%nelp(n)) - blk%piv_z)
            ac_y = -blk%thetaDot**2*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_z =  blk%thetaDDot*(blk%ycent(blk%nelp(n)) - blk%piv_y)
            at_y = -blk%thetaDDot*(blk%zcent(blk%nelp(n)) - blk%piv_z)
        ENDIF
         dpdn = -((ac_z + at_z)*blk%cosGamma(blk%nelp(n)) &
                + (ac_y + at_y)*blk%cosBeta(blk%nelp(n))) &
                -blk%yddot*blk%cosBeta(blk%nelp(n))

         i = blk%interceptedIndexPtr(n, 1)
         j = blk%interceptedIndexPtr(n, 2)
         k = blk%interceptedIndexPtr(n, 3)

         sur2nodeDis = blk%pNormDis(n)

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xp(i) + pt1*blk%cosAlpha(blk%nelp(n))
         pos1_y = blk%yp(j) + pt1*blk%cosBeta(blk%nelp(n))
         pos1_z = blk%zp(k) + pt1*blk%cosGamma(blk%nelp(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xp(il).and.pos1_x<blk%xp(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yp(jl).and.pos1_y<blk%yp(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zp(kl).and.pos1_z<blk%zp(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                           blk%xp, blk%yp, blk%zp, 0, &
                                           blk%p, p_pos1, derivatives)

         dpdn_e  = derivatives(1) * blk%cosAlpha(blk%nelp(n)) &
                 + derivatives(2) * blk%cosBeta(blk%nelp(n))  &
                 + derivatives(3) * blk%cosGamma(blk%nelp(n))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*0.5_dp

         blk%p(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
      !$acc end parallel loop
END SUBROUTINE pressureForcingField

SUBROUTINE velocityForcingField(blk)
      type(Block_t), intent(inout) :: blk
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, &
                         usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
                         dudn_e, dvdn_e, dwdn_e
      real(dp) :: derivatives(3)

 !$acc parallel loop gang vector         &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis, &
 !$acc          usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
 !$acc          dudn_e, dvdn_e, dwdn_e, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc default(present) private(derivatives)
      DO n = 1, blk%ibCellCount

         i = blk%interceptedIndexPtr(n, 1)
         j = blk%interceptedIndexPtr(n, 2)
         k = blk%interceptedIndexPtr(n, 3)

!***********************U(i,j,k)****************************************
         IF (blk%ibSurfID(blk%nelu2(n))==50) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfID(blk%nelu2(n))==51) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfId(blk%nelu2(n))==52) THEN
             usurf = blk%xdot
         ENDIF

         sur2nodeDis = blk%u2NormDis(n)

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xu(i+1) + pt1*blk%cosAlpha(blk%nelu2(n))
         pos1_y = blk%yu(j)   + pt1*blk%cosBeta(blk%nelu2(n))
         pos1_z = blk%zu(k)   + pt1*blk%cosGamma(blk%nelu2(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xu(il).and.pos1_x<blk%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yu(jl).and.pos1_y<blk%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zu(kl).and.pos1_z<blk%zu(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xu, blk%yu, blk%zu, &
                                            1, blk%ut, u_pos1, derivatives)

         dudn_e =   derivatives(1) * blk%cosAlpha(blk%nelu2(n)) &
                  + derivatives(2) * blk%cosBeta(blk%nelu2(n)) &
                  + derivatives(3) * blk%cosGamma(blk%nelu2(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         blk%u(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!******************************U(i-1,j,k)*******************************
         IF (blk%ibSurfID(blk%nelu1(n))==50) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfID(blk%nelu1(n))==51) THEN
             usurf = blk%xdot
         ELSEIF (blk%ibSurfId(blk%nelu1(n))==52) THEN
             usurf = blk%xdot
         ENDIF

         sur2nodeDis = blk%u1NormDis(n)

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 &
               + blk%deltay(j)**2 + blk%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xu(i) + pt1*blk%cosAlpha(blk%nelu1(n))
         pos1_y = blk%yu(j) + pt1*blk%cosBeta(blk%nelu1(n))
         pos1_z = blk%zu(k) + pt1*blk%cosGamma(blk%nelu1(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xu(il).and.pos1_x<blk%xu(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yu(jl).and.pos1_y<blk%yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zu(kl).and.pos1_z<blk%zu(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xu, blk%yu, blk%zu, &
                                            1, blk%ut, u_pos1, derivatives)

         dudn_e =   derivatives(1) * blk%cosAlpha(blk%nelu1(n)) &
                  + derivatives(2) * blk%cosBeta(blk%nelu1(n)) &
                  + derivatives(3) * blk%cosGamma(blk%nelu1(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         blk%u(i-1,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j,k)*************************************
         IF (blk%ibSurfID(blk%nelv2(n))==50) THEN
           vsurf = blk%ydot
         ELSEIF (blk%ibSurfID(blk%nelv2(n))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv2(n)) - blk%piv_z) &
                   + (-blk%alphaDot)*(blk%xcent(blk%nelv2(n)) - blk%piv_x) &
                   + blk%ydot
         ELSEIF (blk%ibSurfId(blk%nelv2(n))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv2(n)) - blk%piv_z) &
                   + (-blk%alphaDot)*( blk%xcent(blk%nelv2(n)) - blk%piv_x) &
                   + blk%ydot
         ENDIF

         sur2nodeDis = blk%v2NormDis(n)

         pt1 = 1.21_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xv(i) + pt1*blk%cosAlpha(blk%nelv2(n))
         pos1_y = blk%yv(j+1) + pt1*blk%cosBeta(blk%nelv2(n))
         pos1_z = blk%zv(k) + pt1*blk%cosGamma(blk%nelv2(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xv(il).and.pos1_x<blk%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yv(jl).and.pos1_y<blk%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zv(kl).and.pos1_z<blk%zv(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xv, blk%yv, blk%zv, &
                                            2, blk%vt, v_pos1, derivatives)
         dvdn_e =  derivatives(1) * blk%cosAlpha(blk%nelv2(n)) &
                 + derivatives(2) * blk%cosBeta(blk%nelv2(n)) &
                 + derivatives(3) * blk%cosGamma(blk%nelv2(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         blk%v(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j-1,k)*************************************
         IF (blk%ibSurfID(blk%nelv1(n))==50) THEN
           vsurf = blk%ydot
         ELSEIF (blk%ibSurfID(blk%nelv1(n))==51) THEN
           blk%thetaDot =  blk%thetaDot1
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv1(n)) - blk%piv_z)  &
                   + (-blk%alphaDot)*(blk%xcent(blk%nelv1(n)) - blk%piv_x) &
                   + blk%ydot
         ELSEIF (blk%ibSurfId(blk%nelv1(n))==52) THEN
           blk%thetaDot =  blk%thetaDot2
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)! + ydot
           vsurf = -blk%thetaDot*(blk%zcent(blk%nelv1(n)) - blk%piv_z)  &
                   + (-blk%alphaDot)*(blk%xcent(blk%nelv1(n)) - blk%piv_x) &
                   + blk%ydot
         ENDIF

         !vsurf =  -block(g)%thetaDot*(block(g)%xcent(block(g)%nelv1(n)) - block(g)%piv_x) + ydot
         sur2nodeDis = blk%v1NormDis(n)

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
              + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xv(i) + pt1*blk%cosAlpha(blk%nelv1(n))
         pos1_y = blk%yv(j) + pt1*blk%cosBeta(blk%nelv1(n))
         pos1_z = blk%zv(k) + pt1*blk%cosGamma(blk%nelv1(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xv(il).and.pos1_x<blk%xv(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yv(jl).and.pos1_y<blk%yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zv(kl).and.pos1_z<blk%zv(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xv, blk%yv, blk%zv, &
                                            2, blk%vt, v_pos1, derivatives)
         dvdn_e = derivatives(1) * blk%cosAlpha(blk%nelv1(n)) &
                + derivatives(2) * blk%cosBeta(blk%nelv1(n)) &
                + derivatives(3) * blk%cosGamma(blk%nelv1(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         blk%v(i,j-1,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k)*************************************
         IF (blk%ibSurfID(blk%nelw2(n))==50) THEN
           wsurf = 0.
         ELSEIF (blk%ibSurfID(blk%nelw2(n))==51) THEN
           blk%thetaDot = blk%thetaDot1
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw2(n)) - blk%piv_y)
         ELSEIF (blk%ibSurfId(blk%nelw2(n))==52) THEN
           blk%thetaDot = blk%thetaDot2
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw2(n)) - blk%piv_y)  ! + ydot
         ENDIF

         sur2nodeDis = blk%w2NormDis(n)
         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xw(i) + pt1*blk%cosAlpha(blk%nelw2(n))
         pos1_y = blk%yw(j) + pt1*blk%cosBeta(blk%nelw2(n))
         pos1_z = blk%zw(k+1) + pt1*blk%cosGamma(blk%nelw2(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xw(il).and.pos1_x<blk%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yw(jl).and.pos1_y<blk%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zw(kl).and.pos1_z<blk%zw(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xw, blk%yw, blk%zw, &
                                            3, blk%wt, w_pos1, derivatives)
         dwdn_e = derivatives(1) * blk%cosAlpha(blk%nelw2(n)) &
                + derivatives(2) * blk%cosBeta(blk%nelw2(n)) &
                + derivatives(3) * blk%cosGamma(blk%nelw2(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         blk%w(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k-1)*************************************
         wsurf = 0._dp
         IF (blk%ibSurfID(blk%nelw1(n))==50) THEN
           wsurf = 0.
         ELSEIF (blk%ibSurfID(blk%nelw1(n))==51) THEN
           blk%thetaDot = blk%thetaDot1
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw1(n)) - blk%piv_y)
         ELSEIF (blk%ibSurfId(blk%nelw1(n))==52) THEN
           blk%thetaDot = blk%thetaDot2
           wsurf    = blk%thetaDot*(blk%ycent(blk%nelw1(n)) - blk%piv_y)  ! + ydot
         ENDIF

         sur2nodeDis = blk%w1NormDis(n)

         pt1 = 1.51_dp*dsqrt(blk%deltax(i)**2 + blk%deltay(j)**2 + blk%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = blk%xw(i) + pt1*blk%cosAlpha(blk%nelw1(n))
         pos1_y = blk%yw(j) + pt1*blk%cosBeta(blk%nelw1(n))
         pos1_z = blk%zw(k) + pt1*blk%cosGamma(blk%nelw1(n))

         !$acc loop seq
         DO il = i-7, i+7
            if(pos1_x>=blk%xw(il).and.pos1_x<blk%xw(il+1)) i_x1 = il
         END DO
         !$acc loop seq
         DO jl = j-7, j+7
            if(pos1_y>=blk%yw(jl).and.pos1_y<blk%yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = k-7, k+7
            if(pos1_z>=blk%zw(kl).and.pos1_z<blk%zw(kl+1)) i_z1 = kl
         END DO

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            blk%xw, blk%yw, blk%zw, &
                                            3, blk%wt, w_pos1, derivatives)
         dwdn_e = derivatives(1) * blk%cosAlpha(blk%nelw1(n)) &
                + derivatives(2) * blk%cosBeta(blk%nelw1(n)) &
                + derivatives(3) * blk%cosGamma(blk%nelw1(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         blk%w(i,j,k-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
      !$acc end parallel loop

END SUBROUTINE velocityForcingField

!> Derivatives are computed using Lagrange polynomials, see
!> https://en.wikipedia.org/wiki/Lagrange_polynomial for more
!> information about how these are constructed.
pure function compute_derivative(x, x2, x1, p_x, p_x2, p_x1) result(out)
  real(dp), intent(in) :: x, x2, x1, p_x, p_x2, p_x1
  real(dp) :: out

  real(dp) :: h1, h2

  h2 = abs(x2 - x)
  h1 = abs(x1 - x)
  out = (h1**2*p_x2 - h2**2*p_x1 + (h2**2- h1**2)*p_x)/(h1*h2*(h1+h2)+1e-16_dp)
end function compute_derivative

subroutine compute_value_and_derivatives(x, y, z, i, j, k, xgrid, ygrid, zgrid, offset_axis, &
                                         var, val, derivatives)
  !$acc routine
  real(dp), intent(in) :: x, y, z
  integer, intent(in) :: i, j, k
  real(dp), intent(in), dimension(:) :: xgrid, ygrid, zgrid
  !> An integer in the range 0-3 (inclusive) determining which axis to
  !> offset A value of 0 corresponds to no offset, 1=x, 2=y, 3=z. This
  !> is used for the calculations of p, u, v, and w.
  integer, intent(in) :: offset_axis
  real(dp), intent(in), dimension(:, :, :) :: var

  real(dp), intent(out) :: val
  real(dp), intent(out) :: derivatives(3)

  ! Internal variables
  ! tmp to hold interpolated results, [[x1, x2], [y1, y2], [z1, z1]]
  real (dp) :: tmp(3, 2)

   !> Indicies for the variable (which may be offset)
  integer :: ii, jj, kk

  ii = i
  jj = j
  kk = k

  if (offset_axis == 1) ii = i - 1
  if (offset_axis == 2) jj = j - 1
  if (offset_axis == 3) kk = k - 1

  ! Compute bilinear interpolation on six faces of a cuboid.
  !
  ! Given the eight corners of a cuboid find the interpolated values on
  ! the faces of the cuboid at some target position e.g. if the corners
  ! are at [(x1,y1,z1), (x2,y1,z1), (x1,y2,z1), (x2,y2,z1), (x1,y1,z2),
  ! (x2,y1,z2), (x1,y2,z2), (x2,y2,z2)], and the target position is (x,
  ! y, z), returns the interpolated values of something at [(x,y,z1),
  ! (x,y,z2), (x,y1,z), (x,y2,z), (x1,y,z), (x2,y,z)].
  !
  ! Note: the multi bilinear interpolation was originally in its own
  ! function, but nvfortran had problems inlining it so it is "by
  ! hand" inlined here
  !
  ! p_x1 y -- z plane
  tmp(1, 1) = bilinear_interpolation(y, z, ygrid(j), ygrid(j+1), zgrid(k), zgrid(k+1), &
              [var(ii, jj, kk), var(ii, jj+1, kk), var(ii, jj, kk+1), var(ii, jj+1, kk+1)])
  ! p_x2 y -- z plane
  tmp(1, 2) = bilinear_interpolation(y, z, ygrid(j), ygrid(j+1), zgrid(k), zgrid(k+1), &
              [var(ii+1, jj, kk), var(ii+1, jj+1, kk), var(ii+1, jj, kk+1), var(ii+1, jj+1, kk+1)])

  ! p_y1 x -- z plane
  tmp(2, 1) = bilinear_interpolation(x, z, xgrid(i), xgrid(i+1), zgrid(k), zgrid(k+1), &
              [var(ii, jj, kk), var(ii+1, jj, kk), var(ii, jj, kk+1), var(ii+1, jj, kk+1)])
  ! p_y2 x -- z plane
  tmp(2, 2) = bilinear_interpolation(x, z, xgrid(i), xgrid(i+1), zgrid(k), zgrid(k+1), &
              [var(ii, jj+1, kk), var(ii+1, jj+1, kk), var(ii, jj+1, kk+1), var(ii+1, jj+1, kk+1)])

  ! p_z1 x -- y plane
  tmp(3, 1) = bilinear_interpolation(x, y, xgrid(i), xgrid(i+1), ygrid(j), ygrid(j+1), &
              [var(ii, jj, kk), var(ii+1, jj, kk), var(ii, jj+1, kk), var(ii+1, jj+1, kk)])
  ! p_z2 x -- y plane
  tmp(3, 2) = bilinear_interpolation(x, y, xgrid(i), xgrid(i+1), ygrid(j), ygrid(j+1), &
              [var(ii, jj, kk+1), var(ii+1, jj, kk+1), var(ii, jj+1, kk+1), var(ii+1, jj+1, kk+1)])

  val = linear_interpolation(x, xgrid(i), xgrid(i+1), tmp(1, 1), tmp(1, 2))

  derivatives(1) = compute_derivative(x, xgrid(i+1), xgrid(i), val, tmp(1, 2), tmp(1, 1))
  derivatives(2) = compute_derivative(y, ygrid(j+1), ygrid(j), val, tmp(2, 2), tmp(2, 1))
  derivatives(3) = compute_derivative(z, zgrid(k+1), zgrid(k), val, tmp(3, 2), tmp(3, 1))

end subroutine compute_value_and_derivatives

end module biocfd_forcing

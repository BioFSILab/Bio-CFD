module biocfd_forcing
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use global, only : block, blk_start, nblocks
  use biocfd_interpolation, only: linear_interpolation, bilinear_interpolation
  use biocfd_blocks,only : Blocks
  implicit none

  private

  public :: pressureForcing1, pressureforcingfield, pressureforcingghost
  public :: velocityforcing1, velocityforcingfield, velocityforcingghost

  contains
SUBROUTINE pressureForcing1

      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1, g
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn, &
                   dpdn_e, ac_y, ac_z, at_y, at_z
      real(dp) :: derivatives(3)
      dpdn = 0._dp
        DO g=blk_start,nblocks


 !$acc parallel loop gang vector                                                                                          &
 !$acc private (n, n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc           dpdn_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc           i_z1,ac_z,ac_y,at_y,at_z)         &
 !$acc default(present)  &
 !$acc private(derivatives)
      DO n = 1, block(g)%ibCellCount
         IF (block(g)%ibSurfId(block(g)%nelp(n))==50) THEN
            !dpdn = 0.
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==51) THEN
            block(g)%thetaDot  = block(g)%thetaDot1
            block(g)%thetaDDot = block(g)%thetaDDot1
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==52) THEN
            block(g)%thetaDot  = block(g)%thetaDot2
            block(g)%thetaDDot = block(g)%thetaDDot2
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
        ENDIF
         dpdn = -((ac_z + at_z)*block(g)%cosGamma(block(g)%nelp(n)) + (ac_y + at_y) * &
                  block(g)%cosBeta(block(g)%nelp(n)))-block(g)%yddot * &
                     block(g)%cosBeta(block(g)%nelp(n))

         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

         sur2nodeDis = block(g)%pNormDis(n)

         pt1 = 1.5_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

        call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                           block(g)%xp, block(g)%yp, block(g)%zp, 0, &
                                           block(g)%p, p_pos1, derivatives)

        dpdn_e = derivatives(1) * block(g)%cosAlpha(block(g)%nelp(n)) &
               + derivatives(2) * block(g)%cosBeta(block(g)%nelp(n)) &
               + derivatives(3) * block(g)%cosGamma(block(g)%nelp(n))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*0.5_dp

         block(g)%p(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
 !$acc end parallel loop
      ENDDO

END SUBROUTINE pressureForcing1

SUBROUTINE velocityForcing1

      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1, g
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, &
                         usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
                         dudn_e, dvdn_e, dwdn_e
       real(dp) :: derivatives(3)

        DO g=blk_start, nblocks
 !$acc parallel loop gang vector         &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis,                       &
 !$acc          usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
 !$acc          dudn_e, dvdn_e, dwdn_e, &
 !$acc          k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc private(derivatives) &
 !$acc default(present)

      DO n = 1, block(g)%ibCellCount

         IF (block(g)%ibSurfID(block(g)%nelu2(n))==50) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu2(n))==51) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfId(block(g)%nelu2(n))==52) THEN
             usurf = 0._dp +block(g)%xdot
         ENDIF
         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

!***********************U(i,j,k)****************************************
         sur2nodeDis = block(g)%u2NormDis(n)

         pt1 = 1.5_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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
         IF(i_x1==block(g)%nx+2) i_x1 = block(g)%nx+1

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xu, block(g)%yu, block(g)%zu, 1, &
                                            block(g)%ut, u_pos1, derivatives)

         dudn_e =  derivatives(1) * block(g)%cosAlpha(block(g)%nelu2(n)) &
                 + derivatives(2) * block(g)%cosBeta(block(g)%nelu2(n)) &
                 + derivatives(3) * block(g)%cosGamma(block(g)%nelu2(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%ut(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval

!******************************U(i-1,j,k)*******************************
         IF (block(g)%ibSurfID(block(g)%nelu1(n))==50) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu1(n))==51) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfId(block(g)%nelu1(n))==52) THEN
             usurf = 0._dp + block(g)%xdot
         ENDIF
         sur2nodeDis = block(g)%u1NormDis(n)

         pt1 = 1.5_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

        call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xu, block(g)%yu, block(g)%zu, 1, &
                                            block(g)%ut, u_pos1, derivatives)

        dudn_e = derivatives(1) * block(g)%cosAlpha(block(g)%nelu1(n)) &
               + derivatives(2) * block(g)%cosBeta(block(g)%nelu1(n)) &
               + derivatives(3) * block(g)%cosGamma(block(g)%nelu1(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%ut(i-1,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j,k)*************************************
         IF (block(g)%ibSurfID(block(g)%nelv2(n))==50) THEN
                 vsurf = 0._dp + block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv2(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) &
                   - block(g)%piv_z) + block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv2(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) &
                   - block(g)%piv_z)+ block(g)%ydot  ! + ydot
         ENDIF

         sur2nodeDis = block(g)%v2NormDis(n)

         pt1 = 1.5_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 &
               + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         IF(i_y1==block(g)%ny+2) i_y1 = block(g)%ny+1

        call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xv, block(g)%yv, block(g)%zv, 2, &
                                            block(g)%vt, v_pos1, derivatives)
        dvdn_e = derivatives(1) * block(g)%cosAlpha(block(g)%nelv2(n)) &
                  + derivatives(2) * block(g)%cosBeta(block(g)%nelv2(n)) &
                  + derivatives(3) * block(g)%cosGamma(block(g)%nelv2(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%vt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j-1,k)*************************************
         IF (block(g)%ibSurfID(block(g)%nelv1(n))==50) THEN
           vsurf = 0._dp+ block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv1(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z) &
                      + block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv1(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z) &
                      + block(g)%ydot  ! + ydot
         ENDIF
         sur2nodeDis = block(g)%v1NormDis(n)

         pt1 = 1.5_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         IF(i_y1==1) i_y1 = 2

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xv, block(g)%yv, block(g)%zv, 2, &
                                            block(g)%vt, v_pos1, derivatives)

           dvdn_e = derivatives(1) * block(g)%cosAlpha(block(g)%nelv1(n)) &
                  + derivatives(2) * block(g)%cosBeta(block(g)%nelv1(n)) &
                  + derivatives(3) * block(g)%cosGamma(block(g)%nelv1(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%vt(i,j-1,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k)*************************************
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

         pt1 = 1.5_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         IF(i_z1==block(g)%nz+2) i_z1 = block(g)%nz+1

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xw, block(g)%yw, block(g)%zw, 3, &
                                            block(g)%wt, w_pos1, derivatives)
          dwdn_e = derivatives(1) *block(g)%cosAlpha(block(g)%nelw2(n)) &
                 + derivatives(2) *block(g)%cosBeta(block(g)%nelw2(n)) &
                 + derivatives(3) *block(g)%cosGamma(block(g)%nelw2(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%wt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k-1)*************************************
         wsurf = 0._dp
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

         pt1 = 1.5_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         IF(i_z1==1) i_z1 = 2

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xw, block(g)%yw, block(g)%zw, 3, &
                                            block(g)%wt, w_pos1, derivatives)

         dwdn_e =    derivatives(1) * block(g)%cosAlpha(block(g)%nelw1(n)) &
                   + derivatives(2) * block(g)%cosBeta(block(g)%nelw1(n)) &
                   + derivatives(3) * block(g)%cosGamma(block(g)%nelw1(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%wt(i,j,k-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
!$acc end parallel loop
      ENDDO

END SUBROUTINE velocityForcing1

SUBROUTINE pressureForcingGhost(blk)
      type(Blocks), intent(inout) :: blk
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

SUBROUTINE velocityForcingGhost

      INTEGER :: g,n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, &
                         usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
                         dudn_e, dvdn_e, dwdn_e
      real(dp) :: derivatives(3)

        DO g=blk_start,nblocks
 !$acc parallel loop gang vector         &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis, &
 !$acc          usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
 !$acc          dudn_e, dvdn_e, dwdn_e, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc default(present) private(derivatives)
      DO n = 1, block(g)%TSCellCount

        i = block(g)%TSIndexPtr(n, 1)
        j = block(g)%TSIndexPtr(n, 2)
        k = block(g)%TSIndexPtr(n, 3)
        !IF (block(g)%cell2(i,j,k).EQ.2) THEN
!***********************U(i,j,k)****************************************
         IF (block(g)%ibSurfID(block(g)%nelu2(block(g)%index_ts(n)))==50) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu2(block(g)%index_ts(n)))==51) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfId(block(g)%nelu2(block(g)%index_ts(n)))==52) THEN
           usurf = 0._dp + block(g)%xdot
         ENDIF

         sur2nodeDis = -block(g)%u2NormDis(block(g)%index_ts(n))

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 &
               + block(g)%deltay(j)**2 &
               + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i+1) - pt1*block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(n)))
         pos1_y = block(g)%yu(j) - pt1*block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(n)))
         pos1_z = block(g)%zu(k) - pt1*block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(n)))

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xu, block(g)%yu, block(g)%zu, &
                                            1, block(g)%ut, u_pos1, derivatives)

         dudn_e = -1 * (derivatives(1) *block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(n))) &
                  + derivatives(2) *block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(n))) &
                  + derivatives(3) *block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         block(g)%u2t_ghost(n) = block(g)%u(i,j,k)
!******************************U(i-1,j,k)*******************************
         IF (block(g)%ibSurfID(block(g)%nelu1(block(g)%index_ts(n)))==50) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu1(block(g)%index_ts(n)))==51) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfId(block(g)%nelu1(block(g)%index_ts(n)))==52) THEN
           usurf = 0._dp + block(g)%xdot
         ENDIF

         sur2nodeDis = -block(g)%u1NormDis(block(g)%index_ts(n))

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xu(i) - pt1*block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(n)))
         pos1_y = block(g)%yu(j) - pt1*block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(n)))
         pos1_z = block(g)%zu(k) - pt1*block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(n)))

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xu, block(g)%yu, block(g)%zu, &
                                            1, block(g)%ut, u_pos1, derivatives)

         dudn_e = -1 * (derivatives(1) * block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(n))) &
                      + derivatives(2) * block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(n))) &
                      + derivatives(3) * block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         block(g)%u1t_ghost(n) = block(g)%u(i-1,j,k)
!**************************V(i,j,k)*************************************
         IF (block(g)%ibSurfID(block(g)%nelv2(block(g)%index_ts(n)))==50) THEN
           vsurf = 0._dp+ block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv2(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(block(g)%index_ts(n))) &
                   - block(g)%piv_z) &
                   + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(block(g)%index_ts(n))) &
                   - block(g)%piv_x)+ block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv2(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(block(g)%index_ts(n))) &
                   - block(g)%piv_z) &
                   + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(block(g)%index_ts(n))) &
                   - block(g)%piv_x) + block(g)%ydot  ! + ydot
         ENDIF

         sur2nodeDis = -block(g)%v2NormDis(block(g)%index_ts(n))

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) - pt1*block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(n)))
         pos1_y = block(g)%yv(j+1) - pt1*block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(n)))
         pos1_z = block(g)%zv(k) - pt1*block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(n)))

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xv, block(g)%yv, block(g)%zv, &
                                            2, block(g)%vt, v_pos1, derivatives)
          dvdn_e = -1 * (derivatives(1) * block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(n))) &
                       + derivatives(2) * block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(n))) &
                       + derivatives(3) * block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%v2_ghost(n) = v(i,j,k)
         !ENDIF
         block(g)%v2t_ghost(n) = block(g)%v(i,j,k)
!**************************V(i,j-1,k)*************************************
         IF (block(g)%ibSurfID(block(g)%nelv1(block(g)%index_ts(n)))==50) THEN
           vsurf = 0._dp + block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv1(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(block(g)%index_ts(n))) &
                      - block(g)%piv_z)  &
                      + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(block(g)%index_ts(n)))&
                      - block(g)%piv_x) + block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv1(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(block(g)%index_ts(n))) &
                      - block(g)%piv_z)  &
                      + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(block(g)%index_ts(n)))&
                      - block(g)%piv_x)+ block(g)%ydot  ! + ydot
         ENDIF

         sur2nodeDis = -block(g)%v1NormDis(block(g)%index_ts(n))

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 &
               + block(g)%deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xv(i) - pt1*block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(n)))
         pos1_y = block(g)%yv(j) - pt1*block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(n)))
         pos1_z = block(g)%zv(k) - pt1*block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(n)))

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xv, block(g)%yv, block(g)%zv, &
                                            2, block(g)%vt, v_pos1, derivatives)
          dvdn_e = -1 * (derivatives(1) * block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(n))) &
                       + derivatives(2) * block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(n))) &
                       + derivatives(3) * block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         block(g)%v1t_ghost(n) = block(g)%v(i,j-1,k)
!**************************W(i,j,k)*************************************
         wsurf = 0._dp
         IF (block(g)%ibSurfID(block(g)%nelw2(block(g)%index_ts(n)))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw2(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           wsurf = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(block(g)%index_ts(n))) &
                   - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw2(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           wsurf = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw2(block(g)%index_ts(n))) &
                   - block(g)%piv_y)  ! + ydot
         ENDIF


         sur2nodeDis = -block(g)%w2NormDis(block(g)%index_ts(n))

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) - pt1*block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(n)))
         pos1_y = block(g)%yw(j) - pt1*block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(n)))
         pos1_z = block(g)%zw(k+1) - pt1*block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(n)))

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xw, block(g)%yw, block(g)%zw, &
                                            3, block(g)%wt, w_pos1, derivatives)
         dwdn_e = -1 * (derivatives(1) * block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(n))) &
                      + derivatives(2) * block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(n))) &
                      + derivatives(3) * block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w2_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !ELSE
         !block(g)%w2_ghost(n) = w(i,j,k)
         !ENDIF
         block(g)%w2t_ghost(n) = block(g)%w(i,j,k)
!**************************W(i,j,k-1)*************************************
         wsurf = 0._dp
         IF (block(g)%ibSurfID(block(g)%nelw1(block(g)%index_ts(n)))==50) THEN
           wsurf = 0.
         ELSEIF (block(g)%ibSurfID(block(g)%nelw1(block(g)%index_ts(n)))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           wsurf = block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(block(g)%index_ts(n))) &
                   - block(g)%piv_y)
         ELSEIF (block(g)%ibSurfId(block(g)%nelw1(block(g)%index_ts(n)))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           wsurf =  block(g)%thetaDot*(block(g)%ycent(block(g)%nelw1(block(g)%index_ts(n))) &
                   - block(g)%piv_y)  ! + ydot
         ENDIF


         sur2nodeDis = -block(g)%w1NormDis(block(g)%index_ts(n))

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

         !coordinates of three points from interceptd cell pressure node
         pos1_x = block(g)%xw(i) - pt1*block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(n)))
         pos1_y = block(g)%yw(j) - pt1*block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(n)))
         pos1_z = block(g)%zw(k) - pt1*block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(n)))

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
         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xw, block(g)%yw, block(g)%zw, &
                                            3, block(g)%wt, w_pos1, derivatives)
         dwdn_e = -1 * (derivatives(1) * block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(n))) &
                      + derivatives(2) * block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(n))) &
                      + derivatives(3) * block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(n))))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w1_ghost(n) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         block(g)%w1t_ghost(n) = block(g)%w(i,j,k-1)
      ENDDO
      !$acc end parallel loop
      ENDDO

END SUBROUTINE velocityForcingGhost
!***********************************************************************

SUBROUTINE pressureForcingField

      INTEGER :: g,n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1, &
                   aval, bval, cval, p_pos1, sur2nodeDis, dpdn, &
                   dpdn_e, ac_y, ac_z, at_y, at_z
       real(dp) :: derivatives(3)

       DO g=blk_start,nblocks
      dpdn = 0._dp
 !$acc parallel loop gang vector                                                                                          &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
 !$acc          dpdn_e, k, j, i, il, jl, kl, i_x1, i_y1,               &
 !$acc          i_z1,ac_z,ac_y,at_y,at_z)         &
 !$acc default(present) private(derivatives)
      DO n = 1, block(g)%ibCellCount
         IF (block(g)%ibSurfId(block(g)%nelp(n))==50) THEN
            ac_z = 0.  !-block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_z)
            ac_y = 0.  !-block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_z = 0.  ! block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(block(g)%index_ts(n))) - block(g)%piv_y)
            at_y = 0.  !
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==51) THEN
            block(g)%thetaDot  = block(g)%thetaDot1
            block(g)%thetaDDot = block(g)%thetaDDot1
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
        ELSEIF (block(g)%ibSurfId(block(g)%nelp(n))==52) THEN
            block(g)%thetaDot  = block(g)%thetaDot2
            block(g)%thetaDDot = block(g)%thetaDDot2
            ac_z = -block(g)%thetaDot**2*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
            ac_y = -block(g)%thetaDot**2*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_z =  block(g)%thetaDDot*(block(g)%ycent(block(g)%nelp(n)) - block(g)%piv_y)
            at_y = -block(g)%thetaDDot*(block(g)%zcent(block(g)%nelp(n)) - block(g)%piv_z)
        ENDIF
         dpdn = -((ac_z + at_z)*block(g)%cosGamma(block(g)%nelp(n)) &
                + (ac_y + at_y)*block(g)%cosBeta(block(g)%nelp(n))) &
                -block(g)%yddot*block(g)%cosBeta(block(g)%nelp(n))

         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

         sur2nodeDis = block(g)%pNormDis(n)

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                           block(g)%xp, block(g)%yp, block(g)%zp, 0, &
                                           block(g)%p, p_pos1, derivatives)

         dpdn_e  = derivatives(1) * block(g)%cosAlpha(block(g)%nelp(n)) &
                 + derivatives(2) * block(g)%cosBeta(block(g)%nelp(n))  &
                 + derivatives(3) * block(g)%cosGamma(block(g)%nelp(n))

         n1 = pt1 + sur2nodeDis

         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*0.5_dp

         block(g)%p(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
      !$acc end parallel loop
      ENDDO
END SUBROUTINE pressureForcingField

SUBROUTINE velocityForcingField

      INTEGER :: g,n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (dp) :: n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, &
                         usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
                         dudn_e, dvdn_e, dwdn_e
      real(dp) :: derivatives(3)

        DO g=blk_start,nblocks
 !$acc parallel loop gang vector         &
 !$acc private (n1, pos1_x, pos1_y, pos1_z, pt1,           &
 !$acc          aval, bval, cval, sur2nodeDis, &
 !$acc          usurf, u_pos1, vsurf, v_pos1, wsurf, w_pos1, &
 !$acc          dudn_e, dvdn_e, dwdn_e, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &
 !$acc default(present) private(derivatives)
      DO n = 1, block(g)%ibCellCount

         i = block(g)%interceptedIndexPtr(n, 1)
         j = block(g)%interceptedIndexPtr(n, 2)
         k = block(g)%interceptedIndexPtr(n, 3)

!***********************U(i,j,k)****************************************
         IF (block(g)%ibSurfID(block(g)%nelu2(n))==50) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu2(n))==51) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfId(block(g)%nelu2(n))==52) THEN
             usurf = 0._dp + block(g)%xdot
         ENDIF

         sur2nodeDis = block(g)%u2NormDis(n)

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xu, block(g)%yu, block(g)%zu, &
                                            1, block(g)%ut, u_pos1, derivatives)

         dudn_e =   derivatives(1) * block(g)%cosAlpha(block(g)%nelu2(n)) &
                  + derivatives(2) * block(g)%cosBeta(block(g)%nelu2(n)) &
                  + derivatives(3) * block(g)%cosGamma(block(g)%nelu2(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!******************************U(i-1,j,k)*******************************
         IF (block(g)%ibSurfID(block(g)%nelu1(n))==50) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfID(block(g)%nelu1(n))==51) THEN
             usurf = 0._dp + block(g)%xdot
         ELSEIF (block(g)%ibSurfId(block(g)%nelu1(n))==52) THEN
             usurf = 0._dp + block(g)%xdot
         ENDIF

         sur2nodeDis = block(g)%u1NormDis(n)

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 &
               + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2) &
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xu, block(g)%yu, block(g)%zu, &
                                            1, block(g)%ut, u_pos1, derivatives)

         dudn_e =   derivatives(1) * block(g)%cosAlpha(block(g)%nelu1(n)) &
                  + derivatives(2) * block(g)%cosBeta(block(g)%nelu1(n)) &
                  + derivatives(3) * block(g)%cosGamma(block(g)%nelu1(n))

         n1 = pt1 + sur2nodeDis

         cval = usurf
         bval = 2._dp/n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         block(g)%u(i-1,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j,k)*************************************
         IF (block(g)%ibSurfID(block(g)%nelv2(n))==50) THEN
           vsurf = 0._dp + block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv2(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z) &
                   + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(n)) - block(g)%piv_x) &
                   + block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv2(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv2(n)) - block(g)%piv_z) &
                   + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv2(n)) - block(g)%piv_x) &
                   + block(g)%ydot
         ENDIF

         sur2nodeDis = block(g)%v2NormDis(n)

         pt1 = 1.21_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xv, block(g)%yv, block(g)%zv, &
                                            2, block(g)%vt, v_pos1, derivatives)
         dvdn_e =  derivatives(1) * block(g)%cosAlpha(block(g)%nelv2(n)) &
                 + derivatives(2) * block(g)%cosBeta(block(g)%nelv2(n)) &
                 + derivatives(3) * block(g)%cosGamma(block(g)%nelv2(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************V(i,j-1,k)*************************************
         IF (block(g)%ibSurfID(block(g)%nelv1(n))==50) THEN
           vsurf = 0._dp + block(g)%ydot
         ELSEIF (block(g)%ibSurfID(block(g)%nelv1(n))==51) THEN
           block(g)%thetaDot =  block(g)%thetaDot1
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)  &
                   + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(n)) - block(g)%piv_x) &
                   + block(g)%ydot
         ELSEIF (block(g)%ibSurfId(block(g)%nelv1(n))==52) THEN
           block(g)%thetaDot =  block(g)%thetaDot2
           !vsurf    = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)! + ydot
           vsurf = -block(g)%thetaDot*(block(g)%zcent(block(g)%nelv1(n)) - block(g)%piv_z)  &
                   + (-block(g)%alphaDot)*( block(g)%xcent(block(g)%nelv1(n)) - block(g)%piv_x) &
                   + block(g)%ydot
         ENDIF

         !vsurf =  -block(g)%thetaDot*(block(g)%xcent(block(g)%nelv1(n)) - block(g)%piv_x) + ydot
         sur2nodeDis = block(g)%v1NormDis(n)

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
              + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xv, block(g)%yv, block(g)%zv, &
                                            2, block(g)%vt, v_pos1, derivatives)
         dvdn_e = derivatives(1) * block(g)%cosAlpha(block(g)%nelv1(n)) &
                + derivatives(2) * block(g)%cosBeta(block(g)%nelv1(n)) &
                + derivatives(3) * block(g)%cosGamma(block(g)%nelv1(n))

         n1 = pt1 + sur2nodeDis

         cval = vsurf
         bval = 2._dp/n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         block(g)%v(i,j-1,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k)*************************************
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
         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xw, block(g)%yw, block(g)%zw, &
                                            3, block(g)%wt, w_pos1, derivatives)
         dwdn_e = derivatives(1) * block(g)%cosAlpha(block(g)%nelw2(n)) &
                + derivatives(2) * block(g)%cosBeta(block(g)%nelw2(n)) &
                + derivatives(3) * block(g)%cosGamma(block(g)%nelw2(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
!**************************W(i,j,k-1)*************************************
         wsurf = 0._dp
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

         pt1 = 1.51_dp*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)&
               + (dabs(sur2nodeDis)-sur2nodeDis)*0.5_dp

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

         call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
                                            block(g)%xw, block(g)%yw, block(g)%zw, &
                                            3, block(g)%wt, w_pos1, derivatives)
         dwdn_e = derivatives(1) * block(g)%cosAlpha(block(g)%nelw1(n)) &
                + derivatives(2) * block(g)%cosBeta(block(g)%nelw1(n)) &
                + derivatives(3) * block(g)%cosGamma(block(g)%nelw1(n))

         n1 = pt1 + sur2nodeDis

         cval = wsurf
         bval = 2._dp/n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         block(g)%w(i,j,k-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
      ENDDO
      !$acc end parallel loop
      ENDDO

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

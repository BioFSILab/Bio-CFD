!***********************************************************************
        SUBROUTINE findDistnode
        USE global
        IMPLICIT NONE
        REAL (KIND=8)      ::  dist, dist1, dist2
        INTEGER (kind = 8) ::  i, g


        DO g=blk_start,nblocks
        dist=0.
        dist1=999999.
        dist2=0.
        DO i = 1, block(g)%ibnodes
        IF (block(g)%ibNodeId(i)==51) THEN
                if (abs(block(g)%znode(i)) > dist)then
                        dist=block(g)%znode(i)
                        block(g)%mk=i
                endif
        ENDIF
        IF (block(g)%ibNodeId(i)==51) THEN
                if (abs(block(g)%xnode(i)) < dist1)then
                        dist1=block(g)%xnode(i)
                        block(g)%mkx1=i
                endif
                if (abs(block(g)%xnode(i)) > dist2)then
                        dist2=block(g)%xnode(i)
                        block(g)%mkx2=i
                endif
        ENDIF
        ENDDO
        ENDDO
        END SUBROUTINE

!***********************************************************************
     SUBROUTINE shiftSurfaceNodesInitial
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  i, g
        REAL (KIND=8)      ::  xr, yr, zr
        REAL (KIND=8)      ::  xr1, yr1, zr1, angt
        REAL (KIND=8)      :: bdy,bdfr, angg

        phase_angle = phase_angle*pi/180_rk
        aoa1 = aoa*pi/180_rk
        aoa2 = -aoa1
        alpha_m = alpha_m*pi/180_rk
        theta_m = theta_m*pi/180_rk
        !freq!/(2.*a0y)
        a0y = 0.  !a0y
        ang_theta = 0.  !2._rk*pi*freq
        alpha_t=(alpha_m*0.5)*(1+cos(ang_theta*(totime+deltat)+phase_angle))
        theta_t       =  theta_m*cos(ang_theta*(totime+deltat))
        DO g=blk_start,nblocks
        ALLOCATE ( block(g)%xnode1(block(g)%ibNodes), block(g)%ynode1(block(g)%ibNodes), block(g)%znode1(block(g)%ibNodes) )
        block(g)%xnode1 = block(g)%xnode
        block(g)%ynode1 = block(g)%ynode
        block(g)%znode1 = block(g)%znode
        xfact=0.05
        bdfr=15
        bdy=15*dxmin
        angt  =  2._rk*pi*bdfr

        ac_x_al=0.
        ac_y_al=0.
        at_x_al=0.
        at_y_al=0.
        ac_x=0.
        ac_y=0.
        ac_z=0.
	  !xnode1 = xnode
         !ynode1 = ynode
         !znode1 = znode

        !IF(ita.eq.0) THEN !!For ita equal to 0
	  block(g)%u_init = 0.
	  block(g)%u_final = 0.
	  block(g)%v_init = 0.
	  block(g)%v_final = 0.
	  block(g)%w_init = 0.
	  block(g)%w_final = 0.
	  block(g)%xmove = 0.
	  block(g)%ymove = 0.
	  block(g)%zmove = 0.
	 !END IF

        block(g)%a0 = block(g)%a0*pi/180_rk
         angt  =  2._rk*pi*block(g)%bfreq
         block(g)%xpth1=block(g)%xshift-(ita*dxmin*xfact)
        block(g)%xpth2=block(g)%xshift-(ita*dxmin*xfact)
        block(g)%ypth1=(block(g)%yamp)*sin(angt*block(g)%xshift)
        block(g)%ypth2=(block(g)%yamp)*sin(angt*block(g)%xshift)
        !block(g)%piv_x = block(g)%xshift + piv_pt
        block(g)%piv_x = block(g)%xshift
        block(g)%piv_y = block(g)%yshift
        block(g)%piv_z = block(g)%zshift
        !piv_y = block(g)%yshift
       !block(g)%z_p=0.;block(g)%y_p=0.; block(g)%v_y_p=0.; block(g)%a_y_p=0.
       !thetaDot  =  0._rk
       !thetaDDot =  0._rk
       !xt =  0._rk
       !xdot = 0._rk
       !xddot = 0._rk

       !yt =  0._rk
       !ydot = 0._rk
       !yddot = 0._rk

        block(g)%thetaDot  =  0.
        block(g)% thetaDDot =  0._rk
        block(g)% thetaDot   =  0.  !ang_theta*a0*cos(2._rk*pi*freq*totime + phase_angle)
        block(g)%thetaDot1  = 0.
        block(g)% thetaDot2  = 0.
        block(g)% alphaDot  = 0.
        block(g)%alphaDDot  = 0.  !-ang_theta*ang_theta*a0*sin(2._rk*pi*freq*totime + phase_angle)
        block(g)% thetaDDot1 = 0.
        block(g)%thetaDDot2 = 0.
        block(g)% thetaDDot  = 0.  !-ang_theta*ang_theta*a0*sin(2._rk*pi*freq*totime + phase_angle)
       block(g)%yt         =  bdy*sin(2*pi*bdfr*totime )
       block(g)%ydot       =  angt*bdy*cos(2*pi*bdfr*totime)
       block(g)%yddot      =  -angt*angt*bdy*sin(2*pi*bdfr*totime)
       block(g)%xt         =  block(g)%xshift- (ita*dxmin*xfact)
       block(g)%xdot       =  -(dxmin*xfact)/deltat
       block(g)%xddot      =  0.
        block(g)%u_prev = 0._rk
        block(g)%u_curr = 0._rk
        block(g)%v_prev = 0._rk
        block(g)%v_curr = 0._rk
        block(g)%w_prev = 0._rk
        block(g)%w_curr = 0._rk
        block(g)%Total_VP_FY = 0._rk
        block(g)%Total_VP_FX = 0._rk
        block(g)%inity_cent=block(g)%yshift
        block(g)%nxty_cent=block(g)%yshift
        block(g)%initx_cent=block(g)%xshift
        block(g)%nxtx_cent=block(g)%xshift
        block(g)%initz_cent=block(g)%zshift
        block(g)%nxtz_cent=block(g)%zshift
        block(g)%ypos=block(g)%yshift
        block(g)%xpos=block(g)%xshift
       !!$acc parallel loop gang vector present (xnode1, ynode1, znode1, xnode, ynode, znode) &
       !!$acc private (xr, yr, zr) firstprivate (xshift, block(g)%yshift, block(g)%zshift)
        DO i = 1, block(g)%ibnodes
        IF (block(g)%ibNodeId(i)==51) THEN
            xr1 =  block(g)%xnode(i)
          !!zr1 =  zShift +(block(g)%ynode(i)-yShift)*sin(aoa1) + (block(g)%znode(i)-zShift)*cos(aoa1)
          !!yr1 =  yShift +(block(g)%ynode(i)-yShift)*cos(aoa1) - (block(g)%znode(i)-zShift)*sin(aoa1)
            zr1 =  block(g)%znode(i)*cos(aoa1) + block(g)%ynode(i)*sin(aoa1) + piv_pt - piv_pt*cos(aoa1)
            yr1 = -block(g)%znode(i)*sin(aoa1) + block(g)%ynode(i)*cos(aoa1) + piv_pt*sin(aoa1)
            !!xr1 =  block(g)%xnode(i)
            !!zr1 = block(g)% znode(i)*cos(theta_t) +block(g)% ynode(i)*sin(theta_t) + piv_pt - piv_pt*cos(theta_t)
            !!yr1 = -block(g)%znode(i)*sin(theta_t) +block(g)% ynode(i)*cos(theta_t) + piv_pt*sin(theta_t)
         ELSEIF (block(g)%ibNodeId(i)==52) THEN
             !xr1 = block(g)% xnode(i)
             !zr1 = block(g)% znode(i)*cos(-theta_t) +block(g)% ynode(i)*sin(-theta_t) + piv_pt - piv_pt*cos(-theta_t)
             !yr1 = -block(g)%znode(i)*sin(-theta_t) +block(g)% ynode(i)*cos(-theta_t) + piv_pt*sin(-theta_t)
            xr1 =  block(g)%xnode(i)
           !zr1 =  zShift +(block(g)%ynode(i)-yShift)*sin(aoa2) + (block(g)%znode(i)-zShift)*cos(aoa2)
           !yr1 =  yShift +(block(g)%ynode(i)-yShift)*cos(aoa2) - (block(g)%znode(i)-zShift)*sin(aoa2)
            zr1 =  block(g)%znode(i)*cos(aoa2) + block(g)%ynode(i)*sin(aoa2)  + piv_pt - piv_pt*cos(aoa2)
            yr1 = -block(g)%znode(i)*sin(aoa2) + block(g)%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
             !xr =  xnode(i)*cos(aoa1) + ynode(i)*sin(aoa1)
             !yr = -xnode(i)*sin(aoa1) + ynode(i)*cos(aoa1)
             !xnode1(i) = xr + piv_pt - piv_pt*cos(aoa1) + xshift
             !ynode1(i) = yr + piv_pt*sin(aoa1) + yt + block(g)%yshift
             !znode1(i) = znode(i) + block(g)%zshift
           ELSE
               xr1 = block(g)% xnode(i)
               zr1 = block(g)% znode(i)  !*cos(aoa1) + ynode(i)*sin(aoa1)
               yr1 = block(g)% ynode(i)  !*sin(aoa1) + ynode(i)*cos(aoa1)
             !xr =  xnode(i)*cos(aoa1) + ynode(i)*sin(aoa1)
             !yr = -xnode(i)*sin(aoa1) + ynode(i)*cos(aoa1)
             !xnode1(i) = xr + piv_pt - piv_pt*cos(aoa1) + xshift
             !ynode1(i) = yr + piv_pt*sin(aoa1) + yt + block(g)%yshift
             !znode1(i) = znode(i) + block(g)%zshift
           ENDIF
           block(g)%xnode1(i) = xr1+ block(g)%xShift
               block(g)%ynode1(i) = yr1+ block(g)%yShift
           block(g)%znode1(i) = zr1+ block(g)%zshift
        ENDDO
       !!$acc end parallel

        END DO

       !!$acc update host (xnode1, ynode1, znode1)
      END SUBROUTINE shiftSurfaceNodesInitial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE computeSurfaceVariables
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  i, g
        REAL (KIND=8)      ::  xr, yr, zr
        REAL (KIND=8)      ::  xr1, yr1, zr1
        REAL (KIND=8)      :: vol, angg, angt
        REAL (KIND=8)      :: bdy,bdfr, yt_prv
        CHARACTER(len=150) :: filename1

        DO g=blk_start,nblocks
        angg=90
        !aoa1       =  (45*pi/180_rk)*sin(2._rk*pi*freq*(totime+deltat) + phase_angle)
        aoa1       =  (block(g)%a0)*sin(2._rk*pi*freq*(totime+deltat) + phase_angle)
        aoa2       = -aoa1
        ang_theta  =  2._rk*pi*freq
        bdfr=block(g)%bfreq
        bdy=block(g)%yamp
        !bdy=a0
        angt  =  2._rk*pi*bdfr
        !angt  =  2._rk*pi*freq

        !alpha_t=(alpha_m*0.5)*(1+cos((ang_theta*(totime+deltat))+phase_angle))
        !aoa1       =  a0*sin(2._rk*pi*freq*(totime+deltat) + phase_angle)
        !theta_t       =  theta_m*cos(ang_theta*(totime+deltat))
        !aoa2       = -aoa1
        !PRINT*, "angles =", aoa1*180./pi, aoa2*180./pi
        !PRINT*,"alpha =", alpha_t*180./pi,"theta=",theta_t*180./pi


        block(g)%xpth2=block(g)%xshift-(ita*dxmin*xfact)
        block(g)%ypth2=(block(g)%yamp)*sin(angt*block(g)%xpth2)
!        block(g)%ychg=block(g)%ypth2-block(g)%ypth1
        block(g)%xchg=block(g)%xpth2-block(g)%xpth1
        block(g)%ypth1=block(g)%ypth2
        block(g)%xpth1=block(g)%xpth2
        !PRINT*, a0, sin(2._rk*pi*freq*(totime+deltat) + phase_angle) ,aoa1
        !PRINT*, freq, totime+deltat, phase_angle, pi
        PRINT*, "angles =", aoa1*180./pi, aoa2*180./pi
        !PAUSE
        block(g)%thetaDot1  =  ang_theta*block(g)%a0*cos(2._rk*pi*freq*(totime+deltat) + phase_angle)
            block(g)%thetaDDot1 = -ang_theta*ang_theta*block(g)%a0*sin(2._rk*pi*freq*(totime+deltat) + phase_angle)
        block(g)%thetaDot2  = -block(g)%thetaDot1
        block(g)%thetaDDot2 = -block(g)%thetaDDot1

                !block(g)% thetaDot1  =  -theta_m*ang_theta*sin(ang_theta*(totime+deltat))
        !block(g)% thetaDDot1 = -ang_theta*ang_theta*theta_m*cos(ang_theta*freq*(totime+deltat))
        !    block(g)% thetaDot2  = -block(g)%thetaDot1
        !        block(g)% thetaDDot2 = -block(g)%thetaDDot1
        !block(g)%alphaDot  = (alpha_m*0.5)*(-ang_theta*sin((ang_theta*(totime+deltat))+phase_angle))
       !block(g)% alphaDDot  = (alpha_m*0.5)*(-ang_theta*ang_theta*cos((ang_theta*(totime+deltat))+phase_angle))

       block(g)%yt         =  bdy*sin(angt*(totime) )
       block(g)%ydot       =  angt*bdy*cos(angt*(totime))
       block(g)%yddot      =  -angt*angt*bdy*sin(angt*(totime))

       block(g)%xt         = block(g)%xpth2
       block(g)%xdot       = -(dxmin*xfact)/deltat
       block(g)%xddot      = 0.

           block(g)%ychg=block(g)%yt - bdy*sin(angt*(totime-deltat) )
        block(g)%ymove = block(g)%yt
        block(g)%xmove = block(g)%xchg
        block(g)%zmove = 0.

        block(g)%ypos =  block(g)%ypos  + block(g)%ychg
        block(g)%xpos =  block(g)%xpos  + block(g)%xmove
        block(g)%piv_y = block(g)%piv_y + block(g)%ychg
        block(g)%piv_x = block(g)%piv_x + block(g)%xmove
        block(g)%piv_z = block(g)%piv_z + block(g)%zmove
        block(g)%nxty_cent= block(g)%nxty_cent + block(g)%ychg
        block(g)%nxtx_cent= block(g)%nxtx_cent + block(g)%xmove
     !   print*, block(g)%ymove, block(g)%v_prev, block(g)%v_curr, block(g)%accnp_Y, block(g)%Total_FY
        WRITE(filename1,1) block(g)%fineg,re, g
      1  FORMAT('d',I4.4,'_index.',F8.2,'.',i3.1,".dat")
        OPEN(UNIT = 17, FILE = filename1,Access='Append', STATUS = 'unknown')
        !open(17,file='index.dat',access='append',status='unknown')
        !write(17,14) totime, block(g)%nxtx_cent, block(g)%xmove , block(g)%xchg, block(g)%xdot, block(g)%xpth2
        write(17,14) totime, block(g)%nxty_cent, block(g)%inity_cent, block(g)%ymove, block(g)%nxtx_cent, block(g)%xmove
!!        write(17,14) totime, block(g)%xpos, block(g)%u_curr,block(g)%Total_VP_FX,block(g)%Total_FX,block(g)%accnp_X
     close(17)
     14      FORMAT(7F15.8)
        print*,'centn',block(g)%nxty_cent,'centi',block(g)%inity_cent,'mv',block(g)%ychg
       !!$acc parallel loop gang vector default(present) private(xr1,yr1,zr1) firstprivate(aoa1,aoa,piv_pt,g, xshift, block(g)%yshift, block(g)%zshift)
      DO i = 1, block(g)%ibnodes
      IF (block(g)%ibNodeId(i)==51) THEN
             xr1 =  block(g)%xnode(i)
             zr1 =  block(g)%znode(i)*cos(aoa1) + block(g)%ynode(i)*sin(aoa1) + piv_pt - piv_pt*cos(aoa1)
             yr1 = -block(g)%znode(i)*sin(aoa1) + block(g)%ynode(i)*cos(aoa1) + piv_pt*sin(aoa1)
            !!xr1 = block(g)% xnode(i)*cos(alpha_t) -block(g)% ynode(i)*sin(alpha_t)*cos(theta_t) +block(g)% znode(i)*sin(alpha_t)*sin(theta_t)
            !!yr1 = block(g)% xnode(i)*sin(alpha_t) +block(g)% ynode(i)*cos(alpha_t)*cos(theta_t) -block(g)% znode(i)*sin(theta_t)*cos(alpha_t)
            !!zr1 = block(g)% znode(i)*cos(theta_t) +block(g)% ynode(i)*sin(theta_t)
             !xr =  xnode(i)*cos(aoa1) + ynode(i)*sin(aoa1)
             !yr = -xnode(i)*sin(aoa1) + ynode(i)*cos(aoa1)
             !xnode1(i) = xr + piv_pt - piv_pt*cos(aoa1) + xshift
             !ynode1(i) = yr + piv_pt*sin(aoa1) + yt + block(g)%yshift
             !znode1(i) = znode(i) + block(g)%zshift
         ELSEIF (block(g)%ibNodeId(i)==52) THEN
             xr1 =  block(g)%xnode(i)
             zr1 =  block(g)%znode(i)*cos(aoa2) + block(g)%ynode(i)*sin(aoa2)  + piv_pt - piv_pt*cos(aoa2)
             yr1 = -block(g)%znode(i)*sin(aoa2) + block(g)%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
            !!xr1 = block(g)% xnode(i)*cos(alpha_t) -block(g)% ynode(i)*sin(alpha_t)*cos(-theta_t) +block(g)% znode(i)*sin(alpha_t)*sin(-theta_t)
            !!yr1 = block(g)% xnode(i)*sin(alpha_t) +block(g)% ynode(i)*cos(alpha_t)*cos(-theta_t) -block(g)% znode(i)*sin(-theta_t)*cos(alpha_t)
            !!zr1 = block(g)% znode(i)*cos(-theta_t) +block(g)% ynode(i)*sin(-theta_t)
             !xr1 =  xnode(i)
             !zr1 =  znode(i)*cos(aoa2) + ynode(i)*sin(aoa2)  + piv_pt - piv_pt*cos(aoa2)
             !yr1 = -znode(i)*sin(aoa2) + ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)

             !xr =  xnode(i)*cos(aoa1) + ynode(i)*sin(aoa1)
             !yr = -xnode(i)*sin(aoa1) + ynode(i)*cos(aoa1)
             !xnode1(i) = xr + piv_pt - piv_pt*cos(aoa1) + xshift
             !ynode1(i) = yr + piv_pt*sin(aoa1) + yt + block(g)%yshift
             !znode1(i) = znode(i) + block(g)%zshift
           ELSE
               xr1 = block(g)% xnode(i)
               zr1 = block(g)% znode(i)
               yr1 = block(g)% ynode(i)
             !xr =  xnode(i)*cos(aoa1) + ynode(i)*sin(aoa1)
             !yr = -xnode(i)*sin(aoa1) + ynode(i)*cos(aoa1)
             !xnode1(i) = xr + piv_pt - piv_pt*cos(aoa1) + xshift
             !ynode1(i) = yr + piv_pt*sin(aoa1) + yt + block(g)%yshift
             !znode1(i) = znode(i) + block(g)%zshift
           ENDIF

           !block(g)% xnode1(i) = xr1 +xshift +block(g)%xmove
           block(g)% xnode1(i) = xr1 +block(g)%xpth2
           block(g)% ynode1(i) = yr1 +block(g)%yshift +block(g)%ymove
           block(g)% znode1(i) = zr1 +block(g)%zshift +block(g)%zmove

         !IF (block(g)%ibNodeId(i).ne.52 .and. block(g)%ibNodeId(i).ne.51 ) THEN
         !  print*,yr1,block(g)%yshift,block(g)%ymove
         !endif
        ENDDO
        !!$acc end parallel
           write(*,*) block(g)%xnode1(1),block(g)% ynode1(1), block(g)% znode1(1)
        ENDDO

       !!$acc update host (xnode1, ynode1, znode1)
      END SUBROUTINE computeSurfaceVariables
!***********************************************************************

!***********************************************************************
      SUBROUTINE computeSurfaceNorm
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, g, nv  !c1, c2, c3, c4
        REAL (KIND=8)      :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL, binor
        REAL (KIND=8)      :: var_xcent, var_ycent, var_zcent, xlim2,zlim1,xlim1, ylim1, ylim2




        DO g=blk_start, nblocks

        !block(g)%ibElemCnt=0.
        ALLOCATE (block(g)%xcent(block(g)%ibElems), block(g)%ycent(block(g)%ibElems), block(g)%zcent(block(g)%ibElems), &
                  block(g)%cosAlpha(block(g)%ibElems), block(g)%cosBeta(block(g)%ibElems), block(g)%cosGamma(block(g)%ibElems), &
                  block(g)%alpha3(block(g)%ibElems), block(g)%beta3(block(g)%ibElems), block(g)%gamma3(block(g)%ibElems))

       ! END DO
        !inor = -1._rk
        !compute centroid and direction cosines
       !!$acc parallel loop gang vector       &
       !!$acc present (xnode1, ynode1, znode1, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, alpha3, beta3, gamma3, ibElP1, ibElP2, ibElP3)   &
       !!$acc private (p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL)  firstprivate (inor)
       !$acc parallel loop gang vector default(present) private (var_xcent, var_ycent, var_zcent,p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL)  firstprivate (inor)
        DO n = 1, block(g)%ibElems
           p1x = block(g)%xnode1(block(g)%ibElP1(n))                       !x coordinate element node 1
           p1y = block(g)%ynode1(block(g)%ibElP1(n))                       !y coordinate element node 1
           p1z = block(g)%znode1(block(g)%ibElP1(n))                       !z coordinate element node 1

           p2x = block(g)%xnode1(block(g)%ibElP2(n))                       !x coordinate element node 2
           p2y = block(g)%ynode1(block(g)%ibElP2(n))                       !y coordinate element node 2
           p2z = block(g)%znode1(block(g)%ibElP2(n))                       !z coordinate element node 2

           p3x = block(g)%xnode1(block(g)%ibElP3(n))                       !x coordinate element node 3
           p3y = block(g)%ynode1(block(g)%ibElP3(n))                       !y coordinate element node 3
           p3z = block(g)%znode1(block(g)%ibElP3(n))                       !z coordinate element node 3


           var_xcent =  (p2x+p1x+p3x)/3._rk                  !centroid x coordinate element
           var_ycent =  (p2y+p1y+p3y)/3._rk                  !centroid y coordinate element
           var_zcent =  (p2z+p1z+p3z)/3._rk                  !centroid z coordinate element


      !   if ( var_ycent.le. ylim1 .and. var_xcent .ge. xlim1 )then
      !                 g=4
      !  elseif( var_ycent .ge. ylim1 .and. var_xcent .ge. xlim2) then
      !                 g=3

      !  elseif(var_xcent .le. xlim2 .and. var_ycent .ge. ylim2)then

      !                 g=2
      ! elseif(var_xcent .le. xlim1 .and. var_ycent .le. ylim2)then
      !
      !                 g=1
      !
      ! end if
      !   if ( var_ycent.le. ylim1 )then
      !                 g=2
      !   else
      !                 g=1
      !   endif
      !    block(g)%ibElemCnt= block(g)%ibElemCnt+1
           !block(g)%ibElemNum( block(g)%ibElemCnt)=n

      !     nv=block(g)%ibElemCnt
           block(g)%xcent(n) =  var_xcent                  !centroid x coordinate element
           block(g)%ycent(n) =  var_ycent                  !centroid y coordinate element
           block(g)%zcent(n) =  var_zcent                  !centroid z coordinate element

           block(g)%cosAlpha(n) = (p2y-p1y)*(p3z-p1z)-(p3y-p1y)*(p2z-p1z)
           block(g)%cosBeta(n)  = (p2z-p1z)*(p3x-p1x)-(p3z-p1z)*(p2x-p1x)
           block(g)%cosGamma(n) = (p2x-p1x)*(p3y-p1y)-(p3x-p1x)*(p2y-p1y)

           block(g)%alpha3(n) = block(g)%cosAlpha(n)
           block(g)%beta3(n)  = block(g)%cosBeta(n)
           block(g)%gamma3(n) = block(g)%cosGamma(n)

           lenEL       = dsqrt(block(g)%cosAlpha(n)**2 + block(g)%cosBeta(n)**2 + block(g)%cosGamma(n)**2)   !length of element

                binor=inor

           block(g)%cosAlpha(n) = block(g)%cosAlpha(n)/lenEl*binor               !direction cosine unit normal along x
           block(g)%cosBeta(n)  = block(g)%cosBeta(n)/lenEl*binor                !direction cosine unit normal along y
           block(g)%cosGamma(n) = block(g)%cosGamma(n)/lenEl*binor               !direction cosine unit normal along z
        ENDDO
       !$acc end parallel

       !!$acc update host (xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, alpha3, beta3, gamma3)
        print*, 'SurfaceNorm done, inor =', inor

        END DO


     END SUBROUTINE computeSurfaceNorm
!**************************************************************************
     SUBROUTINE tagging_th
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) :: g, n, m, i, j, k, n1, n2, n3, n4, nel2Cen, nel2Pnt, sumNodeId

        INTEGER            :: iPt, iPt1, flag_cell, i1, j1, k1
        REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis , minDis1, minDis, &
                              n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
                              cent_x, cent_y, cent_z, dis_cen, dis_pnt

        CHARACTER(LEN=120) :: filename1
        !g=2
        DO g=blk_start, nblocks
        block(g)%ibCellCount = 0
        block(g)%fluidCellCount = 0
        block(g)%solidCellCount = 0
        block(g)%cell = 0
	     block(g)%cell2 = 0
        block(g)%nodeIdTag = 0
        n2dotn = 0

!!$acc enter data copyin(cell, cell2, xp, yp, zp, x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, nodeIdTag)
!!$acc parallel loop collapse(3) present(nodeIdTag, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, cell, xp, yp, zp, x1, y1, z1)
 !$acc parallel loop collapse(3) default(present)
        DO 10 k = block(g)%k_startSearch, block(g)%k_endSearch
        DO 10 j = block(g)%j_startSearch, block(g)%j_endSearch
        DO 10 i = block(g)%i_startSearch, block(g)%i_endSearch
       !DO 10 k = 2, block(g)%nz
       !DO 10 j = 2, block(g)%ny
       !DO 10 i = 2, block(g)%nx
           minDis  = 1e14
           minDis1 = 1e14

           n1x = block(g)%xp(i)
           n1y = block(g)%yp(j)
           n1z = block(g)%zp(k)

           n2x = block(g)%x1(i)
           n2y = block(g)%y1(j)
           n2z = block(g)%z1(k)

           !$acc loop seq
           DO m = 1, block(g)%ibElems
           cent_x = block(g)%xcent(m)
           cent_y = block(g)%ycent(m)
           cent_z = block(g)%zcent(m)
              dis_cen  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2)
              dis_pnt  = dsqrt( (n2y-cent_y)**2 + (n2x-cent_x)**2  + (n2z-cent_z)**2)
              IF (dis_cen<minDis) THEN
                !print*,'1',i,j,k
                 minDis    = dis_cen
                 nel2Cen   = m
              ENDIF
              IF (dis_pnt<minDis1) THEN
                !print*,'2',i,j,k
                 minDis1   = dis_pnt
                 nel2Pnt   = m
              ENDIF
           ENDDO
            !    write(*,*)block(g)%x1(i),xcent(nel2Cen)
 ! 112          FORMAT(' ',I8,' ',I8,' ',I8,' ',I8)
            IF((block(g)%x1(i)<=block(g)%xcent(nel2Cen).AND.block(g)%x1(i+1)>=block(g)%xcent(nel2Cen)).AND. &
                (block(g)%y1(j)<=block(g)%ycent(nel2Cen).AND.block(g)%y1(j+1)>=block(g)%ycent(nel2Cen)).AND. &
                (block(g)%z1(k)<=block(g)%zcent(nel2Cen).AND.block(g)%z1(k+1)>=block(g)%zcent(nel2Cen))) THEN
                !print*,i,j,k
                block(g)%cell(i,j,k) = 2

           ENDIF

                      n2dotn  = (n2x - block(g)%xcent(nel2Pnt))*block(g)%cosAlpha(nel2Pnt) + &
                          (n2y - block(g)%ycent(nel2Pnt))*block(g)%cosBeta(nel2Pnt)  + &
                          (n2z - block(g)%zcent(nel2Pnt))*block(g)%cosGamma(nel2Pnt)

           IF (n2dotn>=-1e-16) THEN
              block(g)%nodeIdTag(i,j,k) = 0
           ELSE
              block(g)%nodeIdTag(i,j,k) = 1
           ENDIF
 10     CONTINUE
!$acc end parallel

!$acc parallel loop collapse(3) default(present)
           DO 20 k = block(g)%k_startSearch, block(g)%k_endSearch
           DO 20 j = block(g)%j_startSearch, block(g)%j_endSearch
           DO 20 i = block(g)%i_startSearch, block(g)%i_endSearch
           IF (block(g)%cell(i,j,k)/=2) THEN
              sumNodeId = 0
              sumNodeId = block(g)%nodeIdTag(i,j,k)      + block(g)%nodeIdTag(i+1,j,k)     &
                        + block(g)%nodeIdTag(i,j+1,k)    + block(g)%nodeIdTag(i+1,j+1,k)   &
                        + block(g)%nodeIdTag(i,j,k+1)    + block(g)%nodeIdTag(i+1,j,k+1)     &
                        + block(g)%nodeIdTag(i,j+1,k+1)  + block(g)%nodeIdTag(i+1,j+1,k+1)
              IF (sumNodeId==8) THEN
                 block(g)%cell(i,j,k) = 1
              ENDIF
           ENDIF
 20     CONTINUE
!$acc end parallel


         WRITE(filename1,1) g
  1      FORMAT('butter_f.',i3.3,".dat")
          OPEN(11,FILE=filename1,status='unknown')
        DO k = 1, block(g)%nz+2
        DO j = 1, block(g)%ny+2
        DO i = 1, block(g)%nx+2
        WRITE(11,*) block(g)%cell(i,j,k), block(g)%nodeIdTag(i,j,k)
        END DO
        END DO
        END DO
        CLOSE(11)



        block(g)%ibCellCount = 0
         block(g)%solidCellCount = 0
         block(g)%fluidCellCount = 0
 !!$acc parallel loop collapse(3) present(cell) reduction(+: solidCellCount, fluidCellCount, ibCellCount)
         DO 30 k = 2, block(g)%nz+1
         DO 30 j = 2, block(g)%ny+1
         DO 30 i = 2, block(g)%nx+1
            !n       = i-1  + nx*(j-2)  + nx*ny*(k-2)
            IF (block(g)%cell(i,j,k)==1) THEN
                block(g)%solidCellCount = block(g)%solidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==0) THEN
               block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==2) THEN
                block(g)%ibCellCount = block(g)%ibCellCount + 1
            ENDIF
 30      CONTINUE
 !!$acc end parallel
 !GOTO 1000
   !!$acc update self(cell)
      WRITE(filename1,22)char_f
 22   FORMAT(A3,'_inter_cell.dat')
      open(82,file=filename1,status='unknown')
      write(82,*)'variables = "x", "y","z", "var"'
    do k = 2,block(g)% nz+1
       do j = 2, block(g)%ny+1
       do i = 2, block(g)%nx+1
    n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
       if(block(g)%cell(i,j,k)==2)then
       write(82,*) block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), block(g)%cell(i,j,k)
       endif
       end do
       end do
    enddo
      close(82)
     ! GOTO 1000
      WRITE(filename1,23)char_f
 23   FORMAT(A3,'_fluid_cell.dat')
      open(83,file=filename1,status='unknown')
      write(83,*)'variables = "x", "y","z","var"'
       do k = 2, block(g)%nz+1
       do j = 2,block(g)% ny+1
       do i = 2, block(g)%nx+1
    n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
       if(block(g)%cell(i,j,k)==0)then
       write(83,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 0
       endif
       end do
       end do
    enddo
      close(83)

      WRITE(filename1,24)char_f
 24   FORMAT(A3,'_solid_cell.dat')
      open(84,file=filename1,status='unknown')
      write(84,*)'variables = "x", "y","z","var"'
       do k = 2, block(g)%nz+1
       do j = 2, block(g)%ny+1
       do i = 2, block(g)%nx+1
    n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
       if(block(g)%cell(i,j,k)==1)then
       write(84,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 1
       endif
       end do
       end do
    enddo
      close(84)

 1000 CONTINUE
         WRITE(filename1,2) g
 2       FORMAT('butter_cellcount_f.',i3.3,".dat")
         OPEN(12,FILE=filename1,FORM='formatted')
        WRITE(12,*) block(g)%solidCellCount, block(g)%fluidCellCount, block(g)%ibCellCount
        CLOSE(12)
         print*, 'search done'
         Print*, 'imms. cells=', block(g)%ibCellCount
        Print*, 'fluid cells=',block(g)%fluidCellCount
        Print*, 'solid cells=', block(g)%solidCellCount
        !stop
        END DO
     END SUBROUTINE tagging_th

!***********************************************************************
!**************************************************************************
     SUBROUTINE tagging_th_move
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) :: g, n, m, i, j, k, n1, n2, n3, n4, nel2Cen, nel2Pnt, sumNodeId

        INTEGER            :: iPt, iPt1, flag_cell, i1, j1, k1, a_blk_no, b_blk_no
        REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis , minDis1, minDis, &
                              n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
                              cent_x, cent_y, cent_z, dis_cen, dis_pnt

        CHARACTER(LEN=120) :: filename1
        !g=2
        !g=2
        DO g=blk_start,nblocks
        if ( block(g)%move_check == 1)then
            block(g)% ibCellCount = 0
        block(g)%fluidCellCount = 0
        block(g)% solidCellCount = 0
        block(g)%cell = 0
	block(g)%cell2 = 0
        block(g)%nodeIdTag = 0
        n2dotn = 0

!!$acc enter data copyin(cell, cell2, xp, yp, zp, x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, nodeIdTag)
!!$acc parallel loop collapse(3) present(nodeIdTag, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, cell, xp, yp, zp, x1, y1, z1)
!$acc parallel loop collapse(3) default(present)
        DO 10 k = block(g)%k_startSearch, block(g)%k_endSearch
        DO 10 j = block(g)%j_startSearch, block(g)%j_endSearch
        DO 10 i = block(g)%i_startSearch, block(g)%i_endSearch
       !DO 10 k = 2, block(g)%nz
       !DO 10 j = 2, block(g)%ny
       !DO 10 i = 2, block(g)%nx
           minDis  = 1e14
           minDis1 = 1e14

           n1x = block(g)%xp(i)
           n1y = block(g)%yp(j)
           n1z = block(g)%zp(k)

           n2x = block(g)%x1(i)
           n2y = block(g)%y1(j)
           n2z = block(g)%z1(k)

           !$acc loop seq
           DO m = 1, block(g)%ibElems
           cent_x = block(g)%xcent(m)
           cent_y = block(g)%ycent(m)
           cent_z = block(g)%zcent(m)
              dis_cen  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2)
              dis_pnt  = dsqrt( (n2y-cent_y)**2 + (n2x-cent_x)**2  + (n2z-cent_z)**2)
              IF (dis_cen<minDis) THEN
                !print*,'1',i,j,k
                 minDis    = dis_cen
                 nel2Cen   = m
              ENDIF
              IF (dis_pnt<minDis1) THEN
                !print*,'2',i,j,k
                 minDis1   = dis_pnt
                 nel2Pnt   = m
              ENDIF
           ENDDO
            !    write(*,*)block(g)%x1(i),xcent(nel2Cen)
 ! 112          FORMAT(' ',I8,' ',I8,' ',I8,' ',I8)
            IF((block(g)%x1(i)<=block(g)%xcent(nel2Cen).AND.block(g)%x1(i+1)>=block(g)%xcent(nel2Cen)).AND. &
                (block(g)%y1(j)<=block(g)%ycent(nel2Cen).AND.block(g)%y1(j+1)>=block(g)%ycent(nel2Cen)).AND. &
                (block(g)%z1(k)<=block(g)%zcent(nel2Cen).AND.block(g)%z1(k+1)>=block(g)%zcent(nel2Cen))) THEN
                !print*,i,j,k
                block(g)%cell(i,j,k) = 2

           ENDIF

                      n2dotn  = (n2x - block(g)%xcent(nel2Pnt))*block(g)%cosAlpha(nel2Pnt) + &
                          (n2y - block(g)%ycent(nel2Pnt))*block(g)%cosBeta(nel2Pnt)  + &
                          (n2z - block(g)%zcent(nel2Pnt))*block(g)%cosGamma(nel2Pnt)

           IF (n2dotn>=-1e-16) THEN
              block(g)%nodeIdTag(i,j,k) = 0
           ELSE
              block(g)%nodeIdTag(i,j,k) = 1
           ENDIF
 10     CONTINUE
!$acc end parallel

!$acc parallel loop collapse(3) default(present)
           DO 20 k = block(g)%k_startSearch, block(g)%k_endSearch
           DO 20 j = block(g)%j_startSearch, block(g)%j_endSearch
           DO 20 i = block(g)%i_startSearch, block(g)%i_endSearch
           IF (block(g)%cell(i,j,k)/=2) THEN
              sumNodeId = 0
              sumNodeId = block(g)%nodeIdTag(i,j,k)      + block(g)%nodeIdTag(i+1,j,k)     &
                        + block(g)%nodeIdTag(i,j+1,k)    + block(g)%nodeIdTag(i+1,j+1,k)   &
                        + block(g)%nodeIdTag(i,j,k+1)    + block(g)%nodeIdTag(i+1,j,k+1)     &
                        + block(g)%nodeIdTag(i,j+1,k+1)  + block(g)%nodeIdTag(i+1,j+1,k+1)
              IF (sumNodeId==8) THEN
                 block(g)%cell(i,j,k) = 1
              ENDIF
           ENDIF
 20     CONTINUE
!$acc end parallel

!    !$acc parallel loop collapse(3) default(present) private(i,j,k)
!       DO  k = k_startSearch, k_endSearch
!       DO  j = j_startSearch, j_endSearch
!       DO  i = i_startSearch, i_endSearch
!               flag_cell=0
!          IF ((block(g)%xp(i) .ge. 0.31 .and. block(g)%xp(i) .le. 0.362) .and. (block(g)%yp(j).ge. 0.485 .and. block(g)%yp(j) .le. 0.52) .and. (block(g)%zp(k) .ge. 0.487 .and. block(g)%zp(k) .le. 0.503) .and. block(g)%cell(i,j,k) .eq. 0)then
!
!         IF ((block(g)%xp(i) .ge. 0.348 .and. block(g)%xp(i) .le. 0.377))then
!             IF ((block(g)%yp(j) .ge. 0.497 .and. block(g)%yp(j) .le. 0.5023))then
!                IF ((block(g)%zp(k) .ge. 0.49704 .and. block(g)%zp(k) .le. 0.5023))then
!                       IF (block(g)%cell(i,j,k) .eq. 0)then

!                               block(g)%cell(i,j,k)=1
!        ENDIF
!        ENDIF
!        ENDIF
!        ENDIF

!
!               !$acc loop seq collapse(3)
!               DO k1=k-1,k+1
!               DO j1=j-1,j+1
!               DO i1=i-1,i+1
!                  if (block(g)% cell(i1,j1,k1) .ne. 0 )then
!                         flag_cell= flag_cell + 1
!                  else
!                         flag_cell= flag_cell + 0
!
!                  end if
!                 !if (block(g)% cell(i1,j1,k1) .eq. 0 .and. (i1 .ne. i) .and. (j1 .ne. j) .and. (k1 .ne. k) )then
!                 !       flag_cell=0
!                 !else
!                 !       flag_cell=1
!                 !end if
!               END DO
!               END DO
!               END DO
!               if (flag_cell .le. 26 .and. flag_cell .ge.23)then
!               !if (flag_cell .eq. 26)then
!                   block(g)%cell(i,j,k)=1
!               end if
!          END IF
!     END DO
!     END DO
!     END DO
!    !$acc end parallel


            !  ENDIF
           !else
           !     print*, i,j,k




        block(g)%ibCellCount = 0
         block(g)%solidCellCount = 0
         block(g)%fluidCellCount = 0
!!$acc parallel loop collapse(3) present(cell) reduction(+: solidCellCount, fluidCellCount, ibCellCount)
         DO 30 k = 2, block(g)%nz+1
         DO 30 j = 2, block(g)%ny+1
         DO 30 i = 2, block(g)%nx+1
            !n       = i-1  + nx*(j-2)  + nx*ny*(k-2)
            IF (block(g)%cell(i,j,k)==1) THEN
                block(g)%solidCellCount = block(g)%solidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==0) THEN
               block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==2) THEN
                block(g)%ibCellCount = block(g)%ibCellCount + 1
            ENDIF
 30      CONTINUE
!!$acc end parallel
      GOTO 1000
!!$acc update self(cell)
      WRITE(filename1,22)char_f
 22   FORMAT(A3,'_inter_cell.dat')
      open(82,file=filename1,status='unknown')
      write(82,*)'variables = "x", "y","z", "var"'
	do k = 2,block(g)% nz+1
       do j = 2, block(g)%ny+1
       do i = 2, block(g)%nx+1
	n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
       if(block(g)%cell(i,j,k)==2)then
       write(82,*) block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), block(g)%cell(i,j,k)
       endif
       end do
       end do
	enddo
      close(82)
     ! GOTO 1000
      WRITE(filename1,23)char_f
 23   FORMAT(A3,'_fluid_cell.dat')
      open(83,file=filename1,status='unknown')
      write(83,*)'variables = "x", "y","z","var"'
       do k = 2, block(g)%nz+1
       do j = 2,block(g)% ny+1
       do i = 2, block(g)%nx+1
	n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
       if(block(g)%cell(i,j,k)==0)then
       write(83,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 0
       endif
       end do
       end do
	enddo
      close(83)

      WRITE(filename1,24)char_f
 24   FORMAT(A3,'_solid_cell.dat')
      open(84,file=filename1,status='unknown')
      write(84,*)'variables = "x", "y","z","var"'
       do k = 2, block(g)%nz+1
       do j = 2, block(g)%ny+1
       do i = 2, block(g)%nx+1
	n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
       if(block(g)%cell(i,j,k)==1)then
       write(84,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 1
       endif
       end do
       end do
	enddo
      close(84)

 1000 CONTINUE
         print*, 'search done'
         Print*, 'imms. cells=', block(g)%ibCellCount
        Print*, 'fluid cells=',block(g)%fluidCellCount
        Print*, 'solid cells=', block(g)%solidCellCount
        !stop
        endif
        ENDDO
        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

        if ( block(b_blk_no)%move_check == 1)then
        !print*,'inside tag_mv'
       call fineUpdate_mv(g)
       call fineUpdate_bd_mv(g)
       !!print*,b_blk_no,'inside_tag_mv'
       !call solidCellBC_move(b_blk_no)
       !call updateVelocity_newv(b_blk_no)
       ! call writeOutput
       ! pause
        endif
        ENDDO
     END SUBROUTINE tagging_th_move

!***********************************************************************
!**************************************************************************

     SUBROUTINE findTScells
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER            :: g,i, j, k, i1, j1, k1, iPt1, m, n, i_x1, i_y1, i_z1, i_x2, i_y2, i_z2, il, jl, kl, tscnt
        REAL (KIND=8)      :: pos1_x, pos1_y, pos1_z, pos2_x, pos2_y, pos2_z, pt1
        CHARACTER(len=70) :: filename1


        print*,'inside findTScells'
        !g=2
        DO g=blk_start, nblocks
        !$acc parallel loop collapse(3) default(present)
                 DO 5 k = 2, block(g)%nz+1
                 DO 5 j = 2, block(g)%ny+1
                 DO 5 i = 2, block(g)%nx+1
        		    block(g)%cell2(i,j,k) = 0
         5      CONTINUE
        !$acc end parallel


        !!$acc parallel loop present(interceptedIndexPtr, ibSurfId, xp, yp, zp, cosAlpha, cosBeta, cosGamma, nelp, cell, cell2, deltax, deltay, deltaz)
        !$acc parallel loop default(present)
        DO 10 n = 1, block(g)%ibCellCount
        i = block(g)%interceptedIndexPtr(n, 1)
        j = block(g)%interceptedIndexPtr(n, 2)
        k = block(g)%interceptedIndexPtr(n, 3)
         GOTO 1
                  !pt1 = 1.51_rk*dsqrt(block(g)%deltax(i)**2 + block(g)%deltay(j)**2 + block(g)%deltaz(k)**2)
                  !pos1_x  = block(g)%xp(i) + pt1*cosAlpha(nelp(n))
                  !pos1_y  = block(g)%yp(j) + pt1*cosBeta(nelp(n))
                  !pos1_z  = block(g)%zp(k) + pt1*cosGamma(nelp(n))
                  !pos2_x  = block(g)%xp(i) - pt1*cosAlpha(nelp(n))
                  !pos2_y  = block(g)%yp(j) - pt1*cosBeta(nelp(n))
                  !pos2_z  = block(g)%zp(k) - pt1*cosGamma(nelp(n))
                  !
                  !!$acc loop seq
                  !DO il = i-7, i+7
                  !   if(pos1_x.ge.block(g)%xp(il).and.pos1_x.lt.block(g)%xp(il+1)) i_x1 = il
                  !   if(pos2_x.ge.block(g)%xp(il).and.pos2_x.lt.block(g)%xp(il+1)) i_x2 = il
                  !END DO
                  !!$acc loop seq
                  !DO jl = j-7, j+7
                  !   if(pos1_y.ge.block(g)%yp(jl).and.pos1_y.lt.block(g)%yp(jl+1)) i_y1 = jl
                  !   if(pos2_y.ge.block(g)%yp(jl).and.pos2_y.lt.block(g)%yp(jl+1)) i_y2 = jl
                  !END DO
                  !!$acc loop seq
                  !DO kl = k-7, k+7
                  !   if(pos1_z.ge.block(g)%zp(kl).and.pos1_z.lt.block(g)%zp(kl+1)) i_z1 = kl
                  !   if(pos2_z.ge.block(g)%zp(kl).and.pos2_z.lt.block(g)%zp(kl+1)) i_z2 = kl
                  !END DO
                  !block(g)%cell2(i, j, k) = 0
                  !IF(block(g)%cell(i_x1, i_y1, i_z1).EQ.0 .AND. block(g)%cell(i_x2, i_y2, i_z2).EQ.0) block(g)%cell2(i, j, k) = 2
         1  CONTINUE
         IF(block(g)%ibSurfID(block(g)%nelp(n))==51.OR.block(g)%ibSurfID(block(g)%nelp(n))==52) block(g)%cell2(i, j, k) = 2

         10     CONTINUE
        !$acc end parallel

        block(g)%TSCellCount = 0
        tscnt=0
        !!$acc parallel loop collapse(3) default(present) reduction(+: TSCellCount)
        !$acc parallel loop collapse(3) default(present) reduction(+:tscnt)
                 DO 20 k = 2, block(g)%nz+1
                 DO 20 j = 2, block(g)%ny+1
                 DO 20 i = 2, block(g)%nx+1
                    IF (block(g)%cell2(i,j,k)==2) THEN
                        !block(g)%TSCellCount = block(g)%TSCellCount + 1
                        tscnt = tscnt + 1
                    ENDIF
         20      CONTINUE
        !$acc end parallel

        block(g)%TSCellCount = tscnt
        !!$acc update self(cell2)
        print*, 'TScell count =', block(g)%TSCellCount

        ALLOCATE(block(g)%TSIndexPtr(block(g)%TSCellCount,3))

        iPt1 = 0
        !$acc loop collapse(3) seq
                 DO 30 k = 0, block(g)%nz+3
                 DO 30 j = 0, block(g)%ny+3
                 DO 30 i = 0, block(g)%nx+3
                    IF (block(g)%cell2(i,j,k)==2) THEN
                       iPt1 = iPt1 + 1
                       block(g)%TSIndexPtr(iPt1, 1) = i
                           block(g)%TSIndexPtr(iPt1, 2) = j
                               block(g)%TSIndexPtr(iPt1, 3) = k
                    ENDIF
         30      CONTINUE

        !!$acc enter data copyin(TSIndexPtr)

        ALLOCATE(block(g)%u2_ghost(block(g)%TSCellCount), block(g)%u2t_ghost(block(g)%TSCellCount), block(g)%v2_ghost(block(g)%TSCellCount), &
            block(g)%v2t_ghost(block(g)%TSCellCount), block(g)%p_ghost(block(g)%TSCellCount), block(g)%pt_ghost(block(g)%TSCellCount), &
                block(g)%u1_ghost(block(g)%TSCellCount), block(g)%u1t_ghost(block(g)%TSCellCount), block(g)%v1_ghost(block(g)%TSCellCount), &
                    block(g)%v1t_ghost(block(g)%TSCellCount), block(g)%w2_ghost(block(g)%TSCellCount), block(g)%w2t_ghost(block(g)%TSCellCount), &
                        block(g)%w1_ghost(block(g)%TSCellCount), block(g)%w1t_ghost(block(g)%TSCellCount), block(g)%index_ts(block(g)%TSCellCount))


        !!$acc enter data create(u2_ghost, u2t_ghost, v2_ghost, v2t_ghost, w2_ghost, w2t_ghost, p_ghost, pt_ghost, u1_ghost, u1t_ghost, v1_ghost, v1t_ghost, w1_ghost, w1t_ghost, index_ts)

        !!$acc parallel loop present(TSIndexPtr, interceptedIndexPtr, cell2, cell, nelp, nelu1, nelu2, nelv1, nelv2, u2_ghost, u2t_ghost, v2_ghost, v2t_ghost, p_ghost, pt_ghost, u1_ghost, u1t_ghost, v1_ghost, v1t_ghost, w2_ghost, w2t_ghost, w1_ghost, w1t_ghost, index_ts)
        !$acc parallel loop default(present)
        DO n = 1, block(g)%TSCellCount
        i = block(g)%TSIndexPtr(n, 1)
        j = block(g)%TSIndexPtr(n, 2)
        k = block(g)%TSIndexPtr(n, 3)
           !$acc loop seq
           DO m = 1, block(g)%ibCellCount
           i1 = block(g)%interceptedIndexPtr(m, 1)
           j1 = block(g)%interceptedIndexPtr(m, 2)
           k1 = block(g)%interceptedIndexPtr(m, 3)
              IF (i1==i .AND. j1==j .AND. k1==k) THEN
                  block(g)%index_ts(n) = m
              ENDIF
           ENDDO
           block(g)%u2_ghost(n)  = 0.
           block(g)%u2t_ghost(n) = 0.
           block(g)%v2_ghost(n)  = 0.
           block(g)%v2t_ghost(n) = 0.
           block(g)%p_ghost(n)   = 0.
           block(g)%pt_ghost(n)  = 0.
           block(g)%u1_ghost(n)  = 0.
           block(g)%u1t_ghost(n) = 0.
           block(g)%v1_ghost(n)  = 0.
           block(g)%v1t_ghost(n) = 0.
           block(g)%w2_ghost(n)  = 0.
           block(g)%w2t_ghost(n) = 0.
           block(g)%w1_ghost(n)  = 0.
           block(g)%w1t_ghost(n) = 0.
        ENDDO
        !$acc end parallel


!               print*, 'TScell done'
         GOTO 1000
        !!$acc update self(cell,cell2)
             WRITE(filename1,111) ita
         111   FORMAT('TS/ts.',i9.9,".dat")
             OPEN(UNIT = 82, FILE = filename1, STATUS = 'unknown')
             open(82,file=filename1,status='unknown')
             write(82,*)'variables = "x", "y","z", "var"'
               do k = 2, block(g)%nz+1
               do j = 2, block(g)%ny+1
               do i = 2, block(g)%nx+1
        	!n = i-1  + nx*(j-2)  + nx*ny*(k-2)
                  if(block(g)%cell2(i,j,k)==2)then
                      write(82,*) block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), block(g)%cell2(i,j,k)
                  endif
               end do
               end do
        	    enddo
              close(82)
         1000 CONTINUE
        	    enddo
         ! STOP
             END SUBROUTINE findTScells
!***********************************************************************
!***********************************************************************

     SUBROUTINE selectiveRetagging_th
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, g, m, i, j, k, n1, n2, n3, n4, i1, j1, k1,i2,j2,k2,flag_cell, nn, &
                               nel2n, sumId, nel2Pnt, nel2Cen, sumNodeID
        INTEGER            :: iPt, iPt1
        INTEGER            :: flcnt, sdcnt, ibcnt
        REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis , minDis, minDis1, dis_cen, dis_pnt, &
                              n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
                              cent_x, cent_y, cent_z

        CHARACTER(len=120) :: filename1
        CHARACTER(len=120) :: filename2

        !g=2
        !g=2
       DO g=blk_start,nblocks
        if( block(g)%blk_mv_tag ==0)then
!!!$acc parallel loop present(interceptedIndexPtr, nodeIdTag, xp, yp, zp, x1, y1, z1, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, cell)
!$acc parallel loop gang vector default(present)
        DO nn = 1, block(g)%ibCellCount
        i1 = block(g)%interceptedIndexPtr(nn, 1)
        j1 = block(g)%interceptedIndexPtr(nn, 2)
        k1 = block(g)%interceptedIndexPtr(nn, 3)
           block(g)%cell(i1,j1,k1) = 0
           !$acc loop collapse(3) seq
           DO 10 k = k1-1, k1+1
           DO 10 j = j1-1, j1+1
           DO 10 i = i1-1, i1+1

           minDis  = 1e14
           minDis1 = 1e14

           n1x = block(g)%xp(i)
           n1y = block(g)%yp(j)
           n1z = block(g)%zp(k)

           n2x = block(g)%x1(i)
           n2y = block(g)%y1(j)
           n2z = block(g)%z1(k)

           !$acc loop seq
           DO m = 1, block(g)%ibElems
           cent_x = block(g)%xcent(m)
           cent_y = block(g)%ycent(m)
           cent_z = block(g)%zcent(m)
              dis_cen  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2)
              dis_pnt  = dsqrt( (n2y-cent_y)**2 + (n2x-cent_x)**2  + (n2z-cent_z)**2)
              IF (dis_cen<minDis) THEN
                 minDis    = dis_cen
                 nel2Cen   = m
              ENDIF
              IF (dis_pnt<minDis1) THEN
                 minDis1   = dis_pnt
                 nel2Pnt   = m
              ENDIF
           ENDDO
           IF((block(g)%x1(i)<=block(g)%xcent(nel2Cen).AND.block(g)%x1(i+1)>=block(g)%xcent(nel2Cen)).AND. &
               (block(g)%y1(j)<=block(g)%ycent(nel2Cen).AND.block(g)%y1(j+1)>=block(g)%ycent(nel2Cen)).AND. &
               (block(g)%z1(k)<=block(g)%zcent(nel2Cen).AND.block(g)%z1(k+1)>=block(g)%zcent(nel2Cen))) THEN
                   block(g)%cell(i,j,k) = 2
           ENDIF

                      n2dotn  = (n2x - block(g)%xcent(nel2Pnt))*block(g)%cosAlpha(nel2Pnt) + &
                          (n2y - block(g)%ycent(nel2Pnt))*block(g)%cosBeta(nel2Pnt)  + &
                          (n2z - block(g)%zcent(nel2Pnt))*block(g)%cosGamma(nel2Pnt)

           IF (n2dotn>=-1e-16) THEN
              block(g)%nodeIdTag(i,j,k) = 0
           ELSE
              block(g)%nodeIdTag(i,j,k) = 1
           ENDIF
 10        CONTINUE
        ENDDO
!$acc end parallel

!!$acc parallel loop present(interceptedIndexPtr, nodeIdTag, cell)
!$acc parallel loop gang vector default(present)
        DO nn = 1, block(g)%ibCellCount
        i1 = block(g)%interceptedIndexPtr(nn, 1)
        j1 = block(g)%interceptedIndexPtr(nn, 2)
        k1 = block(g)%interceptedIndexPtr(nn, 3)
           !$acc loop collapse(3) seq
           DO 15 k = k1-1, k1+1
           DO 15 j = j1-1, j1+1
           DO 15 i = i1-1, i1+1
             IF (block(g)%cell(i,j,k)/=2) THEN
                sumNodeId = 0
                sumNodeId = block(g)%nodeIdTag(i,j,k)      + block(g)%nodeIdTag(i+1,j,k)     &
                          + block(g)%nodeIdTag(i,j+1,k)    + block(g)%nodeIdTag(i+1,j+1,k)   &
                          + block(g)%nodeIdTag(i,j,k+1)    + block(g)%nodeIdTag(i+1,j,k+1)     &
                          + block(g)%nodeIdTag(i,j+1,k+1)  + block(g)%nodeIdTag(i+1,j+1,k+1)
                IF (sumNodeId==8) THEN
                   block(g)%cell(i,j,k) = 1
                ELSE
                   block(g)%cell(i,j,k) = 0
                ENDIF
             ENDIF
 15        CONTINUE
        ENDDO
!$acc end parallel
!    !$acc parallel loop collapse(3) default(present) private(i,j,k)
!       DO  k = k_startSearch, k_endSearch
!       DO  j = j_startSearch, j_endSearch
!       DO  i = i_startSearch, i_endSearch
!               flag_cell=0
!          IF ((block(g)%xp(i) .ge. 0.31 .and. block(g)%xp(i) .le. 0.362) .and. (block(g)%yp(j).ge. 0.485 .and. block(g)%yp(j) .le. 0.52) .and. (block(g)%zp(k) .ge. 0.487 .and. block(g)%zp(k) .le. 0.503) .and. block(g)%cell(i,j,k) .eq. 0)then
!
!         IF ((block(g)%xp(i) .ge. 0.348 .and. block(g)%xp(i) .le. 0.377))then
!             IF ((block(g)%yp(j) .ge. 0.497 .and. block(g)%yp(j) .le. 0.5023))then
!                IF ((block(g)%zp(k) .ge. 0.49704 .and. block(g)%zp(k) .le. 0.5023))then
!                       IF (block(g)%cell(i,j,k) .eq. 0)then

!                               block(g)%cell(i,j,k)=1
!        ENDIF
!        ENDIF
!        ENDIF
!        ENDIF

!
!               !$acc loop seq collapse(3)
!               DO k1=k-1,k+1
!               DO j1=j-1,j+1
!               DO i1=i-1,i+1
!                  if (block(g)% cell(i1,j1,k1) .ne. 0 )then
!                         flag_cell= flag_cell + 1
!                  else
!                         flag_cell= flag_cell + 0
!
!                  end if
!                 !if (block(g)% cell(i1,j1,k1) .eq. 0 .and. (i1 .ne. i) .and. (j1 .ne. j) .and. (k1 .ne. k) )then
!                 !       flag_cell=0
!                 !else
!                 !       flag_cell=1
!                 !end if
!               END DO
!               END DO
!               END DO
!               if (flag_cell .le. 26 .and. flag_cell .ge.23)then
!               !if (flag_cell .eq. 26)then
!                   block(g)%cell(i,j,k)=1
!               end if
!          END IF
!     END DO
!     END DO
!     END DO
!    !$acc end parallel

!    !$acc parallel loop default(present) private(flag_cell)
!       DO nn = 1, ibCellCount
!          i1 = interceptedIndexPtr(nn, 1)
!          j1 = interceptedIndexPtr(nn, 2)
!          k1 = interceptedIndexPtr(nn, 3)
!    !!$acc parallel loop collapse(3) default(present)
!          DO k = k1-1, k1+1
!          DO j = j1-1, j1+1
!          DO i = i1-1, i1+1
!               flag_cell=0
!          IF ((block(g)%xp(i) .ge. 0.31 .and. block(g)%xp(i) .le. 0.362) .and. (block(g)%yp(j).ge. 0.485 .and. block(g)%yp(j) .le. 0.52) .and. (block(g)%zp(k) .ge. 0.487 .and. block(g)%zp(k) .le. 0.503) .and. block(g)%cell(i,j,k) .eq. 0)then
!
!               !$acc loop seq collapse(3)
!               DO k2=k-1,k+1
!               DO j2=j-1,j+1
!               DO i2=i-1,i+1
!                  !if (block(g)% cell(i2,j2,k2) .eq. 0 .and. (i2 .ne. i) .and. (j2 .ne. j) .and. (k2 .ne. k) )then
!                  if (block(g)% cell(i2,j2,k2) .ne. 0 )then
!                         flag_cell= flag_cell + 1
!                  else
!                         flag_cell= flag_cell + 0
!                  end if
!               END DO
!               END DO
!               END DO
!               !print*, flag_cell
!               if (flag_cell .le. 26 .and. flag_cell .ge.23)then
!                   block(g)%cell(i,j,k)=1
!               end if
!          END IF
!               END DO
!               END DO
!               END DO
!    !!$acc end parallel
!     END DO
!    !$acc end parallel
  block(g)%ibCellCount = 0
  block(g)%solidCellCount = 0
block(g)%fluidCellCount = 0
sdcnt=0
flcnt=0
ibcnt=0
!!$acc parallel loop collapse(3) present(cell) reduction(+: solidCellCount, fluidCellCount, ibCellCount)
!$acc parallel loop gang vector collapse(3) default(present) private(i,j,k,n) reduction(+: sdcnt, flcnt, ibcnt)
         DO 20 k = 2, block(g)%nz+1
         DO 20 j = 2, block(g)%ny+1
         DO 20 i = 2, block(g)%nx+1
            n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
            IF (block(g)%cell(i,j,k)==1) THEN
                !block(g)%solidCellCount = block(g)%solidCellCount + 1
                sdcnt = sdcnt + 1
            ELSEIF (block(g)%cell(i,j,k)==0) THEN
               !block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
               flcnt  = flcnt + 1
            ELSEIF (block(g)%cell(i,j,k)==2) THEN
               ! block(g)%ibCellCount = block(g)%ibCellCount + 1
                ibcnt = ibcnt + 1
            ENDIF
 20      CONTINUE
!$acc end parallel
  block(g)%ibCellCount = ibcnt
  block(g)%solidCellCount = sdcnt
block(g)%fluidCellCount = flcnt
!!$acc update self(cell,cell2)
 GOTO 111
!      IF (mod(ita,200).eq.0) THEN
!       WRITE(filename1,991)ita
!991    FORMAT('inter',i9.9,'.dat')
!        !open(82,file='inter.dat',status='unknown')
!        open(82,file=filename1,status='unknown')
!       write(82,*)'variables = "x", "y","z", "var"'
!         do k = 2, block(g)%nz+1
!        do j = 2,block(g)% ny+1
!        do i = 2, block(g)%nx+1
!         n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
!        if(block(g)%cell(i,j,k).eq.2)then
!        write(82,*) block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), block(g)%cell2(i,j,k)
!        endif
!        end do
!        end do
!         enddo
!       close(82)
!goto 111
!       open(83,file='fluid.dat',status='unknown')
!       write(83,*)'variables = "x", "y","z","var"'
!        do k = 2, block(g)%nz+1
!        do j = 2, block(g)%ny+1
!        do i = 2, block(g)%nx+1
!         n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
!        if(block(g)%cell(i,j,k).eq.0)then
!        write(83,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 0
!        endif
!        end do
!        end do
!         enddo
!       close(83)
 !111 continue
!       WRITE(filename2,992)ita
!992    FORMAT('solid',i9.9,'.dat')
!       !open(84,file='solid.dat',status='unknown')
!        open(82,file=filename2,status='unknown')
!       write(84,*)'variables = "x", "y","z","var"'
!        do k = 2, block(g)%nz+1
!        do j = 2, block(g)%ny+1
!        do i = 2, block(g)%nx+1
!         n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
!        if(block(g)%cell(i,j,k).eq.1)then
!        write(84,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 1
!        endif
!        end do
!        end do
!         enddo
!       close(84)
!      ENDIF
 111 CONTINUE
       print*, 'selective retagging', block(g)%ibCellCount, block(g)%fluidCellCount,block(g)%solidCellCount
     END IF
       ENDDO
     END SUBROUTINE selectiveRetagging_th
!***********************************************************************
     SUBROUTINE cellCount_solid
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j, k, g

        print*, "cellCount started"


         DO g=blk_start,nblocks
         iPt  = 0
         iPt1 = 0
         iPt2 = 0

         ALLOCATE(block(g)%interceptedIndexPtr(block(g)%ibCellCount,3),block(g)%fluidIndexPtr(block(g)%fluidCellCount, 3), &
         block(g)%solidIndexPtr(block(g)%solidCellCount, 3))

        !!$acc update host(cell)

         DO 30 k = 2, block(g)%nz+1
         DO 30 j = 2, block(g)%ny+1
         DO 30 i = 2, block(g)%nx+1
            !n  = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)

            IF (block(g)%cell(i,j,k)==0) THEN
               iPt1 = iPt1 + 1
               block(g)%fluidIndexPtr(iPt1, 1) = i
               block(g)%fluidIndexPtr(iPt1, 2) = j
               block(g)%fluidIndexPtr(iPt1, 3) = k
            ELSEIF (block(g)%cell(i,j,k)==1) THEN
               iPt2 = iPt2 + 1
               block(g)%solidIndexPtr(iPt2, 1) = i
               block(g)%solidIndexPtr(iPt2, 2) = j
               block(g)%solidIndexPtr(iPt2, 3) = k
            ELSEIF (block(g)%cell(i,j,k)==2) THEN
               iPt = iPt + 1
               block(g)%interceptedIndexPtr(iPt, 1) = i
               block(g)%interceptedIndexPtr(iPt, 2) = j
               block(g)%interceptedIndexPtr(iPt, 3) = k
            ENDIF
 30      CONTINUE
         block(g)%redCellCount = 0
         block(g)%blackCellCount  = 0

         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
	     k = block(g)%fluidIndexPtr(n, 3)

            IF (mod(i+j+k,2)==1) THEN
               block(g)%redCellCount = block(g)%redCellCount + 1
            ELSE
               block(g)%blackCellCount = block(g)%blackCellCount + 1
            ENDIF
         ENDDO
         ALLOCATE (block(g)%redCellIndexPtr(block(g)%redCellCount,3) ,block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))
         ipt1 = 0
         iPt = 0
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)

            IF (mod(i+j+k,2)==1) THEN
               iPt = iPt + 1
               block(g)%redCellIndexPtr(iPt, 1) = i
               block(g)%redCellIndexPtr(iPt, 2) = j
               block(g)%redCellIndexPtr(iPt, 3) = k

            ELSE
               iPt1 = iPt1 + 1
               block(g)%blackCellIndexPtr(iPt1, 1) = i
               block(g)%blackCellIndexPtr(iPt1, 2) = j
               block(g)%blackCellIndexPtr(iPt1, 3) = k

            ENDIF
         ENDDO

        !!$acc update device(interceptedIndexPtr, block(g)%fluidIndexPtr, solidIndexPtr)
        !!$acc update device(redCellIndexPtr, blackCellIndexPtr)

            print*,g, block(g)%fluidCellCount, block(g)%redCellCount, block(g)%blackCellCount
        END DO
            !print*, "cellCount done"
     END SUBROUTINE cellCount_solid
!**************************************************************************

      SUBROUTINE cellCount
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j
        INTEGER (kind = 8)::  g,ng, k, rccount, bccount


         g=1
         iPt  = 0
         iPt1 = 0
         iPt2 = 0
         !!$acc update self(cell)

    !!     ALLOCATE(interceptedIndexPtr(ibCellCount,2),fluidIndexPtr(block(g)%fluidCellCount, 2), solidIndexPtr(solidCellCount,2))

       !!$acc loop default(present) collapse(3) seq
       ! DO k = 2, block(g)%nz +1
       ! DO j = 2, block(g)%ny +1
       ! DO i = 2, block(g)%nx +1
       !        block(g)%fluidCellCount=block(g)%fluidCellCount +1
       ! end do
       ! end do
       ! end do

         block(g)%fluidCellCount=block(g)%nx * block(g)%ny * block(g)%nz
         ALLOCATE(block(g)%fluidIndexPtr(block(g)%fluidCellCount,3))

         iPt1=0
         DO  k = 2, block(g)%nz +1
         DO  j = 2, block(g)%ny +1
         DO  i = 2, block(g)%nx +1
          !   IF (cell(i,j).eq.0) THEN
              iPt1 = iPt1 + 1
              block(g)%fluidIndexPtr(iPt1, 1) = i
              block(g)%fluidIndexPtr(iPt1, 2) = j
              block(g)%fluidIndexPtr(iPt1, 3) = k
         END DO
         END DO
         END DO

          block(g)%redCellCount = 0
          block(g)%blackCellCount  = 0
          rccount = 0
          bccount  = 0

        !!$acc enter data copyin(solidIndexPtr,interceptedIndexPtr,fluidIndexPtr)
         !$acc parallel loop default(present) reduction(+:rccount,bccount)
         !DO g = fl_blk(1), fl_blk(nfl_blk)
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
            IF (mod(i+j,2)==1) THEN
               !block(g)%redCellCount = block(g)%redCellCount + 1
               rccount = rccount + 1
            ELSE
               !block(g)%blackCellCount = block(g)%blackCellCount +1
               bccount = bccount + 1

            ENDIF
         ENDDO
        !$acc end parallel
          block(g)%redCellCount = rccount
          block(g)%blackCellCount  = bccount
         ALLOCATE(block(g)%redCellIndexPtr(block(g)%redCellCount,3),block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))

         !DO g = fl_blk(1), fl_blk(nfl_blk)
         !$acc loop seq
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
            IF (mod(i+j+k,2)==1) THEN
               iPt = iPt + 1
               block(g)%redCellIndexPtr(iPt, 1) = i
               block(g)%redCellIndexPtr(iPt, 2) = j
               block(g)%redCellIndexPtr(iPt, 3) = k
            ELSE
               iPt1 = iPt1 + 1
               block(g)%blackCellIndexPtr(iPt1, 1) = i
               block(g)%blackCellIndexPtr(iPt1, 2) = j
               block(g)%blackCellIndexPtr(iPt1, 3) = k
            ENDIF
         ENDDO
        !!$acc enter data copyin(redCellIndexPtr, blackCellIndexPtr)
         !DO g = fl_blk(1), fl_blk(nfl_blk)
            print*, g,block(g)%fluidCellCount,block(g)%redCellCount, block(g)%blackCellCount

       END SUBROUTINE cellCount


!***********************************************************************************************************************
     SUBROUTINE computeNormDistance
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER            ::  nel2u1, nel2u2, nel2v1, nel2v2, nel2w1, nel2w2
        INTEGER            :: i, j, k, n, nel2p, g, ibxx, m
        REAL (KIND=8)      :: n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, dis, dis1, dis2, dis3, dis4, dis5, dis6, minDis, minDis1, minDis2, minDis3, minDis4, minDis5, minDis6, cent_x, cent_y,cent_z
        CHARACTER(len=70)  :: filename1

        print*, 'computeNormDistance started'
        DO g=blk_start,nblocks

        ibxx=block(g)%ibCellCount
        print*,ibxx
    	ALLOCATE(block(g)%pNormDis(ibxx), block(g)%nelp(ibxx), &
        block(g)%nelu1(ibxx), block(g)%nelu2(ibxx), block(g)%nelv1(ibxx), &
        block(g)%nelv2(ibxx), block(g)%nelw1(ibxx), block(g)%nelw2(ibxx), &
        block(g)%u1NormDis(ibxx), block(g)%u2NormDis(ibxx) , block(g)%v1NormDis(ibxx), &
        block(g)%v2NormDis(ibxx) , block(g)%w1NormDis(ibxx), block(g)%w2NormDis(ibxx))

       !!$acc parallel loop gang vector  &
       !!$acc private (i, j, k, n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, nel2p)          &
       !!$acc present (xp, yp, zp, x1, y1, z1, interceptedIndexPtr, minElemcell, nelp, nelu1,          &
       !!$acc          nelu2, nelv1, nelv2, nelw1, nelw2, pNormDis, u1NormDis, u2NormDis, v1NormDis,   &
       !!$acc          v2NormDis, w1NormDis, w2NormDis, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, cell)
       ! print*,'1'
        !$acc parallel loop gang vector default(present) &
        !$acc private (i, j, k, n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, nel2p)
        DO 10 k = 1, block(g)%ibCellCount

           m = 0
           minDis  = 1e14
           minDis1 = 1e14
           minDis2 = 1e14
           minDis3 = 1e14
           minDis4 = 1e14
           minDis5 = 1e14
           minDis6 = 1e14
           n1x = block(g)%xp(block(g)%interceptedIndexPtr(k, 1))
           n2x = block(g)%x1(block(g)%interceptedIndexPtr(k, 1))
           n3x = block(g)%x1(block(g)%interceptedIndexPtr(k, 1)+1)
           n1y = block(g)%yp(block(g)%interceptedIndexPtr(k, 2))
           n2y = block(g)%y1(block(g)%interceptedIndexPtr(k, 2))
           n3y = block(g)%y1(block(g)%interceptedIndexPtr(k, 2)+1)
           n1z = block(g)%zp(block(g)%interceptedIndexPtr(k, 3))
           n2z = block(g)%z1(block(g)%interceptedIndexPtr(k, 3))
           n3z = block(g)%z1(block(g)%interceptedIndexPtr(k, 3)+1)
          !print*,  normDisPtr(k, 2), n3x, normDisPtr(k, 2), n3y
          !$acc loop seq
        DO m = 1, block(g)%ibElems
        cent_x = block(g)%xcent(m)
        cent_y = block(g)%ycent(m)
        cent_z = block(g)%zcent(m)
              dis   = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2 + (n1z-cent_z)**2 )
              dis1  = dsqrt( (n1y-cent_y)**2 + (n2x-cent_x)**2 + (n1z-cent_z)**2 )
              dis2  = dsqrt( (n1y-cent_y)**2 + (n3x-cent_x)**2 + (n1z-cent_z)**2 )
              dis3  = dsqrt( (n2y-cent_y)**2 + (n1x-cent_x)**2 + (n1z-cent_z)**2 )
              dis4  = dsqrt( (n3y-cent_y)**2 + (n1x-cent_x)**2 + (n1z-cent_z)**2 )
              dis5  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2 + (n2z-cent_z)**2 )
              dis6  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2 + (n3z-cent_z)**2 )
	       IF (dis<minDis) THEN
                 minDis   = dis
                 nel2p    = m
              ENDIF
              IF (dis1<minDis1) THEN
                 minDis1   = dis1
                 nel2u1    = m
              ENDIF
              IF (dis2<minDis2) THEN
                 minDis2   = dis2
                 nel2u2    = m
              ENDIF
              IF (dis3<minDis3) THEN
                 minDis3   = dis3
                 nel2v1    = m
              ENDIF
              IF (dis4<minDis4) THEN
                 minDis4   = dis4
                 nel2v2    = m
              ENDIF
              IF (dis5<minDis5) THEN
                 minDis5   = dis5
                 nel2w1    = m
              ENDIF
              IF (dis6<minDis6) THEN
                 minDis6   = dis6
                 nel2w2    = m
              ENDIF
           ENDDO

           block(g)% nelp(k)  = nel2p
           block(g)% nelu1(k) = nel2u1
           block(g)% nelu2(k) = nel2u2
           block(g)% nelv1(k) = nel2v1
           block(g)% nelv2(k) = nel2v2
           block(g)% nelw1(k) = nel2w1
           block(g)% nelw2(k) = nel2w2
           !nu1nu2 = cosAlpha(nel2u1)*cosAlpha(nel2u2) + cosBeta(nel2u1)*cosBeta(nel2u2) + cosGamma(nel2u1)*cosGamma(nel2u2)
           !nu1nv1 = cosAlpha(nel2u1)*cosAlpha(nel2v1) + cosBeta(nel2u1)*cosBeta(nel2v1) + cosGamma(nel2u1)*cosGamma(nel2v1)
           !nu1nv2 = cosAlpha(nel2u1)*cosAlpha(nel2v2) + cosBeta(nel2u1)*cosBeta(nel2v2) + cosGamma(nel2u1)*cosGamma(nel2v2)
           !nu1nw1 = cosAlpha(nel2u1)*cosAlpha(nel2w1) + cosBeta(nel2u1)*cosBeta(nel2w1) + cosGamma(nel2u1)*cosGamma(nel2w1)
           !nu1nw2 = cosAlpha(nel2w2)*cosAlpha(nel2w2) + cosBeta(nel2u1)*cosBeta(nel2w2) + cosGamma(nel2u1)*cosGamma(nel2w2)
           !block(g)%cell2(i, j, k) = 0
           !IF (nu1nu2.LT.0 .OR. nu1nv1.LT.0 .OR. nu1nv2.LT.0 .OR. nu1nw1.LT.0 .OR. nu1nw2.LT.0 ) block(g)%cell2(i,j,k) = 2
          block(g)% pNormDis(k)  = (n1x -block(g)% xcent(nel2p))*block(g)%cosAlpha(nel2p) + &
                   (n1y -block(g)% ycent(nel2p))*block(g)%cosBeta(nel2p)  + &
                   (n1z -block(g)% zcent(nel2p))*block(g)%cosGamma(nel2p)

          block(g)% u1NormDis(k) = (n2x -block(g)% xcent(nel2u1))*block(g)%cosAlpha(nel2u1) + &
                       (n1y -block(g)% ycent(nel2u1))*block(g)%cosBeta(nel2u1)  + &
                       (n1z - block(g)%zcent(nel2u1))*block(g)%cosGamma(nel2u1)

          block(g)% u2NormDis(k) = (n3x - block(g)%xcent(nel2u2))*block(g)%cosAlpha(nel2u2) + &
              (n1y -block(g)% ycent(nel2u2))*block(g)%cosBeta(nel2u2)  + &
              (n1z -block(g)% zcent(nel2u2))*block(g)%cosGamma(nel2u2)

          block(g)% v1NormDis(k) = (n1x - block(g)%xcent(nel2v1))*block(g)%cosAlpha(nel2v1) + &
              (n2y - block(g)%ycent(nel2v1))*block(g)%cosBeta(nel2v1)  + &
              (n1z - block(g)%zcent(nel2v1))*block(g)%cosGamma(nel2v1)

          block(g)% v2NormDis(k) = (n1x - block(g)%xcent(nel2v2))*block(g)%cosAlpha(nel2v2) + &
                   (n3y -block(g)% ycent(nel2v2))*block(g)%cosBeta(nel2v2)  + &
                   (n1z -block(g)% zcent(nel2v2))*block(g)%cosGamma(nel2v2)

          block(g)% w1NormDis(k) = (n1x -block(g)% xcent(nel2w1))*block(g)%cosAlpha(nel2w1) + &
                       (n1y -block(g)% ycent(nel2w1))*block(g)%cosBeta(nel2w1)  + &
                       (n2z -block(g)% zcent(nel2w1))*block(g)%cosGamma(nel2w1)

          block(g)% w2NormDis(k) = (n1x -block(g)% xcent(nel2w2))*block(g)%cosAlpha(nel2w2) + &
                           (n1y -block(g)% ycent(nel2w2))*block(g)%cosBeta(nel2w2)  + &
                           (n3z -block(g)% zcent(nel2w2))*block(g)%cosGamma(nel2w2)
 10      CONTINUE
        !$acc end parallel

        ! DEALLOCATE (block(g)%minElemcell)
         END DO
         print*, 'computeNormDistance done'

     END SUBROUTINE computeNormDistance

!******************************************************************

        SUBROUTINE fine_block_cell


        USE global
        IMPLICIT NONE
        INTEGER (KIND=8) :: i, j, k, g, f, factor, a_blk_no, b_blk_no
        integer (kind=4) :: nx_var,ny_var, nx_var_r,ny_var_r,nx_var_t,ny_var_t,nx_var_tn,ny_var_tn
        integer (kind=4) :: nz_var,nz_var_r,nz_var_t,nz_var_tn
        integer (kind=4) :: st_rc_x, en_rc_x, st_rc_y, en_rc_y
        integer (kind=4) :: st_rc_z, en_rc_z


        block(1)%cell_n=0
        block(1)%cell=0
        block(1)%cell_pr=0

        DO g=1, intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk
           factor=intfr(g)%b_msh/intfr(g)%a_msh


       !nx_var_r=block(a_blk_no)%irc
       !ny_var_r=block(a_blk_no)%jrc
       !nx_var_t=block(a_blk_no)%itc
       !ny_var_t=block(a_blk_no)%jtc
       !nx_var_tn=block(a_blk_no)%itn
       !ny_var_tn=block(a_blk_no)%jtn
       !st_rc_x=block(a_blk_no)%irc_st
       !en_rc_x=block(a_blk_no)%irc_en
       !st_rc_y=block(a_blk_no)%jrc_st
       !en_rc_y=block(a_blk_no)%jrc_en


       ! DO j = st_rc_y, en_rc_y
       ! DO i = st_rc_x, en_rc_x
       !!$acc parallel loop gang vector default(present) firstprivate(a_blk_no,b_blk_no,factor)
         DO k = 2, block(a_blk_no)%nz +1
         DO j = 2, block(a_blk_no)%ny +1
         DO i = 2, block(a_blk_no)%nx +1

        if ( block(a_blk_no)%xp(i) >= intfr(g)%xintf_start+ block(b_blk_no)%cintp*block(a_blk_no)%dx .and. &
             block(a_blk_no)%xp(i) <= intfr(g)%xintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dx   .and. &
             block(a_blk_no)%zp(k) >= intfr(g)%zintf_start+ block(b_blk_no)%cintp*block(a_blk_no)%dz .and. &
             block(a_blk_no)%zp(k) <= intfr(g)%zintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dz   .and. &
             block(a_blk_no)%yp(j) >= intfr(g)%yintf_start + block(b_blk_no)%cintp*block(a_blk_no)%dy .and. &
             block(a_blk_no)%yp(j) <= intfr(g)%yintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dy     )then

                block(a_blk_no)%cell_n(i,j,k)=1
                block(a_blk_no)%cell(i,j,k)=1
                !print*,block(a_blk_no)%xp(i),block(a_blk_no)%yp(j),block(a_blk_no)%cell(i,j)

        endif
        if ( block(a_blk_no)%xp(i) >= intfr(g)%xintf_start .and. &
             block(a_blk_no)%xp(i) <= intfr(g)%xintf_end   .and. &
             block(a_blk_no)%zp(k) >= intfr(g)%zintf_start .and. &
             block(a_blk_no)%zp(k) <= intfr(g)%zintf_end   .and. &
             block(a_blk_no)%yp(j) >= intfr(g)%yintf_start .and. &
             block(a_blk_no)%yp(j) <= intfr(g)%yintf_end   )then

                block(a_blk_no)%cell_pr(i,j,k)=1
                !print*,block(a_blk_no)%xp(i),block(a_blk_no)%yp(j),block(a_blk_no)%cell(i,j)

        endif
        ENDDO
        ENDDO
        ENDDO
        !!$acc end parallel
        ENDDO
       ! DO j = st_rc_y, en_rc_y
       ! DO i = st_rc_x, en_rc_x
      !! DO j = 2, block(1)%ny +1
      !! DO i = 2, block(1)%nx +1
      !!        block(1)%cell(i,j)=block(1)%cell_n(i,j)
      !!ENDDO
      !!ENDDO


        END SUBROUTINE
!**************************************************************************
!**************************************************************************
        SUBROUTINE cellCount_coarse_new
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j, k
        INTEGER (kind = 8)::  g,ng
        integer (kind=4) :: nx_var,ny_var, nx_var_r,ny_var_r,nx_var_t,ny_var_t,nx_var_tn,ny_var_tn
        integer (kind=4) :: nz_var,nz_var_r,nz_var_t,nz_var_tn
        integer (kind=4) :: st_rc_x, en_rc_x, st_rc_y, en_rc_y
        integer (kind=4) :: st_rc_z, en_rc_z


        g=1
       !nx_var_r=block(g)%irc
       !ny_var_r=block(g)%jrc
       !nx_var_t=block(g)%itc
       !ny_var_t=block(g)%jtc
       !nx_var_tn=block(g)%itn
       !ny_var_tn=block(g)%jtn
       !st_rc_x=block(g)%irc_st
       !en_rc_x=block(g)%irc_en
       !st_rc_y=block(g)%jrc_st
       !en_rc_y=block(g)%jrc_en


         iPt  = 0
         iPt1 = 0
         iPt2 = 0
!         !$acc update self(cell)
!         ALLOCATE(interceptedIndexPtr(ibCellCount,2), fluidIndexPtr(block(g)%fluidCellCount, 2), solidIndexPtr(solidCellCount, 2))
!!$acc loop collapse(2) seq
        !DO  ng = 1,nfl_blk
        !       g=fl_blk(ng)
         DO k = 2, block(g)%nz +1
         DO j = 2, block(g)%ny +1
         DO i = 2, block(g)%nx +1
        !DO j = st_rc_y, en_rc_y
        !DO i = st_rc_x, en_rc_x
                if( block(g)%cell(i,j,k) == 0)then
                block(g)%fluidCellCount=block(g)%fluidCellCount +1
                endif
         end do
         end do
         end do
        if ( block(g)%move_check == 1)then
         DEALLOCATE (block(g)%fluidIndexPtr)
        endif

         ALLOCATE (block(g)%fluidIndexPtr(block(g)%fluidCellCount,3))
        ! END DO

        !DO 30 ng = 1,nfl_blk
        !       g=fl_blk(ng)
                iPt1=0
         DO 30 k = 2, block(g)%nz +1
         DO 30 j = 2, block(g)%ny +1
         DO 30 i = 2, block(g)%nx +1
        !DO 30 j = st_rc_y, en_rc_y
        !DO 30 i = st_rc_x, en_rc_x
          !   IF (cell(i,j).eq.0) THEN
                if( block(g)%cell(i,j,k) == 0)then
              iPt1 = iPt1 + 1
              block(g)%fluidIndexPtr(iPt1, 1) = i
              block(g)%fluidIndexPtr(iPt1, 2) = j
              block(g)%fluidIndexPtr(iPt1, 3) = k
                endif
  30     CONTINUE
         !DO g = fl_blk(1), fl_blk(nfl_blk)
        !DO ng = 1,nfl_blk
        !       g=fl_blk(ng)
          block(g)%redCellCount = 0
          block(g)%blackCellCount  = 0
        !end do

!!$acc enter data copyin(solidIndexPtr,interceptedIndexPtr, fluidIndexPtr)
!!$acc parallel loop present(fluidIndexPtr) reduction(+: block(g)%redCellCount, block(g)%blackCellCount)
         !DO g = fl_blk(1), fl_blk(nfl_blk)
       ! DO  ng = 1,nfl_blk
       !        g=fl_blk(ng)
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
                if( block(g)%cell(i,j,k) == 0)then
            IF (mod(i+j+k,2)==1) THEN
               block(g)%redCellCount = block(g)%redCellCount + 1
            ELSE
               block(g)%blackCellCount = block(g)%blackCellCount + 1
            ENDIF
            ENDIF
         ENDDO
!!$acc end parallel
         ALLOCATE (block(g)%redCellIndexPtr(block(g)%redCellCount,3) ,block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))
       ! END DO
         !DO g = fl_blk(1), fl_blk(nfl_blk)
       ! DO ng = 1,nfl_blk
       !        g=fl_blk(ng)
         ipt1 = 0
         iPt = 0

!!$acc loop seq
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
                !if( block(g)%cell(i,j,k) .eq. 0)then
            IF (mod(i+j+k,2)==1) THEN
               iPt = iPt + 1
               block(g)%redCellIndexPtr(iPt, 1) = i
               block(g)%redCellIndexPtr(iPt, 2) = j
               block(g)%redCellIndexPtr(iPt, 3) = k
            ELSE
               iPt1 = iPt1 + 1
               block(g)%blackCellIndexPtr(iPt1, 1) = i
               block(g)%blackCellIndexPtr(iPt1, 2) = j
               block(g)%blackCellIndexPtr(iPt1, 3) = k
            ENDIF
            !ENDIF
         ENDDO
       !  ENDDO
!!$acc enter data copyin(redCellIndexPtr, blackCellIndexPtr)
         !DO g = fl_blk(1), fl_blk(nfl_blk)
       ! DO ng = 1,nfl_blk
       !        g=fl_blk(ng)
            print*, g,block(g)%fluidCellCount, block(g)%redCellCount, block(g)%blackCellCount
       ! ENDDO
      END SUBROUTINE cellCount_coarse_new
!**************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE cellCount_solid_coarse_mv
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j, k
        !INTEGER (kind = 8),Intent(in) ::  g
        INTEGER (kind = 8) ::  g
        integer (kind=4) :: nx_var,ny_var, nx_var_r,ny_var_r,nx_var_t,ny_var_t,nx_var_tn,ny_var_tn
        integer (kind=4) :: st_rc_x, en_rc_x, st_rc_y, en_rc_y
        g=1

        if ( coarse_flcnt_check == 1)then
      ! nx_var_r=block(g)%irc
         !end do
      ! ny_var_r=block(g)%jrc
      ! nx_var_t=block(g)%itc
      ! ny_var_t=block(g)%jtc
      ! nx_var_tn=block(g)%itn
      ! ny_var_tn=block(g)%jtn
      ! st_rc_x=block(g)%irc_st
      ! en_rc_x=block(g)%irc_en
      ! st_rc_y=block(g)%jrc_st
      ! en_rc_y=block(g)%jrc_en
        block(g)%fluidCellCount=0
         iPt  = 0
         iPt1 = 0
         iPt2 = 0
         DO k = 2, block(g)%nz +1
         DO j = 2, block(g)%ny +1
         DO i = 2, block(g)%nx +1
            IF (block(g)%cell(i,j,k)==0) THEN
                block(g)%fluidCellCount=block(g)%fluidCellCount +1
            ENDIF
         end do
         end do
         end do
         print*,'bef deall'
         DEALLOCATE (block(g)%fluidIndexPtr)
         print*,'aft deall'
         ALLOCATE (block(g)%fluidIndexPtr(block(g)%fluidCellCount,3))
         !print*,'aft all'
         !!$acc update self(cell)
       ! !$acc loop collapse(2) seq
         DO 30 k = 2, block(g)%nz +1
         DO 30 j = 2, block(g)%ny +1
         DO 30 i = 2, block(g)%nx +1
        !DO 30 j = st_rc_y, en_rc_y
        !DO 30 i = st_rc_x, en_rc_x
            !n   = (j-st_rc_y)*nx_var_r + (i-st_rc_x+1)
            n=    i-1  + (block(g)%nx)*(j-2)
            IF (block(g)%cell(i,j,k)==0) THEN
               iPt1 = iPt1 + 1
               block(g)%fluidIndexPtr(iPt1, 1) = i
               block(g)%fluidIndexPtr(iPt1, 2) = j
               block(g)%fluidIndexPtr(iPt1, 3) = k
       !    ELSEIF (block(g)%cell(i,j).eq.1) THEN
       !       iPt2 = iPt2 + 1
       !       block(g)%solidIndexPtr(iPt2, 1) = i
       !       block(g)%solidIndexPtr(iPt2, 2) = j
       !    ELSEIF (block(g)%cell(i,j).eq.2) THEN
       !       iPt = iPt + 1
       !       block(g)%interceptedIndexPtr(iPt, 1) = i
       !       block(g)%interceptedIndexPtr(iPt, 2) = j
            ENDIF
 30      CONTINUE
         block(g)%redCellCount = 0
         block(g)%blackCellCount  = 0

!!$acc enter data copyin(solidIndexPtr,interceptedIndexPtr, fluidIndexPtr)
!!$acc parallel loop present(fluidIndexPtr) reduction(+: block(g)%redCellCount, block(g)%blackCellCount)
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
            IF (mod(i+j+k,2)==1) THEN
               block(g)%redCellCount = block(g)%redCellCount + 1
            ELSE
               block(g)%blackCellCount = block(g)%blackCellCount + 1
            ENDIF
         ENDDO
!!$acc end parallel
         DEALLOCATE (block(g)%redCellIndexPtr ,block(g)%blackCellIndexPtr)
         ALLOCATE (block(g)%redCellIndexPtr(block(g)%redCellCount,3) ,block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))
         ipt1 = 0
         iPt = 0

!!$acc loop seq
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
            IF (mod(i+j+k,2)==1) THEN
               iPt = iPt + 1
               block(g)%redCellIndexPtr(iPt, 1) = i
               block(g)%redCellIndexPtr(iPt, 2) = j
               block(g)%redCellIndexPtr(iPt, 3) = k
            ELSE
               iPt1 = iPt1 + 1
               block(g)%blackCellIndexPtr(iPt1, 1) = i
               block(g)%blackCellIndexPtr(iPt1, 2) = j
               block(g)%blackCellIndexPtr(iPt1, 3) = k
            ENDIF
         ENDDO
!!$acc enter data copyin(redCellIndexPtr, blackCellIndexPtr)
            print*, g, block(g)%fluidCellCount, block(g)%redCellCount, block(g)%blackCellCount
          ENDIF
          coarse_flcnt_check=0

        END SUBROUTINE cellCount_solid_coarse_mv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************

        SUBROUTINE fine_block_cell_mv


        USE global
        IMPLICIT NONE
        INTEGER (KIND=8) :: i, j, k, g, f, factor, a_blk_no, b_blk_no
        integer (kind=4) :: nx_var,ny_var, nx_var_r,ny_var_r,nx_var_t,ny_var_t,nx_var_tn,ny_var_tn
        integer (kind=4) :: st_rc_x, en_rc_x, st_rc_y, en_rc_y


        !block(1)%cell_n=0
        block(1)%cell=0

        DO g=1, intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk
           factor=intfr(g)%b_msh/intfr(g)%a_msh


       !nx_var_r=block(a_blk_no)%irc
       !ny_var_r=block(a_blk_no)%jrc
       !nx_var_t=block(a_blk_no)%itc
       !ny_var_t=block(a_blk_no)%jtc
       !nx_var_tn=block(a_blk_no)%itn
       !ny_var_tn=block(a_blk_no)%jtn
       !st_rc_x=block(a_blk_no)%irc_st
       !en_rc_x=block(a_blk_no)%irc_en
       !st_rc_y=block(a_blk_no)%jrc_st
       !en_rc_y=block(a_blk_no)%jrc_en


       ! DO j = st_rc_y, en_rc_y
       ! DO i = st_rc_x, en_rc_x
         DO k = 2, block(a_blk_no)%nz +1
         DO j = 2, block(a_blk_no)%ny +1
         DO i = 2, block(a_blk_no)%nx +1

        if ( block(a_blk_no)%xp(i) >= intfr(g)%xintf_start+ block(b_blk_no)%cintp*block(a_blk_no)%dx .and. &
             block(a_blk_no)%xp(i) <= intfr(g)%xintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dx   .and. &
             block(a_blk_no)%zp(k) >= intfr(g)%zintf_start+ block(b_blk_no)%cintp*block(a_blk_no)%dz .and. &
             block(a_blk_no)%zp(k) <= intfr(g)%zintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dz   .and. &
             block(a_blk_no)%yp(j) >= intfr(g)%yintf_start + block(b_blk_no)%cintp*block(a_blk_no)%dy .and. &
             block(a_blk_no)%yp(j) <= intfr(g)%yintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dy     )then

                !block(a_blk_no)%cell_n(i,j)=1
                block(a_blk_no)%cell(i,j,k)=1
                !print*,block(a_blk_no)%xp(i),block(a_blk_no)%yp(j),block(a_blk_no)%cell(i,j)

        endif
        ENDDO
        ENDDO
        ENDDO
        ENDDO

    !! ! DO j = st_rc_y, en_rc_y
    !! ! DO i = st_rc_x, en_rc_x
    !!   DO j = 2, block(1)%ny +1
    !!   DO i = 2, block(1)%nx +1
    !!          block(1)%cell(i,j)=block(1)%cell_n(i,j)
    !!  ENDDO
    !!  ENDDO


        END SUBROUTINE
!**************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE block_move_check
        use global
        IMPLICIT NONE
        INTEGER(KIND=8) :: i,j,k,g, a_blk_no, b_blk_no, factor
        REAL(KIND=8) :: margin, ydisp1, xdisp1, zdisp1, mg1
        REAL(KIND=8) :: marginx, marginy, marginz, y_up_lt, y_dw_lt, yval_up, yval_dw, xval_lt, xval_rt

        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        mg1=block(b_blk_no)%cintp*block(1)%dx
        !margin=0.008 + (block(b_blk_no)%cintp*block(1)%dx)
       !marginx=dmin1((0.365-0.300-mg1),(0.490-0.400-mg1))
       !marginy=dmin1((0.600-0.500-mg1),(0.500-0.410-mg1))
       !marginz=dmin1((0.600-0.500-mg1),(0.500-0.410-mg1))
        marginx=4.5*mg1
        marginy=5*mg1
        marginz=1.52*mg1
        factor=intfr(g)%b_msh/intfr(g)%a_msh
        yval_up=dmax1(block(b_blk_no)%ynode1(block(b_blk_no)%mk),block(b_blk_no)%nxty_cent+0.02)
        yval_dw=dmin1(block(b_blk_no)%ynode1(block(b_blk_no)%mk),block(b_blk_no)%nxty_cent-0.02)
        xval_lt=(block(b_blk_no)%xnode1(block(b_blk_no)%mkx1))
        xval_rt=(block(b_blk_no)%xnode1(block(b_blk_no)%mkx2))
        !ydisp1= dmin1(abs(block(b_blk_no)%nxty_cent-y_dw_lt - intfr(g)%yintf_start),abs(block(b_blk_no)%nxty_cent+y_up_lt - intfr(g)%yintf_end))
        ydisp1= dmin1(abs(yval_dw-intfr(g)%yintf_start),abs(yval_up - intfr(g)%yintf_end))
        !xdisp1= dmin1(abs(xval_lt - intfr(g)%xintf_start),abs(xval_rt - intfr(g)%xintf_end))
        xdisp1= (abs(xval_lt - intfr(g)%xintf_start))
        zdisp1= dmin1(abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start),abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end))
        print*,'blk_check_cond', marginx, marginy
        print*,'ydsip',abs(yval_dw-intfr(g)%yintf_start),abs(yval_up - intfr(g)%yintf_end)
        print*,block(b_blk_no)%xnode1(block(b_blk_no)%mkx1),block(b_blk_no)%nxtx_cent
        print*,block(b_blk_no)%ynode1(block(b_blk_no)%mk),block(b_blk_no)%nxty_cent
       ! print*,'xdisp',g,xdisp1
       ! print*,'zdisp',g,zdisp1
        block(b_blk_no)%move_amty=0
        block(b_blk_no)%move_amtx=0
        block(b_blk_no)%move_amtz=0
       !print*,block(b_blk_no)%ynode(block(b_blk_no)%bt_pt),intfr(g)%yintf_start,abs(block(b_blk_no)%ynode(block(b_blk_no)%bt_pt) - intfr(g)%yintf_start)
      ! print*,(abs(block(b_blk_no)%piv_y - intfr(g)%yintf_start) )
      ! print*,(abs(block(b_blk_no)%piv_y - intfr(g)%yintf_end)   )
        if ( (abs(xval_lt - intfr(g)%xintf_start)      <= marginx )  .or. &
             !(abs(xval_rt - intfr(g)%xintf_end)        .le. marginx ) .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start)      <= marginz )  .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end)        <= marginz )  .or. &
             ((abs(yval_dw - intfr(g)%yintf_start) <= marginy ) )   .or. &
             ((abs(yval_up - intfr(g)%yintf_end)  <= marginy ) ) )then
            ! ((abs(block(b_blk_no)%nxty_cent-y_dw_lt - intfr(g)%yintf_start) .le. marginy ) .and. block(g)%ydot .le.0)   .or. &
            ! ((abs(block(b_blk_no)%nxty_cent+y_up_lt - intfr(g)%yintf_end)  .le. marginy ) .and. block(g)%ydot .ge.0) )then

        !if ( (xdisp1 .lt. marginx) .or. (ydisp1 .lt. marginy) .or. (zdisp1 .lt. marginz))then

                block(b_blk_no)%move_check=1
             block(b_blk_no)%blk_mv_tag=1.
                 coarse_flcnt_check=1
             !   if( ydisp1 .lt. marginy)then
            !if( (abs(block(b_blk_no)%nxty_cent - y_dw_lt - intfr(g)%yintf_start) .le. marginy  ) .or. &
            !    (abs(block(b_blk_no)%nxty_cent + y_up_lt - intfr(g)%yintf_end) .le. marginy))then
            if( abs(yval_dw - intfr(g)%yintf_start) <= marginy    .or. &
                (abs(yval_up - intfr(g)%yintf_end)  <= marginy ) )then

                block(b_blk_no)%move_amty=floor((block(b_blk_no)%nxty_cent -block(b_blk_no)%inity_cent)/block(b_blk_no)%dy)
               !print*,block(b_blk_no)%move_amty,block(b_blk_no)%nxty_cent -block(b_blk_no)%inity_cent,block(b_blk_no)%dy
               !print*,block(b_blk_no)%nxty_cent ,block(b_blk_no)%inity_cent
                if ( abs(block(b_blk_no)%move_amty) < factor)then
                        if ( block(b_blk_no)%move_amty < 0)then
                        block(b_blk_no)%move_amty=-factor
                        else
                        block(b_blk_no)%move_amty=factor
                        endif

                else
                        if ( block(b_blk_no)%move_amty < 0)then
                        block(b_blk_no)%move_amty=block(b_blk_no)%move_amty + mod(abs(block(b_blk_no)%move_amty),factor)
                        else
                        block(b_blk_no)%move_amty=block(b_blk_no)%move_amty - mod(abs(block(b_blk_no)%move_amty),factor)
                        endif

                endif
                block(b_blk_no)%inity_cent=block(b_blk_no)%nxty_cent

                endif
                !if( xdisp1 .lt. dmin1((0.390-0.320-mg1),(0.490-0.390-mg1)))then
        if ( (abs(xval_lt- intfr(g)%xintf_start) < marginx ) )then
!             (abs(xval_rt- intfr(g)%xintf_end)   .lt. marginx )) then
                block(b_blk_no)%move_amtx=floor((block(b_blk_no)%nxtx_cent -block(b_blk_no)%initx_cent)/block(b_blk_no)%dx)-factor
                print*,'move_amtx',block(b_blk_no)%move_amtx
                if ( abs(block(b_blk_no)%move_amtx) < factor)then
                        if ( block(b_blk_no)%move_amtx < 0)then
                        block(b_blk_no)%move_amtx=-factor
                        else
                        block(b_blk_no)%move_amtx=factor
                        endif

                else
                        if ( block(b_blk_no)%move_amtx < 0)then
                        block(b_blk_no)%move_amtx=block(b_blk_no)%move_amtx + mod(abs(block(b_blk_no)%move_amtx),factor)
                        else
                        block(b_blk_no)%move_amtx=block(b_blk_no)%move_amtx - mod(abs(block(b_blk_no)%move_amtx),factor)
                        endif

                endif
                block(b_blk_no)%initx_cent=block(b_blk_no)%nxtx_cent
                endif
                !if( zdisp1 .lt. dmin1((0.600-0.500-mg1),(0.500-0.410-mg1)))then
            if( (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start) < (0.500-0.410-mg1))  .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end)   < (0.600-0.500-mg1) ))then
                block(b_blk_no)%move_amtz=floor((block(b_blk_no)%nxtz_cent -block(b_blk_no)%initz_cent)/block(b_blk_no)%dz)
                if ( abs(block(b_blk_no)%move_amtz) < factor)then
                        if ( block(b_blk_no)%move_amtz < 0)then
                        block(b_blk_no)%move_amtz=-factor
                        else
                        block(b_blk_no)%move_amtz=factor
                        endif

                else
                        if ( block(b_blk_no)%move_amtz < 0)then
                        block(b_blk_no)%move_amtz=block(b_blk_no)%move_amtz + mod(abs(block(b_blk_no)%move_amtz),factor)
                        else
                        block(b_blk_no)%move_amtz=block(b_blk_no)%move_amtz - mod(abs(block(b_blk_no)%move_amtz),factor)
                        endif

                endif
                block(b_blk_no)%initz_cent=block(b_blk_no)%nxtz_cent
                endif
                 print*,'blk_movez',block(b_blk_no)%move_amtz
                 print*,'blk_movey',block(b_blk_no)%move_amty
                 print*,'blk_movex',block(b_blk_no)%move_amtx
               !intfr(g)%yintf_start=intfr(g)%yintf_start + block(b_blk_no)%move_amty * block(b_blk_no)%dy
               !intfr(g)%yintf_end=intfr(g)%yintf_end + block(b_blk_no)%move_amty * block(b_blk_no)%dy
       !! endif
!!!     if ( ita .eq. 2) then
!!!             block(b_blk_no)%move_check=1
!!!             block(b_blk_no)%move_amty=-10
!!!             block(b_blk_no)%move_amtx=0
!!!             block(b_blk_no)%blk_mv_tag=1.
!!!     endif

      !!if ( ita .le. 3) then
      !!!if ( mod(ita,10) .eq. 0)then
      !!        block(b_blk_no)%move_amty=32
      !!        block(b_blk_no)%move_amtx=-32
      !!endif
        !dist1=(block(b_blk_no)%ynode(block(b_blk_)%tp_pt) - block(b_blk_no)%yintf_end)**2
        !if (dist1 .lt. 1.5) then

      !!if ( mod(ita,10) .eq. 0)then
      !!         block(b_blk_no)%move_check=1
      !!         coarse_flcnt_check=1
        !if( block(b_blk_no)%move_amty .lt. 0) then
       !if ( block(b_blk_no)%y1(3) .lt. 18)then
       !        block(b_blk_no)%move_amty=50
       !elseif ( block(b_blk_no)%y1(3) .gt. 18)then
       !        block(b_blk_no)%move_amty=-50
       !elseif ( block(b_blk_no)%y1(3) .eq. 18 .and. block(b_blk_no)%move_amty .eq.50)then
       !        block(b_blk_no)%move_amty=50
       !elseif ( block(b_blk_no)%y1(3) .eq. 18 .and. block(b_blk_no)%move_amty .eq. -50)then
       !        block(b_blk_no)%move_amty=-50
       !endif
       !if ( block(b_blk_no)%x1(3) .lt. 8)then
       !        block(b_blk_no)%move_amtx=50
       !elseif ( block(b_blk_no)%x1(3) .gt. 8)then
       !        block(b_blk_no)%move_amtx=-50
       !elseif ( block(b_blk_no)%x1(3) .eq. 8 .and. block(b_blk_no)%move_amtx .eq.50)then
       !        block(b_blk_no)%move_amtx=50
       !elseif ( block(b_blk_no)%x1(3) .eq. 8 .and. block(b_blk_no)%move_amtx .eq. -50)then
       !        block(b_blk_no)%move_amtx=-50
       !endif

       !!    block(b_blk_no)%move_amty = block(b_blk_no)%move_amty * ((-1))
       !!    block(b_blk_no)%move_amtx = block(b_blk_no)%move_amtx * ((-1))
       !!    block(b_blk_no)%blk_mv_tag=1.
             intfr(g)%xintf_st_new=dmax1(intfr(g)%xintf_start,intfr(g)%xintf_start + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx))
             intfr(g)%xintf_en_new=dmin1(intfr(g)%xintf_end,intfr(g)%xintf_end + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx))
             intfr(g)%yintf_st_new=dmax1(intfr(g)%yintf_start,intfr(g)%yintf_start + (block(b_blk_no)%move_amty * block(b_blk_no)%dy))
             intfr(g)%yintf_en_new=dmin1(intfr(g)%yintf_end,intfr(g)%yintf_end + (block(b_blk_no)%move_amty * block(b_blk_no)%dy))
             intfr(g)%zintf_st_new=dmax1(intfr(g)%zintf_start,intfr(g)%zintf_start + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz))
             intfr(g)%zintf_en_new=dmin1(intfr(g)%zintf_end,intfr(g)%zintf_end + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz))
            print*,'zstart_org',intfr(g)%zintf_start,'zend_org',intfr(g)%zintf_end
            print*,'ystart_org',intfr(g)%yintf_start,'yend_org',intfr(g)%yintf_end
            print*,'xstart_org',intfr(g)%xintf_start,'xend_org',intfr(g)%xintf_end

             intfr(g)%xintf_start=intfr(g)%xintf_start + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx)
             intfr(g)%xintf_end=intfr(g)%xintf_end + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx)
             intfr(g)%yintf_start=intfr(g)%yintf_start + (block(b_blk_no)%move_amty * block(b_blk_no)%dy)
             intfr(g)%yintf_end=intfr(g)%yintf_end + (block(b_blk_no)%move_amty * block(b_blk_no)%dy)
             intfr(g)%zintf_start=intfr(g)%zintf_start + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz)
             intfr(g)%zintf_end=intfr(g)%zintf_end + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz)
            print*,'zstart_mv',intfr(g)%zintf_start,'zend_mv',intfr(g)%zintf_end
            print*,'ystart_mv',intfr(g)%yintf_start,'yend_mv',intfr(g)%yintf_end
            print*,'xstart_mv',intfr(g)%xintf_start,'xend_mv',intfr(g)%xintf_end
            print*,'zstart_nw',intfr(g)%zintf_st_new,'zend_nw',intfr(g)%zintf_en_new
            print*,'ystart_nw',intfr(g)%yintf_st_new,'yend_nw',intfr(g)%yintf_en_new
            print*,'xstart_nw',intfr(g)%xintf_st_new,'xend_nw',intfr(g)%xintf_en_new

       !block(b_blk_no)%cell_n=0
       !DO i=1,block(b_blk_no)%nx+2
       !DO j=1,block(b_blk_no)%ny+2
       !   if( block(b_blk_no)%xp(i) .gt. intfr(g)%xintf_st_new .and. &
       !       block(b_blk_no)%xp(i) .lt. intfr(g)%xintf_en_new .and. &
       !       block(b_blk_no)%yp(i) .gt. intfr(g)%yintf_st_new .and. &
       !       block(b_blk_no)%yp(i) .lt. intfr(g)%yintf_en_new )then

       !        block(b_blk_no)%cell_n=1
       !
       !   endif
       !ENDDO
       !ENDDO
        DO i=1,block(b_blk_no)%nx+2
           if( block(b_blk_no)%xp(i) > intfr(g)%xintf_st_new)then
                block(b_blk_no)%cpy_x_start=i
                print*,'xpst',block(b_blk_no)%xp(i),'i',i
                exit
          endif
        ENDDO

        DO i=block(b_blk_no)%cpy_x_start,block(b_blk_no)%nx+2
           if( block(b_blk_no)%xp(i) > intfr(g)%xintf_en_new)then
                block(b_blk_no)%cpy_x_end=i-1
                print*,'xpen',block(b_blk_no)%xp(i-1),'i',i-1
                exit
          endif
        ENDDO
        DO j=1,block(b_blk_no)%ny+2
           if( block(b_blk_no)%yp(j) > intfr(g)%yintf_st_new)then
                block(b_blk_no)%cpy_y_start=j
                print*,'ypst',block(b_blk_no)%yp(j),'j',j
                exit
          endif
        ENDDO

        DO j=block(b_blk_no)%cpy_y_start,block(b_blk_no)%ny+2
           if( block(b_blk_no)%yp(j) > intfr(g)%yintf_en_new)then
                block(b_blk_no)%cpy_y_end=j-1
                print*,'ypen',block(b_blk_no)%yp(j-1),'j',j-1
                exit
          endif
        ENDDO

        DO k=1,block(b_blk_no)%nz+2
           if( block(b_blk_no)%zp(k) > intfr(g)%zintf_st_new)then
                block(b_blk_no)%cpy_z_start=k
                print*,'zpst',block(b_blk_no)%zp(k),'k',k
                exit
          endif
        ENDDO

        DO k=block(b_blk_no)%cpy_z_start,block(b_blk_no)%nz+2
           if( block(b_blk_no)%zp(k) > intfr(g)%zintf_en_new)then
                block(b_blk_no)%cpy_z_end=k-1
                print*,'zpen',block(b_blk_no)%zp(k-1),'k',k-1
                exit
          endif
        ENDDO

        DO k=1,block(b_blk_no)%nz+2
        DO i=1,block(b_blk_no)%nx+2
        DO j=1,block(b_blk_no)%ny+2
            block(b_blk_no)%u_dum(i,j,k)=block(b_blk_no)%u(i,j,k)
            block(b_blk_no)%v_dum(i,j,k)=block(b_blk_no)%v(i,j,k)
            block(b_blk_no)%w_dum(i,j,k)=block(b_blk_no)%w(i,j,k)
            block(b_blk_no)%p_dum(i,j,k)=block(b_blk_no)%p(i,j,k)
        ENDDO
        ENDDO
        ENDDO
        block(b_blk_no)%xp_dum=block(b_blk_no)%xp
        block(b_blk_no)%yp_dum=block(b_blk_no)%yp
        block(b_blk_no)%zp_dum=block(b_blk_no)%zp
        !pause
        endif

        ENDDO
        END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************

        SUBROUTINE change_block_coords

        use global
        IMPLICIT NONE
        INTEGER(KIND=8) :: i,j,k,g, a_blk_no, b_blk_no,countx_st,countz_st,county_st
        REAL(KIND=8) :: change_y_f,change_x_f
        REAL(KIND=8) :: change_z_f
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        DO g=blk_start,nblocks

        if ( block(g)% move_check == 1) then


       change_z_f= block(g)%move_amtz * block(g)%dz
       change_y_f= block(g)%move_amty * block(g)%dy
       change_x_f= block(g)%move_amtx * block(g)%dx
        print*, 'chz',block(g)%move_amtz , block(g)%dz
        print*, 'chy',block(g)%move_amty , block(g)%dy
        print*, 'chx',block(g)%move_amtx , block(g)%dx
       !block(g)%yshift_move=yshift+change_y_f
       !block(g)%xshift_move=xshift+change_x_f

        DO i = 1, block(g)%nx+3
           block(g)%x1(i) = block(g)%x1(i) + change_x_f
        ENDDO

        DO i = 1, block(g)%ny+3
                  ! print*,'y1p',block(g)%y1(i)
           block(g)%y1(i) = block(g)%y1(i) + change_y_f
                  ! print*,'y1a',block(g)%y1(i)
        ENDDO
       ! pause
	 DO i = 1, block(g)%nz+3
           block(g)%z1(i) = block(g)%z1(i)+ change_z_f
        ENDDO

        DO i = 1, block(g)%nx+3
           block(g)%xu(i) = block(g)%x1(i)
        ENDDO

        DO i = 1, block(g)%ny+3
           block(g)%yv(i) = block(g)%y1(i)
        ENDDO

	 DO i = 1, block(g)%nz+3
           block(g)%zw(i) = block(g)%z1(i)
        ENDDO

        DO i = 1, block(g)%ny+2
           block(g)%yu(i) = 0.5_rk*(block(g)%y1(i)+block(g)%y1(i+1))
	   block(g)%yw(i) = block(g)%yu(i)
           block(g)%yp(i) = block(g)%yu(i)
        END DO

        DO i = 1, block(g)%nx+2
           block(g)%xv(i) = 0.5_rk*(block(g)%x1(i)+block(g)%x1(i+1))
           block(g)%xw(i) = block(g)%xv(i)
           block(g)%xp(i) = block(g)%xv(i)
        END DO

	 DO i = 1, block(g)%nz+2
           block(g)%zu(i) = 0.5_rk*(block(g)%z1(i)+block(g)%z1(i+1))
	   block(g)%zv(i) = block(g)%zu(i)
           block(g)%zp(i) = block(g)%zu(i)
        END DO
        DO k=1,block(g)%nz+1
        DO j=1,block(g)%ny+1
        DO i=1,block(g)%nx+1

        block(g)%xpn1(i,j,k)=block(g)%x1(i)
        block(g)%ypn1(i,j,k)=block(g)%y1(j)
        block(g)%zpn1(i,j,k)=block(g)%z1(k)


        END DO
        END DO
        END DO

        DO k=2,block(g)%nz+1
        DO j=2,block(g)%ny+1
        DO i=2,block(g)%nx+1

        block(g)%xp1(i,j,k)=block(g)%xp(i)
        block(g)%yp1(i,j,k)=block(g)%yp(j)
        block(g)%zp1(i,j,k)=block(g)%zp(k)


        END DO
        END DO
        END DO

        ENDIF

        ENDDO
      ! block(g)%u_dum=0
      ! block(g)%v_dum=0
      ! block(g)%p_dum=0

      ! DO i=1,block(g)%nx+2
      ! DO j=1,block(g)%ny+2
      !     if (block(1)%cell_n(i,j) .eq. 1 .and. block(1)%cell(i,j) .eq. 1)then
      !     i_new=i+block(g)%move_amtx
      !     j_new=j+block(g)%move_amty
      !     block(g)%u_dum(i,j)=block(g)%u(i_new,j_new)
      !     block(g)%v_dum(i,j)=block(g)%v(i_new,j_new)
      !     block(g)%p_dum(i,j)=block(g)%p(i_new,j_new)
      !     endif
      ! ENDDO
      ! ENDDO

!!      DO i = block(g)%itc_st, block(g)%itc_en
!!         block(g)%xu(i) = 0.5_rk*(block(g)%x1(i)+block(g)%x1(i+1))
!!         block(g)%xp(i) = block(g)%xu(i)
!!             ! print*,i,block(g)%y1(i), block(g)%yu(i)
!!      END DO
!!      DO i = block(g)%itn_st, block(g)%itn_en
!!         block(g)%xv(i) = block(g)%x1(i)
!!      ENDDO

!!!!y

!!      DO i = block(g)%jtc_st, block(g)%jtc_en
!!         block(g)%yu(i) = 0.5_rk*(block(g)%y1(i)+block(g)%y1(i+1))
!!         block(g)%yp(i) = block(g)%yu(i)
!!             ! print*,i,block(g)%y1(i), block(g)%yu(i)
!!      END DO
!!      DO i = block(g)%jtn_st, block(g)%jtn_en
!!         block(g)%yv(i) = block(g)%y1(i)
!!      ENDDO
        !block(g)%move_check=0

!!!!search
      !! DO i = block(g)%irn_st, block(g)%irn_en
      !!       if (block(g)%x1(i) .gt. 9.5)then
      !!        block(g)%i_startSearch= i-10
      !!        exit
      !!       endif
      !! ENDDO
      !! DO i = block(g)%i_startSearch, block(g)%irn_en
      !!       if (block(g)%x1(i) .gt. 10.5)then
      !!        block(g)%i_endSearch= i+10
      !!        exit
      !!       endif
      !! ENDDO

      !! DO i = block(g)%jrn_st, block(g)%jrn_en
      !!       if (block(g)%y1(i) .gt. 19.5)then
      !!        block(g)%j_startSearch= i-10
      !!        exit
      !!       endif
      !! ENDDO
      !! DO i = block(g)%j_startSearch, block(g)%jrn_en
      !!       if (block(g)%y1(i) .gt. 20.5)then
      !!        block(g)%j_endSearch= i+10
      !!        exit
      !!       endif
      !! ENDDO
       ! print*, block(g)%y1(1), change_y_f
     !     DEALLOCATE(block(g)%xcent, block(g)%ycent, block(g)%cosAlpha, block(g)%cosBeta, block(g)%area)
     !     CALL computeSurfaceVariables_move(g)
     !     CALL computeSurfaceNorm_move(g)
        !print*, 'before tag mv'
        !call tagging_block_move(g)
!!          DEALLOCATE(block(g)%interceptedIndexPtr, block(g)%pNormDis,block(g)%nelp, &
!!             block(g)%nelu1, block(g)%nelu2, block(g)%nelv1, block(g)%nelv2, block(g)%u1NormDis, block(g)%u2NormDis, &
!!             block(g)%v1NormDis, block(g)%v2NormDis,block(g)%solidIndexPtr)
!!     !print*,'After deallocate'
!!         !DO i=1,nblocks
!!              DEALLOCATE(block(g)%fluidIndexPtr,block(g)%redCellIndexPtr,block(g)%blackCellIndexPtr)
!!       CALL cellCount_solid_move(g)
!!      call computeNormDistance_move(g)

        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        if ( block(b_blk_no)% move_check == 1) then
        block(b_blk_no)%cell_n=0
        DO k=1,block(b_blk_no)%nz+2
        DO j=1,block(b_blk_no)%ny+2
        DO i=1,block(b_blk_no)%nx+2
        if( block(b_blk_no)%xp(i) > (intfr(g)%xintf_st_new  + block(b_blk_no)%dx) .and. &
               block(b_blk_no)%xp(i) < (intfr(g)%xintf_en_new  - block(b_blk_no)%dx).and. &
               block(b_blk_no)%zp(k) > (intfr(g)%zintf_st_new  + block(b_blk_no)%dz) .and. &
               block(b_blk_no)%zp(k) < (intfr(g)%zintf_en_new  - block(b_blk_no)%dz).and. &
               block(b_blk_no)%yp(j) > (intfr(g)%yintf_st_new  + block(b_blk_no)%dy).and. &
               block(b_blk_no)%yp(j) < (intfr(g)%yintf_en_new  - block(b_blk_no)%dy))then

                block(b_blk_no)%cell_n(i,j,k)=1

          !else
          !     print*,'not',i,j
           endif
        ENDDO
        ENDDO
        ENDDO

        DO i=1,block(b_blk_no)%nx+2
           if( block(b_blk_no)%xp(i) > intfr(g)%xintf_st_new)then
                block(b_blk_no)%cpy_x_start_mv=i
                print*,'xpst_f',block(b_blk_no)%xp(i),'i',i
                exit
          endif
        ENDDO

        DO i=block(b_blk_no)%cpy_x_start_mv,block(b_blk_no)%nx+2
           if( block(b_blk_no)%xp(i) > intfr(g)%xintf_en_new)then
                block(b_blk_no)%cpy_x_end_mv=i-1
                print*,'xpen_f',block(b_blk_no)%xp(i-1),'i',i-1
                exit
          endif
        ENDDO
        DO j=1,block(b_blk_no)%ny+2
           if( block(b_blk_no)%yp(j) > intfr(g)%yintf_st_new)then
                block(b_blk_no)%cpy_y_start_mv=j
                print*,'ypst_f',block(b_blk_no)%yp(j),'j',j
                exit
          endif
        ENDDO

        DO j=block(b_blk_no)%cpy_y_start_mv,block(b_blk_no)%ny+2
           if( block(b_blk_no)%yp(j) > intfr(g)%yintf_en_new)then
                block(b_blk_no)%cpy_y_end_mv=j-1
                print*,'ypen_f',block(b_blk_no)%yp(j-1),'j',j-1
                exit
          endif
        ENDDO
        DO j=1,block(b_blk_no)%nz+2
           if( block(b_blk_no)%zp(j) > intfr(g)%zintf_st_new)then
                block(b_blk_no)%cpy_z_start_mv=j
                print*,'zpst_f',block(b_blk_no)%zp(j),'k',j
                exit
          endif
        ENDDO

        DO j=block(b_blk_no)%cpy_z_start_mv,block(b_blk_no)%nz+2
           if( block(b_blk_no)%zp(j) > intfr(g)%zintf_en_new)then
                block(b_blk_no)%cpy_z_end_mv=j-1
                print*,'zpen_f',block(b_blk_no)%zp(j-1),'k',j-1
                exit
          endif
        ENDDO
        print*,'bef'
        countx_st=block(b_blk_no)%cpy_x_start
        county_st=block(b_blk_no)%cpy_y_start
        countz_st=block(b_blk_no)%cpy_z_start
       !print*,block(b_blk_no)%cpy_x_start_mv,block(b_blk_no)%cpy_x_end_mv
       !print*,block(b_blk_no)%cpy_y_start_mv,block(b_blk_no)%cpy_y_end_mv
        OPEN(UNIT=12,FILE='log.dat',STATUS='unknown',access='append')
 11     FORMAT(2F13.5)
 111     FORMAT(6I6)
        block(b_blk_no)%u=0
        block(b_blk_no)%v=0
        block(b_blk_no)%w=0
        block(b_blk_no)%p=0
        DO k=block(b_blk_no)%cpy_z_start_mv,block(b_blk_no)%cpy_z_end_mv
        countx_st=block(b_blk_no)%cpy_x_start
        DO i=block(b_blk_no)%cpy_x_start_mv,block(b_blk_no)%cpy_x_end_mv
        county_st=block(b_blk_no)%cpy_y_start
        DO j=block(b_blk_no)%cpy_y_start_mv,block(b_blk_no)%cpy_y_end_mv
                !print*,'inside'
                block(b_blk_no)%u(i,j,k)=block(b_blk_no)%u_dum(countx_st,county_st,countz_st)
                block(b_blk_no)%v(i,j,k)=block(b_blk_no)%v_dum(countx_st,county_st,countz_st)
                block(b_blk_no)%w(i,j,k)=block(b_blk_no)%w_dum(countx_st,county_st,countz_st)
                block(b_blk_no)%p(i,j,k)=block(b_blk_no)%p_dum(countx_st,county_st,countz_st)
!               write(12,111)i,countx_st,j,county_st,k,countz_st
!               write(12,11)block(b_blk_no)%xp(i),block(b_blk_no)%xp_dum(countx_st)
!               write(12,11)block(b_blk_no)%yp(j),block(b_blk_no)%yp_dum(county_st)
!               write(12,11)block(b_blk_no)%zp(k),block(b_blk_no)%zp_dum(countz_st)
                county_st=county_st+1
        ENDDO
                countx_st=countx_st+1
        ENDDO
                countz_st=countz_st+1
        ENDDO
        print*,'aft'
        close(12)
       ! pause
        ENDIF
        ENDDO



        END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE change_block_interface

        use global
        IMPLICIT NONE
        INTEGER(KIND=8) :: i,j,k,g, a_blk_no, b_blk_no, factor
        REAL(KIND=8) :: change_y_f
        INTEGER, PARAMETER :: rk = selected_real_kind(8)



        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
         factor=intfr(g)%b_msh/intfr(g)%a_msh
        if ( block(b_blk_no)% move_check == 1) then
        DO j=1,intfr(g)%counterxp

        intfr(g)%px_interface_det(1,j)= intfr(g)%px_interface_det(1,j) + (block(b_blk_no)%move_amtx/factor)
        !print*,j,intfr(g)%py_interface_det(1,j)
        ENDDO
        DO j=1,intfr(g)%counterxu

        intfr(g)%ux_interface_det(1,j)= intfr(g)%ux_interface_det(1,j) + (block(b_blk_no)%move_amtx/factor)

        ENDDO
        DO j=1,intfr(g)%counterxv

        intfr(g)%vx_interface_det(1,j)= intfr(g)%vx_interface_det(1,j) + (block(b_blk_no)%move_amtx/factor)

        ENDDO
        DO j=1,intfr(g)%counterxw

        intfr(g)%wx_interface_det(1,j)= intfr(g)%wx_interface_det(1,j) + (block(b_blk_no)%move_amtx/factor)

        ENDDO

        DO j=1,intfr(g)%counteryp

        intfr(g)%py_interface_det(1,j)= intfr(g)%py_interface_det(1,j) + (block(b_blk_no)%move_amty/factor)
        !print*,j,intfr(g)%py_interface_det(1,j)
        ENDDO
        DO j=1,intfr(g)%counteryu

        intfr(g)%uy_interface_det(1,j)= intfr(g)%uy_interface_det(1,j) + (block(b_blk_no)%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counteryv

        intfr(g)%vy_interface_det(1,j)= intfr(g)%vy_interface_det(1,j) + (block(b_blk_no)%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counteryw

        intfr(g)%wy_interface_det(1,j)= intfr(g)%wy_interface_det(1,j) + (block(b_blk_no)%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counterzp

        intfr(g)%pz_interface_det(1,j)= intfr(g)%pz_interface_det(1,j) + (block(b_blk_no)%move_amtz/factor)
        !print*,j,intfr(g)%py_interface_det(1,j)
        ENDDO
        DO j=1,intfr(g)%counterzu

        intfr(g)%uz_interface_det(1,j)= intfr(g)%uz_interface_det(1,j) + (block(b_blk_no)%move_amtz/factor)

        ENDDO
        DO j=1,intfr(g)%counterzv

        intfr(g)%vz_interface_det(1,j)= intfr(g)%vz_interface_det(1,j) + (block(b_blk_no)%move_amtz/factor)

        ENDDO
        DO j=1,intfr(g)%counterzw

        intfr(g)%wz_interface_det(1,j)= intfr(g)%wz_interface_det(1,j) + (block(b_blk_no)%move_amtz/factor)

        ENDDO
      !!call fineUpdate
      !!call fineUpdate_p
      !!call solidCellBC_move(b_blk_no)
        !block(b_blk_no)%move_check=0
        DO j=1,intflines
        print*,'**********************px***************************'
        DO i=1,intfr(j)%counterxp
        WRITE(*,33)'px',i,intfr(j)%px_interface_det(1,i),intfr(j)%px_interface_det(2,i),intfr(j)%px_interface_det(3,i),block(a_blk_no)%xp(intfr(j)%px_interface_det(1,i)),block(b_blk_no)%xp(intfr(j)%px_interface_det(2,i)),block(b_blk_no)%xp(intfr(j)%px_interface_det(3,i))
 33       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************py***************************'
        DO i=1,intfr(j)%counteryp
        WRITE(*,133)'py',i,intfr(j)%py_interface_det(1,i),intfr(j)%py_interface_det(2,i),intfr(j)%py_interface_det(3,i),block(a_blk_no)%yp(intfr(j)%py_interface_det(1,i)),block(b_blk_no)%yp(intfr(j)%py_interface_det(2,i)),block(b_blk_no)%yp(intfr(j)%py_interface_det(3,i))
 133       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intflines
        print*,'**********************pz***************************'
        DO i=1,intfr(j)%counterzp
        WRITE(*,933)'pz',i,intfr(j)%pz_interface_det(1,i),intfr(j)%pz_interface_det(2,i),intfr(j)%pz_interface_det(3,i),block(a_blk_no)%zp(intfr(j)%pz_interface_det(1,i)),block(b_blk_no)%zp(intfr(j)%pz_interface_det(2,i)),block(b_blk_no)%zp(intfr(j)%pz_interface_det(3,i))
 933       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do

        DO j=1,intflines
        print*,'**********************ux***************************'
        DO i=1,intfr(j)%counterxu
        WRITE(*,331)'ux',i,intfr(j)%ux_interface_det(1,i),intfr(j)%ux_interface_det(2,i),intfr(j)%ux_interface_det(3,i)  ,block(a_blk_no)%xu(intfr(j)%ux_interface_det(1,i)),block(b_blk_no)%xu(intfr(j)%ux_interface_det(2,i)),block(b_blk_no)%xu(intfr(j)%ux_interface_det(3,i))
 331       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************uy***************************'
        DO i=1,intfr(j)%counteryu
        WRITE(*,1331)'uy',i,intfr(j)%uy_interface_det(1,i),intfr(j)%uy_interface_det(2,i),intfr(j)%uy_interface_det(3,i) ,block(a_blk_no)%yu(intfr(j)%uy_interface_det(1,i)),block(b_blk_no)%yu(intfr(j)%uy_interface_det(2,i)),block(b_blk_no)%yu(intfr(j)%uy_interface_det(3,i))
 1331       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************uz***************************'
        DO i=1,intfr(j)%counterzu
        WRITE(*,91331)'uz',i,intfr(j)%uz_interface_det(1,i),intfr(j)%uz_interface_det(2,i),intfr(j)%uz_interface_det(3,i) ,block(a_blk_no)%zu(intfr(j)%uz_interface_det(1,i)),block(b_blk_no)%zu(intfr(j)%uz_interface_det(2,i)),block(b_blk_no)%zu(intfr(j)%uz_interface_det(3,i))
91331       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do

        DO j=1,intflines
        print*,'**********************vx***************************'
        DO i=1,intfr(j)%counterxv
        WRITE(*,332)'vx',i,intfr(j)%vx_interface_det(1,i),intfr(j)%vx_interface_det(2,i),intfr(j)%vx_interface_det(3,i)  ,block(a_blk_no)%xv(intfr(j)%vx_interface_det(1,i)),block(b_blk_no)%xv(intfr(j)%vx_interface_det(2,i)),block(b_blk_no)%xv(intfr(j)%vx_interface_det(3,i))
 332       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************vy***************************'
        DO i=1,intfr(j)%counteryv
        WRITE(*,1332)'vy',i,intfr(j)%vy_interface_det(1,i),intfr(j)%vy_interface_det(2,i),intfr(j)%vy_interface_det(3,i) ,block(a_blk_no)%yv(intfr(j)%vy_interface_det(1,i)),block(b_blk_no)%yv(intfr(j)%vy_interface_det(2,i)),block(b_blk_no)%yv(intfr(j)%vy_interface_det(3,i))
 1332       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************vz***************************'
        DO i=1,intfr(j)%counterzv
        WRITE(*,91332)'vz',i,intfr(j)%vz_interface_det(1,i),intfr(j)%vz_interface_det(2,i),intfr(j)%vz_interface_det(3,i) ,block(a_blk_no)%zv(intfr(j)%vz_interface_det(1,i)),block(b_blk_no)%zv(intfr(j)%vz_interface_det(2,i)),block(b_blk_no)%zv(intfr(j)%vz_interface_det(3,i))
91332       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do


        DO j=1,intflines
        print*,'**********************wx***************************'
        DO i=1,intfr(j)%counterxw
        WRITE(*,3329)'wx',i,intfr(j)%wx_interface_det(1,i),intfr(j)%wx_interface_det(2,i),intfr(j)%wx_interface_det(3,i),block(a_blk_no)%xw(intfr(j)%wx_interface_det(1,i)),block(b_blk_no)%xw(intfr(j)%wx_interface_det(2,i)),block(b_blk_no)%xw(intfr(j)%wx_interface_det(3,i))
 3329       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************wy***************************'
        DO i=1,intfr(j)%counteryw
        WRITE(*,13329)'wy',i,intfr(j)%wy_interface_det(1,i),intfr(j)%wy_interface_det(2,i),intfr(j)%wy_interface_det(3,i),block(a_blk_no)%yw(intfr(j)%wy_interface_det(1,i)),block(b_blk_no)%yw(intfr(j)%wy_interface_det(2,i)),block(b_blk_no)%yw(intfr(j)%wy_interface_det(3,i))
13329       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************wz***************************'
        DO i=1,intfr(j)%counterzw
        WRITE(*,93329)'wz',i,intfr(j)%wz_interface_det(1,i),intfr(j)%wz_interface_det(2,i),intfr(j)%wz_interface_det(3,i),block(a_blk_no)%zw(intfr(j)%wz_interface_det(1,i)),block(b_blk_no)%zw(intfr(j)%wz_interface_det(2,i)),block(b_blk_no)%zw(intfr(j)%wz_interface_det(3,i))
93329       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        endif
        ENDDO

        END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************

     SUBROUTINE tagging_block_move
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, m, i, j, k, n1, n2, n3, n4, nel2n, g, ibElems_cnt

        INTEGER            :: iPt, iPt1, a_blk_no, b_blk_no
        REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, &
                              n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
                              cent_x, cent_y, cent_z
        CHARACTER(LEN=100) :: cLine
        CHARACTER(len=150)  :: filename1
        DO g=blk_start,nblocks
        if ( block(g)%move_check == 1)then
	 ALLOCATE(block(g)%minElemcell(block(g)%nx+2,block(g)%ny+2,block(g)%nz+2))


        block(g)%ibCellCount = 0
        block(g)%fluidCellCount = 0
        block(g)%solidCellCount = 0
        block(g)%cell = 0
        block(g)%minElemcell = 0
        n1dotn = 0
        n2dotn = 0
        n3dotn = 0
        n4dotn = 0
        n5dotn = 0
        n6dotn = 0
        n7dotn = 0
        n8dotn = 0
        n9dotn = 0
      !  GOTO 1010
       !!$acc parallel loop gang vector collapse(3)        &
       !!$acc present (cosAlpha, cosBeta, cosGamma, xcent, ycent, zcent, cell, minElemcell, x1, y1, z1, xp, yp, zp)       &
       !!$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, nel2n, m, cent_x, cent_y, cent_z,         &
       !!$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn)
      !$acc parallel loop gang vector collapse(3)        &
      !$acc default(present)       &
      !$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, nel2n, m, cent_x, cent_y, cent_z,         &
      !$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn)

        DO 10 k = block(g)%k_startSearch, block(g)%k_endSearch
        DO 10 j = block(g)%j_startSearch, block(g)%j_endSearch
        DO 10 i = block(g)%i_startSearch, block(g)%i_endSearch

           m = 0
           minDis = 1e14

           n1x = block(g)%xp(i)
           n1y = block(g)%yp(j)
           n1z = block(g)%zp(k)

           n2x = block(g)%x1(i)
           n3x = block(g)%x1(i+1)

           n2y = block(g)%y1(j)
           n3y = block(g)%y1(j+1)

           n2z = block(g)%z1(k)
           n3z = block(g)%z1(k+1)

           !ibElems_cnt=block(g)%ibElemCnt


           !$acc loop seq
           !DO m = 1, ibElems_cnt
           DO m = 1, block(g)%ibElems
              cent_x = block(g)%xcent(m)
              cent_y = block(g)%ycent(m)
              cent_z = block(g)%zcent(m)
              dis  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2)
              IF (dis<minDis) THEN
                 minDis   = dis
                 nel2n    = m
              ENDIF
           ENDDO
           block(g)%minElemcell(i,j,k) = nel2n

             !n1dotn  = (n1x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                     !(n1y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                     !(n1z - block(g)%zcent(nel2n))*cosGamma(nel2n)

             n2dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

             n3dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

             n4dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

             n5dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

            n6dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

             n7dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

             n8dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

             n9dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
                      (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
                      (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

           !n       = i-1  + nx*(j-2)  + block(g)%nx*ny*(k-2)

           IF ( n2dotn<=-1e-16 .AND. n3dotn<=-1e-16 .AND. n4dotn<=-1e-16 .AND. n5dotn<=-1e-16  &
               .AND. n6dotn<=-1e-16 .AND. n7dotn<=-1e-16 .AND. n8dotn<=-1e-16 .AND. n9dotn<=-1e-16) THEN
        	    block(g)%cell(i,j,k) = 1
           ELSEIF (n2dotn>-1e-16 .AND. n3dotn>-1e-16 .AND. n4dotn>-1e-16 .AND. n5dotn>-1e-16 &
        	    .AND. n6dotn>-1e-16 .AND. n7dotn>-1e-16 .AND. n8dotn>-1e-16 .AND. n9dotn>-1e-16) THEN
                  block(g)%cell(i,j,k) = 0
           ELSE
                  block(g)%cell(i,j,k) = 2

           ENDIF

 10      CONTINUE

         !$acc end parallel
         !!$acc update host(cell)

         block(g)%ibCellCount = 0
         block(g)%solidCellCount = 0
         block(g)%fluidCellCount = 0

        !!$acc parallel loop gang vector collapse(3) default(present) reduction(+: solidCellCount, fluidCellCount, ibCellCount)
         DO 20 k = 2, block(g)%nz+1
         DO 20 j = 2, block(g)%ny+1
         DO 20 i = 2, block(g)%nx+1
            !n       = i-1  + nx*(j-2)  + block(g)%nx*block(g)%block(g)%ny*(k-2)
            IF (block(g)%cell(i,j,k)==1) THEN
               block(g)%solidCellCount = block(g)%solidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==0) THEN
               block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==2) THEN
               block(g)%ibCellCount = block(g)%ibCellCount + 1
            ENDIF
 20      CONTINUE
        !!$acc end parallel

      GOTO 1000
        !1010 open(82,file='inter.dat',status='unknown')
        !  1010  continue
       open(82,file='inter.dat',form='formatted')
      !write(82,*)'variables = "x", "y","z", "var"'
        read(82,*) cLine
        do k = 2, block(g)%nz+1
       do j = 2, block(g)%ny+1
       do i = 2, block(g)%nx+1
        !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
       if(block(g)%cell(i,j,k)==2)then
      ! write(82,*) block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 2
       read(82,*) block(g)%xp(i),block(g)%yp(j),block(g)%zp(k)
       read(82,*) block(g)%cell(i,j,k)
               block(g)%ibCellCount = block(g)%ibCellCount + 1
       endif
       end do
       end do
        enddo
      close(82)

      !open(83,file='fluid.dat',status='unknown')
      open(83,file='fluid.dat',form='formatted')
      !write(83,*)'variables = "x", "y","z","var"'
        read(83,*)cLine
       do k = 2, block(g)%nz+1
       do j = 2, block(g)%ny+1
       do i = 2, block(g)%nx+1
        !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
       if(block(g)%cell(i,j,k)==0)then
      ! write(83,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 0
       read(83,*) block(g)%xp(i),block(g)%yp(j),block(g)%zp(k)
       read(83,*) block(g)%cell(i,j,k)
               block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
       endif
       end do
       end do
        enddo
      close(83)

      !open(84,file='solid.dat',status='unknown')
      open(84,file='solid.dat',form='formatted')
     ! write(84,*)'variables = "x", "y","z","var"'
        read(84,*)cLine
       do k = 2, block(g)%nz+1
       do j = 2, block(g)%ny+1
       do i = 2, block(g)%nx+1
        !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
       if(block(g)%cell(i,j,k)==1)then
      ! write(84,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 1
       read(84,*) block(g)%xp(i),block(g)%yp(j),block(g)%zp(k)
       read(84,*) block(g)%cell(i,j,k)
               block(g)%solidCellCount = block(g)%solidCellCount + 1
       endif
       end do
       end do
        enddo
      close(84)

 1000   CONTINUE
        print*, 'block move tag', block(g)%ibCellCount, block(g)%fluidCellCount, block(g)%solidCellCount
        !block(g)%move_check=0
        ENDIF
        END DO
        !block(g)%move_check=0
        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

        if ( block(b_blk_no)%move_check == 1)then
        !print*,'inside tag_mv'
       call fineUpdate_mv(g)
       call fineUpdate_bd_mv(g)
       !!print*,b_blk_no,'inside_tag_mv'
       !call solidCellBC_move(b_blk_no)
       !call updateVelocity_newv(b_blk_no)
       ! call writeOutput
       ! pause
        endif
        ENDDO




     END SUBROUTINE tagging_block_move

!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE cellCount_solid_coarse
        USE global
        IMPLICIT NONE
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j, k
        !INTEGER (kind = 8),Intent(in) ::  g
        INTEGER (kind = 8) ::  g
        integer (kind=4) :: nx_var,ny_var, nx_var_r,ny_var_r,nx_var_t,ny_var_t,nx_var_tn,ny_var_tn
        integer (kind=4) :: st_rc_x, en_rc_x, st_rc_y, en_rc_y
        g=1

!        if ( coarse_flcnt_check .eq. 1)then
      ! nx_var_r=block(g)%irc
         !end do
      ! ny_var_r=block(g)%jrc
      ! nx_var_t=block(g)%itc
      ! ny_var_t=block(g)%jtc
      ! nx_var_tn=block(g)%itn
      ! ny_var_tn=block(g)%jtn
      ! st_rc_x=block(g)%irc_st
      ! en_rc_x=block(g)%irc_en
      ! st_rc_y=block(g)%jrc_st
      ! en_rc_y=block(g)%jrc_en
        block(g)%fluidCellCount=0
         iPt  = 0
         iPt1 = 0
         iPt2 = 0
         DO k = 2, block(g)%nz +1
         DO j = 2, block(g)%ny +1
         DO i = 2, block(g)%nx +1
            IF (block(g)%cell(i,j,k)==0) THEN
                block(g)%fluidCellCount=block(g)%fluidCellCount +1
            ENDIF
         end do
         end do
         end do
         ALLOCATE (block(g)%fluidIndexPtr(block(g)%fluidCellCount,3))
         !print*,'aft all'
         !!$acc update self(cell)
       ! !$acc loop collapse(2) seq
         DO 30 k = 2, block(g)%nz +1
         DO 30 j = 2, block(g)%ny +1
         DO 30 i = 2, block(g)%nx +1
        !DO 30 j = st_rc_y, en_rc_y
        !DO 30 i = st_rc_x, en_rc_x
            !n   = (j-st_rc_y)*nx_var_r + (i-st_rc_x+1)
           ! n=    i-1  + (block(g)%nx)*(j-2)
            IF (block(g)%cell(i,j,k)==0) THEN
               iPt1 = iPt1 + 1
               block(g)%fluidIndexPtr(iPt1, 1) = i
               block(g)%fluidIndexPtr(iPt1, 2) = j
               block(g)%fluidIndexPtr(iPt1, 3) = k
       !    ELSEIF (block(g)%cell(i,j).eq.1) THEN
       !       iPt2 = iPt2 + 1
       !       block(g)%solidIndexPtr(iPt2, 1) = i
       !       block(g)%solidIndexPtr(iPt2, 2) = j
       !    ELSEIF (block(g)%cell(i,j).eq.2) THEN
       !       iPt = iPt + 1
       !       block(g)%interceptedIndexPtr(iPt, 1) = i
       !       block(g)%interceptedIndexPtr(iPt, 2) = j
            ENDIF
 30      CONTINUE
         block(g)%redCellCount = 0
         block(g)%blackCellCount  = 0

!!$acc enter data copyin(solidIndexPtr,interceptedIndexPtr, fluidIndexPtr)
!!$acc parallel loop present(fluidIndexPtr) reduction(+: block(g)%redCellCount, block(g)%blackCellCount)
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
            IF (mod(i+j+k,2)==1) THEN
               block(g)%redCellCount = block(g)%redCellCount + 1
            ELSE
               block(g)%blackCellCount = block(g)%blackCellCount + 1
            ENDIF
         ENDDO
!!$acc end parallel
         !DEALLOCATE (block(g)%redCellIndexPtr ,block(g)%blackCellIndexPtr)
         ALLOCATE (block(g)%redCellIndexPtr(block(g)%redCellCount,3) ,block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))
         ipt1 = 0
         iPt = 0

!!$acc loop seq
         DO n = 1, block(g)%fluidCellCount
            i = block(g)%fluidIndexPtr(n, 1)
            j = block(g)%fluidIndexPtr(n, 2)
            k = block(g)%fluidIndexPtr(n, 3)
            IF (mod(i+j+k,2)==1) THEN
               iPt = iPt + 1
               block(g)%redCellIndexPtr(iPt, 1) = i
               block(g)%redCellIndexPtr(iPt, 2) = j
               block(g)%redCellIndexPtr(iPt, 3) = k
            ELSE
               iPt1 = iPt1 + 1
               block(g)%blackCellIndexPtr(iPt1, 1) = i
               block(g)%blackCellIndexPtr(iPt1, 2) = j
               block(g)%blackCellIndexPtr(iPt1, 3) = k
            ENDIF
         ENDDO
!!$acc enter data copyin(redCellIndexPtr, blackCellIndexPtr)
            print*, g, block(g)%fluidCellCount, block(g)%redCellCount, block(g)%blackCellCount
         ! ENDIF
         ! coarse_flcnt_check=0

        END SUBROUTINE cellCount_solid_coarse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!******************************************************************

!    SUBROUTINE tagging
!       USE global
!       IMPLICIT NONE
!       INTEGER, PARAMETER :: rk = selected_real_kind(8)
!       INTEGER (kind = 8) ::  n, m, i, j, k, n1, n2, n3, n4, nel2n, g, ibElems_cnt
!
!       INTEGER            :: iPt, iPt1
!       REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, &
!                             n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
!                             cent_x, cent_y, cent_z
!       CHARACTER(LEN=100) :: cLine
!       CHARACTER*150  filename1
!       DO g=blk_start,nblocks
!       ALLOCATE(block(g)%minElemcell(block(g)%nx+2,block(g)%ny+2,block(g)%nz+2))
!
!       !print*,'inside tagging'
!
!       block(g)%ibCellCount = 0
!       block(g)%fluidCellCount = 0
!       block(g)%solidCellCount = 0
!       block(g)%cell = 0
!       block(g)%minElemcell = 0
!       n1dotn = 0
!       n2dotn = 0
!       n3dotn = 0
!       n4dotn = 0
!       n5dotn = 0
!       n6dotn = 0
!       n7dotn = 0
!       n8dotn = 0
!       n9dotn = 0
!     !  GOTO 1010
!      !!$acc parallel loop gang vector collapse(3)        &
!      !!$acc present (cosAlpha, cosBeta, cosGamma, xcent, ycent, zcent, cell, minElemcell, x1, y1, z1, xp, yp, zp)       &
!      !!$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, nel2n, m, cent_x, cent_y, cent_z,         &
!      !!$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn)
!     !$acc parallel loop gang vector collapse(3)        &
!     !$acc default(present)       &
!     !$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, nel2n, m, cent_x, cent_y, cent_z,         &
!     !$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn)
!
!       DO 10 k = block(g)%k_startSearch, block(g)%k_endSearch
!       DO 10 j = block(g)%j_startSearch, block(g)%j_endSearch
!       DO 10 i = block(g)%i_startSearch, block(g)%i_endSearch
!
!          m = 0
!          minDis = 1e14
!
!          n1x = block(g)%xp(i)
!          n1y = block(g)%yp(j)
!          n1z = block(g)%zp(k)
!
!          n2x = block(g)%x1(i)
!          n3x = block(g)%x1(i+1)
!
!          n2y = block(g)%y1(j)
!          n3y = block(g)%y1(j+1)
!
!          n2z = block(g)%z1(k)
!          n3z = block(g)%z1(k+1)
!
!          !ibElems_cnt=block(g)%ibElemCnt

!         !print*,i,j,k
!          !$acc loop seq
!          !DO m = 1, ibElems_cnt
!          DO m = 1, block(g)%ibElems
!             cent_x = block(g)%xcent(m)
!             cent_y = block(g)%ycent(m)
!             cent_z = block(g)%zcent(m)
!             dis  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2)
!             IF (dis.LT.minDis) THEN
!                minDis   = dis
!                nel2n    = m
!             ENDIF
!          ENDDO
!          block(g)%minElemcell(i,j,k) = nel2n
!
!            !n1dotn  = (n1x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                    !(n1y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                    !(n1z - block(g)%zcent(nel2n))*cosGamma(nel2n)

!            n2dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!            n3dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!            n4dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)
!
!            n5dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!           n6dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!            n7dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!            n8dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)
!
!            n9dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                     (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                     (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)
!
!          !n       = i-1  + nx*(j-2)  + block(g)%nx*ny*(k-2)
!
!          IF ( n2dotn.LE.-1e-16 .AND. n3dotn.LE.-1e-16 .AND. n4dotn.LE.-1e-16 .AND. n5dotn.LE.-1e-16  &
!              .AND. n6dotn.LE.-1e-16 .AND. n7dotn.LE.-1e-16 .AND. n8dotn.LE.-1e-16 .AND. n9dotn.LE.-1e-16) THEN
!       	    block(g)%cell(i,j,k) = 1
!          ELSEIF (n2dotn.GT.-1e-16 .AND. n3dotn.GT.-1e-16 .AND. n4dotn.GT.-1e-16 .AND. n5dotn.GT.-1e-16 &
!       	    .AND. n6dotn.GT.-1e-16 .AND. n7dotn.GT.-1e-16 .AND. n8dotn.GT.-1e-16 .AND. n9dotn.GT.-1e-16) THEN
!                 block(g)%cell(i,j,k) = 0
!          ELSE
!                 block(g)%cell(i,j,k) = 2
!
!          ENDIF
!
!10      CONTINUE
!
!        !$acc end parallel
!        !!$acc update host(cell)
!        WRITE(filename1,1) g
! 1      FORMAT('sphere_f.',i3.3,".dat")
!         OPEN(11,FILE=filename1,status='unknown')
!       DO k = 1, block(g)%nz+2
!       DO j = 1, block(g)%ny+2
!       DO i = 1, block(g)%nx+2
!       WRITE(11,*) block(g)%cell(i,j,k), block(g)%minElemcell(i,j,k)
!       END DO
!       END DO
!       END DO
!       CLOSE(11)





!
!        block(g)%ibCellCount = 0
!        block(g)%solidCellCount = 0
!        block(g)%fluidCellCount = 0
!
!       !!$acc parallel loop gang vector collapse(3) default(present) reduction(+: solidCellCount, fluidCellCount, ibCellCount)
!        DO 20 k = 2, block(g)%nz+1
!        DO 20 j = 2, block(g)%ny+1
!        DO 20 i = 2, block(g)%nx+1
!           !n       = i-1  + nx*(j-2)  + block(g)%nx*block(g)%block(g)%ny*(k-2)
!           IF (block(g)%cell(i,j,k).eq.1) THEN
!              block(g)%solidCellCount = block(g)%solidCellCount + 1
!           ELSEIF (block(g)%cell(i,j,k).eq.0) THEN
!              block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
!           ELSEIF (block(g)%cell(i,j,k).eq.2) THEN
!              block(g)%ibCellCount = block(g)%ibCellCount + 1
!           ENDIF
!20      CONTINUE
!       !!$acc end parallel
!
!     GOTO 1000
!       !1010 open(82,file='inter.dat',status='unknown')
!       !  1010  continue
!      open(82,file='inter.dat',form='formatted')
!     !write(82,*)'variables = "x", "y","z", "var"'
!       read(82,*) cLine
!       do k = 2, block(g)%nz+1
!      do j = 2, block(g)%ny+1
!      do i = 2, block(g)%nx+1
!       !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
!      if(block(g)%cell(i,j,k).eq.2)then
!     ! write(82,*) block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 2
!      read(82,*) block(g)%xp(i),block(g)%yp(j),block(g)%zp(k)
!      read(82,*) block(g)%cell(i,j,k)
!              block(g)%ibCellCount = block(g)%ibCellCount + 1
!      endif
!      end do
!      end do
!       enddo
!     close(82)
!
!     !open(83,file='fluid.dat',status='unknown')
!     open(83,file='fluid.dat',form='formatted')
!     !write(83,*)'variables = "x", "y","z","var"'
!       read(83,*)cLine
!      do k = 2, block(g)%nz+1
!      do j = 2, block(g)%ny+1
!      do i = 2, block(g)%nx+1
!       !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
!      if(block(g)%cell(i,j,k).eq.0)then
!     ! write(83,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 0
!      read(83,*) block(g)%xp(i),block(g)%yp(j),block(g)%zp(k)
!      read(83,*) block(g)%cell(i,j,k)
!              block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
!      endif
!      end do
!      end do
!       enddo
!     close(83)
!
!     !open(84,file='solid.dat',status='unknown')
!     open(84,file='solid.dat',form='formatted')
!    ! write(84,*)'variables = "x", "y","z","var"'
!       read(84,*)cLine
!      do k = 2, block(g)%nz+1
!      do j = 2, block(g)%ny+1
!      do i = 2, block(g)%nx+1
!       !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
!      if(block(g)%cell(i,j,k).eq.1)then
!     ! write(84,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 1
!      read(84,*) block(g)%xp(i),block(g)%yp(j),block(g)%zp(k)
!      read(84,*) block(g)%cell(i,j,k)
!              block(g)%solidCellCount = block(g)%solidCellCount + 1
!      endif
!      end do
!      end do
!       enddo
!     close(84)

!1000   CONTINUE

!        WRITE(filename1,2) g
!2       FORMAT('sphere_cellcount_f.',i3.3,".dat")
!        OPEN(12,FILE=filename1,FORM='formatted')
!       WRITE(12,*) block(g)%solidCellCount, block(g)%fluidCellCount, block(g)%ibCellCount
!       CLOSE(12)


!
!       print*, 'search done'
!       Print*, 'imms. cells=', block(g)%ibCellCount
!       Print*, 'fluid cells=', block(g)%fluidCellCount
!        Print*, 'solid cells=', block(g)%solidCellCount

!       END DO
!    END SUBROUTINE tagging
!
!***********************************************************************
        SUBROUTINE readTagging
        use global
        implicit none
        INTEGER (kind = 8) ::   i, j, k,g
        CHARACTER(len=150)  :: filename1

        DO g=blk_start,nblocks
         WRITE(filename1,1) g
  1       FORMAT('butter_f.',i3.3,".dat")
          OPEN(11,FILE=filename1,status='unknown')
        DO k = 1, block(g)%nz+2
        DO j = 1, block(g)%ny+2
        DO i = 1, block(g)%nx+2
        READ(11,*) block(g)%cell(i,j,k), block(g)%nodeIdTag(i,j,k)
        END DO
        END DO
        END DO
        CLOSE(11)


         WRITE(filename1,2) g
 2        FORMAT('butter_cellcount_f.',i3.3,".dat")
         OPEN(12,FILE=filename1,FORM='formatted')
        READ(12,*) block(g)%solidCellCount, block(g)%fluidCellCount, block(g)%ibCellCount
        CLOSE(12)
        END DO

        END SUBROUTINE

!***********************************************************************

!    SUBROUTINE selectiveRetagging
!       USE global
!       IMPLICIT NONE
!       INTEGER, PARAMETER :: rk = selected_real_kind(8)
!       INTEGER (kind = 8) ::  n, m, i, j, k, n1, n2, n3, n4, i1, j1, k1, nn, &
!                              nel2n, sumId,g, tag_flag, ibElems_cnt
!       INTEGER            :: iPt, iPt1, ibcnt, flcnt, sldcnt
!       REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis , minDis, &
!                             n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
!                             cent_x, cent_y, cent_z
!
!       DO g=blk_start,nblocks
!        if( block(g)%blk_mv_tag .eq.0)then
!    ALLOCATE(block(g)%minElemcell(block(g)%nx+2,block(g)%ny+2,block(g)%nz+2))
!
!
!       tag_flag=1
!    block(g)%minElemcell(:,:,:) = 0
!    n1dotn = 0
!       n2dotn = 0
!       n3dotn = 0
!       n4dotn = 0
!       n5dotn = 0
!       n6dotn = 0
!       n7dotn = 0
!       n8dotn = 0
!       n9dotn = 0

!       !!$acc parallel loop gang vector       &
!       !!$acc present(interceptedIndexPtr, xp, yp, zp, x1, y1, z1, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, cell, minElemcell)  &
!       !!$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, m, nel2n, cent_x, cent_y, cent_z,   &
!       !!$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn, i1, j1, k1, i, j, k)
!
!     !$acc parallel loop gang vector         &
!     !$acc default(present)       &
!     !$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, nel2n, m, cent_x, cent_y, cent_z,         &
!     !$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn)
!       DO nn = 1, block(g)%ibCellCount
!           i1 = block(g)%interceptedIndexPtr(nn, 1)
!           j1 = block(g)%interceptedIndexPtr(nn, 2)
!           k1 = block(g)%interceptedIndexPtr(nn, 3)
!
!          !$acc loop collapse(3) seq
!           DO 10 k = k1-2, k1+3
!           DO 10 j = j1-2, j1+3
!           DO 10 i = i1-2, i1+3
!              m = 0
!           minDis = 1e14
!   	 n1x = block(g)%xp(i)
!   	 n1y = block(g)%yp(j)
!   	 n1z = block(g)%zp(k)
!
!   	 n2x = block(g)%x1(i)
!   	 n3x = block(g)%x1(i+1)
!
!   	 n2y = block(g)%y1(j)
!   	 n3y = block(g)%y1(j+1)
!
!   	 n2z = block(g)%z1(k)
!   	 n3z = block(g)%z1(k+1)
!
!          !ibElems_cnt=block(g)%ibElemCnt
!   	 !$acc loop seq
!   	 !DO m = 1, ibElems_cnt
!   	 DO m = 1, block(g)%ibElems
!   	    cent_x = block(g)%xcent(m)
!   	    cent_y = block(g)%ycent(m)
!   	    cent_z = block(g)%zcent(m)
!   	    dis  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2)
!   	    IF (dis.LT.minDis) THEN
!   		 minDis   = dis
!   		 nel2n    = m
!              ENDIF
!           ENDDO
!
!           block(g)%minElemcell(i,j,k) = nel2n
!
!   	 !n1dotn  = (n1x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!              !          (n1y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!              !          (n1z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!   	 n2dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!   	 n3dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!   	 n4dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)
!
!   	 n5dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n2z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!   	 n6dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!   	 n7dotn  = (n2x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)

!   	 n8dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n2y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)
!
!   	 n9dotn  = (n3x - block(g)%xcent(nel2n))*block(g)%cosAlpha(nel2n) + &
!                        (n3y - block(g)%ycent(nel2n))*block(g)%cosBeta(nel2n)  + &
!                        (n3z - block(g)%zcent(nel2n))*block(g)%cosGamma(nel2n)
!
!   	 !n       = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
!
!           IF (n2dotn.LE.-1e-16 .AND. n3dotn.LE.-1e-16 .AND. n4dotn.LE.-1e-16 .AND. n5dotn.LE.-1e-16  &
!   	    .AND. n6dotn.LE.-1e-16 .AND. n7dotn.LE.-1e-16 .AND. n8dotn.LE.-1e-16 .AND. n9dotn.LE.-1e-16) THEN
!   		   block(g)%cell(i,j,k) = 1
!           ELSEIF (n2dotn.GT.-1e-16 .AND. n3dotn.GT.-1e-16 .AND. n4dotn.GT.-1e-16 .AND. n5dotn.GT.-1e-16 &
!   	    .AND. n6dotn.GT.-1e-16 .AND. n7dotn.GT.-1e-16 .AND. n8dotn.GT.-1e-16 .AND. n9dotn.GT.-1e-16) THEN
!   		   block(g)%cell(i,j,k) = 0
!           ELSE
!   		   block(g)%cell(i,j,k) = 2
!   	 ENDIF
!10         CONTINUE
!     ENDDO
!       !$acc end parallel

!       !!$acc update host(cell)

!     block(g)%ibCellCount = 0
!     block(g)%solidCellCount = 0
!     block(g)%fluidCellCount = 0
!     ibcnt=0.
!     flcnt=0.
!     sldcnt=0.
!       !!$acc parallel loop collapse(3) present(cell) reduction(+: solidCellCount, fluidCellCount, ibCellCount)
!       !$acc parallel loop collapse(3) default(present) reduction(+: sldcnt, flcnt, ibcnt)
!        DO 20 k = 2, block(g)%nz+1
!        DO 20 j = 2, block(g)%ny+1
!        DO 20 i = 2, block(g)%nx+1
!        !n       = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
!           IF (block(g)%cell(i,j,k).eq.1) THEN
!   	 !block(g)%solidCellCount = block(g)%solidCellCount + 1
!   	    sldcnt = sldcnt + 1
!           ELSEIF (block(g)%cell(i,j,k).eq.0) THEN
!            !  block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
!   	    flcnt = flcnt + 1
!           ELSEIF (block(g)%cell(i,j,k).eq.2) THEN
!             ! block(g)%ibCellCount = block(g)%ibCellCount + 1
!   	    ibcnt = ibcnt + 1
!           ENDIF
!20      CONTINUE
!       !$acc end parallel
!
!     block(g)%ibCellCount = ibcnt
!     block(g)%solidCellCount = sldcnt
!     block(g)%fluidCellCount = flcnt
!      GOTO 1001
!
!      IF (ita.eq.500) THEN
!    open(82,file='inter.dat',status='unknown')
!       write(82,*)'variables = "x", "y","z", "var"'
!     do k = 2, block(g)%nz+1
!        do j = 2, block(g)%ny+1
!        do i = 2, block(g)%nx+1
!     !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
!        if(block(g)%cell(i,j,k).eq.2)then
!        write(82,*) block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 2
!        endif
!        end do
!        end do
!     enddo
!       close(82)
!
!       open(83,file='fluid.dat',status='unknown')
!       write(83,*)'variables = "x", "y","z","var"'
!        do k = 2, block(g)%nz+1
!        do j = 2, block(g)%ny+1
!        do i = 2, block(g)%nx+1
!     !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
!        if(block(g)%cell(i,j,k).eq.0)then
!        write(83,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 0
!        endif
!        end do
!        end do
!     enddo
!       close(83)
!
!       open(84,file='solid.dat',status='unknown')
!       write(84,*)'variables = "x", "y","z","var"'
!        do k = 2, block(g)%nz+1
!        do j = 2, block(g)%ny+1
!        do i = 2, block(g)%nx+1
!     !n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*ny*(k-2)
!        if(block(g)%cell(i,j,k).eq.1)then
!        write(84,*)block(g)%xp(i),block(g)%yp(j), block(g)%zp(k), 1
!        endif
!        end do
!        end do
!     enddo
!       close(84)
!      ENDIF
!1001  CONTINUE
!      print*, 'selective retagging', block(g)%ibCellCount, block(g)%fluidCellCount, block(g)%solidCellCount
!       ENDIF
!       END DO
!    END SUBROUTINE selectiveRetagging
!***********************************************************************


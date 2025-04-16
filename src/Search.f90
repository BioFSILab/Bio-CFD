module biocfd_search
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64, int32
  ! allow(use-all) - TODO: Aim to fix this in the future
  use global
  use biocfd_fine_interp, only: fineUpdate_mv
  use biocfd_fine_interp_bound, only : fineUpdate_bd_mv
  implicit NONE

  private

  public :: findDistnode, shiftSurfaceNodesInitial, computeSurfaceNorm
  public :: tagging_th, tagging_th_move, block_move_check, cellcount_solid
  public :: cellcount_solid_coarse, cellcount_solid_coarse_mv, change_block_coords
  public :: change_block_interface, computenormdistance, computesurfacevariables, findtscells
  public :: fine_block_cell, selectiveretagging_th

  contains
        SUBROUTINE findDistnode
        REAL(dp)      ::  dist, dist1, dist2
        INTEGER(int64) ::  i, g


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
        end subroutine findDistnode

     SUBROUTINE shiftSurfaceNodesInitial
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER(int64) ::  i, g
        REAL(dp)      ::  xr1, yr1, zr1, angt
        REAL(dp)      :: bdy,bdfr

        phase_angle = phase_angle*pi/180_rk
        aoa1 = aoa*pi/180_rk
        aoa2 = -aoa1
        alpha_m = alpha_m*pi/180_rk
        theta_m = theta_m*pi/180_rk
        a0y = 0.  !a0y
        ang_theta = 0.  !2._rk*pi*freq
        alpha_t=(alpha_m*0.5_dp)*(1+cos(ang_theta*(totime+deltat)+phase_angle))
        theta_t       =  theta_m*cos(ang_theta*(totime+deltat))
        DO g=blk_start,nblocks
        ALLOCATE (block(g)%xnode1(block(g)%ibNodes), block(g)%ynode1(block(g)%ibNodes), &
                  block(g)%znode1(block(g)%ibNodes) )
        block(g)%xnode1 = block(g)%xnode
        block(g)%ynode1 = block(g)%ynode
        block(g)%znode1 = block(g)%znode
        xfact=0.05_dp
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
        block(g)%u_init = 0.
        block(g)%u_final = 0.
        block(g)%v_init = 0.
        block(g)%v_final = 0.
        block(g)%w_init = 0.
        block(g)%w_final = 0.
        block(g)%xmove = 0.
        block(g)%ymove = 0.
        block(g)%zmove = 0.

        block(g)%a0 = block(g)%a0*pi/180_rk
         angt  =  2._rk*pi*block(g)%bfreq
         block(g)%xpth1=block(g)%xshift-(ita*dxmin*xfact)
        block(g)%xpth2=block(g)%xshift-(ita*dxmin*xfact)
        block(g)%ypth1=(block(g)%yamp)*sin(angt*block(g)%xshift)
        block(g)%ypth2=(block(g)%yamp)*sin(angt*block(g)%xshift)
        block(g)%piv_x = block(g)%xshift
        block(g)%piv_y = block(g)%yshift
        block(g)%piv_z = block(g)%zshift

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
        DO i = 1, block(g)%ibnodes
        IF (block(g)%ibNodeId(i)==51) THEN
            xr1 =  block(g)%xnode(i)
            zr1 =  block(g)%znode(i)*cos(aoa1) + block(g)%ynode(i)*sin(aoa1) + piv_pt &
                   - piv_pt*cos(aoa1)
            yr1 = -block(g)%znode(i)*sin(aoa1) + block(g)%ynode(i)*cos(aoa1) + piv_pt*sin(aoa1)
         ELSEIF (block(g)%ibNodeId(i)==52) THEN
            xr1 =  block(g)%xnode(i)
            zr1 =  block(g)%znode(i)*cos(aoa2) + block(g)%ynode(i)*sin(aoa2)  + piv_pt &
                   - piv_pt*cos(aoa2)
            yr1 = -block(g)%znode(i)*sin(aoa2) + block(g)%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
           ELSE
               xr1 = block(g)% xnode(i)
               zr1 = block(g)% znode(i)  !*cos(aoa1) + ynode(i)*sin(aoa1)
               yr1 = block(g)% ynode(i)  !*sin(aoa1) + ynode(i)*cos(aoa1)
           ENDIF
           block(g)%xnode1(i) = xr1+ block(g)%xShift
               block(g)%ynode1(i) = yr1+ block(g)%yShift
           block(g)%znode1(i) = zr1+ block(g)%zshift
        ENDDO

        END DO

      END SUBROUTINE shiftSurfaceNodesInitial

      SUBROUTINE computeSurfaceVariables

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER(int64) ::  i, g
        REAL(dp)      ::  xr1, yr1, zr1
        REAL(dp)      :: angg, angt
        REAL(dp)      :: bdy,bdfr
        CHARACTER(len=150) :: filename1

        DO g=blk_start,nblocks
        angg=90
        aoa1       =  (block(g)%a0)*sin(2._rk*pi*freq*(totime+deltat) + phase_angle)
        aoa2       = -aoa1
        ang_theta  =  2._rk*pi*freq
        bdfr=block(g)%bfreq
        bdy=block(g)%yamp
        angt  =  2._rk*pi*bdfr

        block(g)%xpth2=block(g)%xshift-(ita*dxmin*xfact)
        block(g)%ypth2=(block(g)%yamp)*sin(angt*block(g)%xpth2)
        block(g)%xchg=block(g)%xpth2-block(g)%xpth1
        block(g)%ypth1=block(g)%ypth2
        block(g)%xpth1=block(g)%xpth2
        PRINT*, "angles =", aoa1*180._dp/pi, aoa2*180._dp/pi
        block(g)%thetaDot1 = ang_theta*block(g)%a0*cos(2._rk*pi*freq*(totime+deltat) + phase_angle)
        block(g)%thetaDDot1 = -ang_theta*ang_theta*block(g)%a0*sin(2._rk*pi*freq*(totime+deltat) &
                              + phase_angle)
        block(g)%thetaDot2  = -block(g)%thetaDot1
        block(g)%thetaDDot2 = -block(g)%thetaDDot1

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
        WRITE(filename1,1) block(g)%fineg,re, g
      1  FORMAT('d',I4.4,'_index.',F8.2,'.',i3.1,".dat")
        OPEN(UNIT = 17, FILE = filename1,Access='Append', STATUS = 'unknown')
        write(17,14) totime, block(g)%nxty_cent, block(g)%inity_cent, block(g)%ymove, &
                     block(g)%nxtx_cent, block(g)%xmove
     close(17)
     14      FORMAT(7F15.8)
        print*,'centn',block(g)%nxty_cent,'centi',block(g)%inity_cent,'mv',block(g)%ychg
      DO i = 1, block(g)%ibnodes
      IF (block(g)%ibNodeId(i)==51) THEN
             xr1 =  block(g)%xnode(i)
             zr1 =  block(g)%znode(i)*cos(aoa1) + block(g)%ynode(i)*sin(aoa1) + piv_pt &
                    - piv_pt*cos(aoa1)
             yr1 = -block(g)%znode(i)*sin(aoa1) + block(g)%ynode(i)*cos(aoa1) + piv_pt*sin(aoa1)
         ELSEIF (block(g)%ibNodeId(i)==52) THEN
             xr1 =  block(g)%xnode(i)
             zr1 =  block(g)%znode(i)*cos(aoa2) + block(g)%ynode(i)*sin(aoa2)  + piv_pt &
                    - piv_pt*cos(aoa2)
             yr1 = -block(g)%znode(i)*sin(aoa2) + block(g)%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
           ELSE
               xr1 = block(g)% xnode(i)
               zr1 = block(g)% znode(i)
               yr1 = block(g)% ynode(i)
           ENDIF

           block(g)% xnode1(i) = xr1 +block(g)%xpth2
           block(g)% ynode1(i) = yr1 +block(g)%yshift +block(g)%ymove
           block(g)% znode1(i) = zr1 +block(g)%zshift +block(g)%zmove
        ENDDO
           write(*,*) block(g)%xnode1(1),block(g)% ynode1(1), block(g)% znode1(1)
        ENDDO

      END SUBROUTINE computeSurfaceVariables

      SUBROUTINE computeSurfaceNorm

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER(int64) ::  n, g  !c1, c2, c3, c4
        REAL(dp)      :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL, binor
        REAL(dp)      :: var_xcent, var_ycent, var_zcent

        DO g=blk_start, nblocks

        ALLOCATE (block(g)%xcent(block(g)%ibElems), block(g)%ycent(block(g)%ibElems), &
                  block(g)%zcent(block(g)%ibElems), &
                  block(g)%cosAlpha(block(g)%ibElems), block(g)%cosBeta(block(g)%ibElems), &
                  block(g)%cosGamma(block(g)%ibElems), &
                  block(g)%alpha3(block(g)%ibElems), block(g)%beta3(block(g)%ibElems), &
                  block(g)%gamma3(block(g)%ibElems))

        !compute centroid and direction cosines
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

           block(g)%xcent(n) =  var_xcent                  !centroid x coordinate element
           block(g)%ycent(n) =  var_ycent                  !centroid y coordinate element
           block(g)%zcent(n) =  var_zcent                  !centroid z coordinate element

           block(g)%cosAlpha(n) = (p2y-p1y)*(p3z-p1z)-(p3y-p1y)*(p2z-p1z)
           block(g)%cosBeta(n)  = (p2z-p1z)*(p3x-p1x)-(p3z-p1z)*(p2x-p1x)
           block(g)%cosGamma(n) = (p2x-p1x)*(p3y-p1y)-(p3x-p1x)*(p2y-p1y)

           block(g)%alpha3(n) = block(g)%cosAlpha(n)
           block(g)%beta3(n)  = block(g)%cosBeta(n)
           block(g)%gamma3(n) = block(g)%cosGamma(n)

           lenEL = dsqrt(block(g)%cosAlpha(n)**2 + block(g)%cosBeta(n)**2 + block(g)%cosGamma(n)**2)   !length of element

                binor=inor

           block(g)%cosAlpha(n) = block(g)%cosAlpha(n)/lenEl*binor               !direction cosine unit normal along x
           block(g)%cosBeta(n)  = block(g)%cosBeta(n)/lenEl*binor                !direction cosine unit normal along y
           block(g)%cosGamma(n) = block(g)%cosGamma(n)/lenEl*binor               !direction cosine unit normal along z
        ENDDO
       !$acc end parallel

        print*, 'SurfaceNorm done, inor =', inor

        END DO


     END SUBROUTINE computeSurfaceNorm

     SUBROUTINE tagging_th

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER(int64) :: g, n, m, i, j, k,  nel2Cen, nel2Pnt, sumNodeId
        REAL(dp)      :: n1x, n1y, n1z, n2x, n2y, n2z, minDis1, minDis, &
                         n2dotn, cent_x, cent_y, cent_z, dis_cen, dis_pnt

        CHARACTER(LEN=120) :: filename1
        DO g=blk_start, nblocks
        block(g)%ibCellCount = 0
        block(g)%fluidCellCount = 0
        block(g)%solidCellCount = 0
        block(g)%cell = 0
        block(g)%cell2 = 0
        block(g)%nodeIdTag = 0
        n2dotn = 0

 !$acc parallel loop collapse(3) default(present)
        DO k = block(g)%k_startSearch, block(g)%k_endSearch
        DO j = block(g)%j_startSearch, block(g)%j_endSearch
        DO i = block(g)%i_startSearch, block(g)%i_endSearch
            minDis  = 1e14_dp
            minDis1 = 1e14_dp

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
            IF((block(g)%x1(i)<=block(g)%xcent(nel2Cen) .AND. &
                block(g)%x1(i+1)>=block(g)%xcent(nel2Cen)).AND. &
               (block(g)%y1(j)<=block(g)%ycent(nel2Cen) .AND. &
                block(g)%y1(j+1)>=block(g)%ycent(nel2Cen)).AND. &
               (block(g)%z1(k)<=block(g)%zcent(nel2Cen) .AND. &
                block(g)%z1(k+1)>=block(g)%zcent(nel2Cen))) THEN
               block(g)%cell(i,j,k) = 2

            ENDIF

                        n2dotn  = (n2x - block(g)%xcent(nel2Pnt))*block(g)%cosAlpha(nel2Pnt) + &
                           (n2y - block(g)%ycent(nel2Pnt))*block(g)%cosBeta(nel2Pnt)  + &
                           (n2z - block(g)%zcent(nel2Pnt))*block(g)%cosGamma(nel2Pnt)

            IF (n2dotn>=-1e-16_dp) THEN
               block(g)%nodeIdTag(i,j,k) = 0
            ELSE
               block(g)%nodeIdTag(i,j,k) = 1
            ENDIF
         END DO
         END DO
         END DO
!$acc end parallel

!$acc parallel loop collapse(3) default(present)
           DO k = block(g)%k_startSearch, block(g)%k_endSearch
           DO j = block(g)%j_startSearch, block(g)%j_endSearch
           DO i = block(g)%i_startSearch, block(g)%i_endSearch
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
           END DO
           END DO
           END DO
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
         DO k = 2, block(g)%nz+1
         DO j = 2, block(g)%ny+1
         DO i = 2, block(g)%nx+1
               IF (block(g)%cell(i,j,k)==1) THEN
                  block(g)%solidCellCount = block(g)%solidCellCount + 1
               ELSEIF (block(g)%cell(i,j,k)==0) THEN
                  block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
               ELSEIF (block(g)%cell(i,j,k)==2) THEN
                  block(g)%ibCellCount = block(g)%ibCellCount + 1
               ENDIF
         END DO
         END DO
         END DO
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

         WRITE(filename1,2) g
 2       FORMAT('butter_cellcount_f.',i3.3,".dat")
         OPEN(12,FILE=filename1,FORM='formatted')
        WRITE(12,*) block(g)%solidCellCount, block(g)%fluidCellCount, block(g)%ibCellCount
        CLOSE(12)
         print*, 'search done'
         Print*, 'imms. cells=', block(g)%ibCellCount
        Print*, 'fluid cells=',block(g)%fluidCellCount
        Print*, 'solid cells=', block(g)%solidCellCount
        END DO
     END SUBROUTINE tagging_th

     SUBROUTINE tagging_th_move

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER(int64) :: g, n, m, i, j, k, nel2Cen, nel2Pnt, sumNodeId

        INTEGER            :: a_blk_no, b_blk_no
        REAL(dp)      :: n1x, n1y, n1z, n2x,n2y,n2z, minDis1, minDis, &
                              n2dotn, cent_x, cent_y, cent_z, dis_cen, dis_pnt

        CHARACTER(LEN=120) :: filename1
        DO g=blk_start,nblocks
        if ( block(g)%move_check == 1)then
            block(g)% ibCellCount = 0
        block(g)%fluidCellCount = 0
        block(g)% solidCellCount = 0
        block(g)%cell = 0
        block(g)%cell2 = 0
        block(g)%nodeIdTag = 0
        n2dotn = 0

!$acc parallel loop collapse(3) default(present)
        DO k = block(g)%k_startSearch, block(g)%k_endSearch
        DO j = block(g)%j_startSearch, block(g)%j_endSearch
        DO i = block(g)%i_startSearch, block(g)%i_endSearch
            minDis  = 1e14_dp
            minDis1 = 1e14_dp

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
            IF((block(g)%x1(i)<=block(g)%xcent(nel2Cen).AND. &
                block(g)%x1(i+1)>=block(g)%xcent(nel2Cen)).AND. &
               (block(g)%y1(j)<=block(g)%ycent(nel2Cen).AND. &
                block(g)%y1(j+1)>=block(g)%ycent(nel2Cen)).AND. &
               (block(g)%z1(k)<=block(g)%zcent(nel2Cen).AND. &
                block(g)%z1(k+1)>=block(g)%zcent(nel2Cen))) THEN
               block(g)%cell(i,j,k) = 2

            ENDIF

                        n2dotn  = (n2x - block(g)%xcent(nel2Pnt))*block(g)%cosAlpha(nel2Pnt) + &
                           (n2y - block(g)%ycent(nel2Pnt))*block(g)%cosBeta(nel2Pnt)  + &
                           (n2z - block(g)%zcent(nel2Pnt))*block(g)%cosGamma(nel2Pnt)

            IF (n2dotn>=-1e-16_dp) THEN
               block(g)%nodeIdTag(i,j,k) = 0
            ELSE
               block(g)%nodeIdTag(i,j,k) = 1
            ENDIF
         END DO
         END DO
         END DO
!$acc end parallel

!$acc parallel loop collapse(3) default(present)
           DO k = block(g)%k_startSearch, block(g)%k_endSearch
           DO j = block(g)%j_startSearch, block(g)%j_endSearch
           DO i = block(g)%i_startSearch, block(g)%i_endSearch
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
           END DO
           END DO
           END DO
!$acc end parallel

        block(g)%ibCellCount = 0
         block(g)%solidCellCount = 0
         block(g)%fluidCellCount = 0
         DO k = 2, block(g)%nz+1
         DO j = 2, block(g)%ny+1
         DO i = 2, block(g)%nx+1
            IF (block(g)%cell(i,j,k)==1) THEN
                block(g)%solidCellCount = block(g)%solidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==0) THEN
               block(g)%fluidCellCount  = block(g)%fluidCellCount + 1
            ELSEIF (block(g)%cell(i,j,k)==2) THEN
                block(g)%ibCellCount = block(g)%ibCellCount + 1
            ENDIF
         END DO
         END DO
         END DO

         print*, 'search done'
         Print*, 'imms. cells=', block(g)%ibCellCount
        Print*, 'fluid cells=',block(g)%fluidCellCount
        Print*, 'solid cells=', block(g)%solidCellCount
        endif
        ENDDO
        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

        if ( block(b_blk_no)%move_check == 1)then
       call fineUpdate_mv(g)
       call fineUpdate_bd_mv(g)
        endif
        ENDDO
     END SUBROUTINE tagging_th_move

     SUBROUTINE findTScells

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER            :: g,i, j, k, i1, j1, k1, iPt1, m, n, tscnt
        CHARACTER(len=70) :: filename1


        print*,'inside findTScells'

        DO g=blk_start, nblocks
        !$acc parallel loop collapse(3) default(present)
                 DO k = 2, block(g)%nz+1
                 DO j = 2, block(g)%ny+1
                 DO i = 2, block(g)%nx+1
                    block(g)%cell2(i,j,k) = 0
                 END DO
                 END DO
                 END DO
        !$acc end parallel

        !$acc parallel loop default(present)
        DO n = 1, block(g)%ibCellCount
        i = block(g)%interceptedIndexPtr(n, 1)
        j = block(g)%interceptedIndexPtr(n, 2)
        k = block(g)%interceptedIndexPtr(n, 3)
         IF(block(g)%ibSurfID(block(g)%nelp(n))==51.OR.block(g)%ibSurfID(block(g)%nelp(n))==52) &
            block(g)%cell2(i, j, k) = 2

         END DO
        !$acc end parallel

        block(g)%TSCellCount = 0
        tscnt=0
        !$acc parallel loop collapse(3) default(present) reduction(+:tscnt)
                 DO k = 2, block(g)%nz+1
                 DO j = 2, block(g)%ny+1
                 DO i = 2, block(g)%nx+1
                     IF (block(g)%cell2(i,j,k)==2) THEN
                           !block(g)%TSCellCount = block(g)%TSCellCount + 1
                           tscnt = tscnt + 1
                     ENDIF
                 END DO
                 END DO
                 END DO
        !$acc end parallel

        block(g)%TSCellCount = tscnt
        print*, 'TScell count =', block(g)%TSCellCount

        ALLOCATE(block(g)%TSIndexPtr(block(g)%TSCellCount,3))

        iPt1 = 0
        !$acc loop collapse(3) seq
                 DO k = 0, block(g)%nz+3
                 DO j = 0, block(g)%ny+3
                 DO i = 0, block(g)%nx+3
                    IF (block(g)%cell2(i,j,k)==2) THEN
                       iPt1 = iPt1 + 1
                       block(g)%TSIndexPtr(iPt1, 1) = i
                           block(g)%TSIndexPtr(iPt1, 2) = j
                               block(g)%TSIndexPtr(iPt1, 3) = k
                    ENDIF
                 END DO
                 END DO
                 END DO

        ALLOCATE(&
          block(g)%u2_ghost(block(g)%TSCellCount),  &
          block(g)%u2t_ghost(block(g)%TSCellCount), &
          block(g)%v2_ghost(block(g)%TSCellCount), &
          block(g)%v2t_ghost(block(g)%TSCellCount), &
          block(g)%p_ghost(block(g)%TSCellCount), &
          block(g)%pt_ghost(block(g)%TSCellCount), &
          block(g)%u1_ghost(block(g)%TSCellCount),  &
          block(g)%u1t_ghost(block(g)%TSCellCount), &
          block(g)%v1_ghost(block(g)%TSCellCount), &
          block(g)%v1t_ghost(block(g)%TSCellCount), &
          block(g)%w2_ghost(block(g)%TSCellCount), &
          block(g)%w2t_ghost(block(g)%TSCellCount), &
          block(g)%w1_ghost(block(g)%TSCellCount),  &
          block(g)%w1t_ghost(block(g)%TSCellCount), &
          block(g)%index_ts(block(g)%TSCellCount))

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

           enddo
             END SUBROUTINE findTScells

     SUBROUTINE selectiveRetagging_th

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER(int64) ::  n, g, m, i, j, k, i1, j1, k1, nn, &
                               nel2Pnt, nel2Cen, sumNodeID
        INTEGER            :: flcnt, sdcnt, ibcnt
        REAL(dp)      :: n1x, n1y, n1z, n2x, n2y, n2z, minDis, minDis1, dis_cen, dis_pnt, &
                              n2dotn, cent_x, cent_y, cent_z

       DO g=blk_start,nblocks
        if( block(g)%blk_mv_tag ==0)then
!$acc parallel loop gang vector default(present)
        DO nn = 1, block(g)%ibCellCount
        i1 = block(g)%interceptedIndexPtr(nn, 1)
        j1 = block(g)%interceptedIndexPtr(nn, 2)
        k1 = block(g)%interceptedIndexPtr(nn, 3)
           block(g)%cell(i1,j1,k1) = 0
           !$acc loop collapse(3) seq
           DO k = k1-1, k1+1
           DO j = j1-1, j1+1
           DO i = i1-1, i1+1

               minDis  = 1e14_dp
               minDis1 = 1e14_dp

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
               IF((block(g)%x1(i)<=block(g)%xcent(nel2Cen).AND. &
                   block(g)%x1(i+1)>=block(g)%xcent(nel2Cen)).AND. &
                  (block(g)%y1(j)<=block(g)%ycent(nel2Cen).AND. &
                   block(g)%y1(j+1)>=block(g)%ycent(nel2Cen)).AND. &
                  (block(g)%z1(k)<=block(g)%zcent(nel2Cen).AND. &
                   block(g)%z1(k+1)>=block(g)%zcent(nel2Cen))) THEN
                        block(g)%cell(i,j,k) = 2
               ENDIF

                           n2dotn  = (n2x - block(g)%xcent(nel2Pnt))*block(g)%cosAlpha(nel2Pnt) + &
                              (n2y - block(g)%ycent(nel2Pnt))*block(g)%cosBeta(nel2Pnt)  + &
                              (n2z - block(g)%zcent(nel2Pnt))*block(g)%cosGamma(nel2Pnt)

               IF (n2dotn>=-1e-16_dp) THEN
                  block(g)%nodeIdTag(i,j,k) = 0
               ELSE
                  block(g)%nodeIdTag(i,j,k) = 1
               ENDIF
           END DO
           END DO
           END DO
        ENDDO
!$acc end parallel

!$acc parallel loop gang vector default(present)
        DO nn = 1, block(g)%ibCellCount
        i1 = block(g)%interceptedIndexPtr(nn, 1)
        j1 = block(g)%interceptedIndexPtr(nn, 2)
        k1 = block(g)%interceptedIndexPtr(nn, 3)
           !$acc loop collapse(3) seq
           DO k = k1-1, k1+1
           DO j = j1-1, j1+1
           DO i = i1-1, i1+1
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
           END DO
           END DO
           END DO
        ENDDO
  block(g)%ibCellCount = 0
  block(g)%solidCellCount = 0
block(g)%fluidCellCount = 0
sdcnt=0
flcnt=0
ibcnt=0
!$acc parallel loop gang vector collapse(3) default(present) private(i,j,k,n) reduction(+: sdcnt, flcnt, ibcnt)
         DO k = 2, block(g)%nz+1
         DO j = 2, block(g)%ny+1
         DO i = 2, block(g)%nx+1
            n = i-1  + block(g)%nx*(j-2)  + block(g)%nx*block(g)%ny*(k-2)
            IF (block(g)%cell(i,j,k)==1) THEN
                sdcnt = sdcnt + 1
            ELSEIF (block(g)%cell(i,j,k)==0) THEN
               flcnt  = flcnt + 1
            ELSEIF (block(g)%cell(i,j,k)==2) THEN
                ibcnt = ibcnt + 1
            ENDIF
         END DO
         END DO
         END DO
!$acc end parallel
  block(g)%ibCellCount = ibcnt
  block(g)%solidCellCount = sdcnt
block(g)%fluidCellCount = flcnt
       print*, 'selective retagging', block(g)%ibCellCount, block(g)%fluidCellCount, &
                block(g)%solidCellCount
     END IF
       ENDDO
     END SUBROUTINE selectiveRetagging_th

     SUBROUTINE cellCount_solid

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (int64) ::  n, iPt, iPt1, iPt2, i, j, k, g

        print*, "cellCount started"


         DO g=blk_start,nblocks
         iPt  = 0
         iPt1 = 0
         iPt2 = 0

         ALLOCATE(block(g)%interceptedIndexPtr(block(g)%ibCellCount,3), &
                  block(g)%fluidIndexPtr(block(g)%fluidCellCount, 3), &
                  block(g)%solidIndexPtr(block(g)%solidCellCount, 3))

         DO k = 2, block(g)%nz+1
         DO j = 2, block(g)%ny+1
         DO i = 2, block(g)%nx+1

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
         END DO
         END DO
         END DO
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
         ALLOCATE (block(g)%redCellIndexPtr(block(g)%redCellCount,3), &
                   block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))
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

            print*,g, block(g)%fluidCellCount, block(g)%redCellCount, block(g)%blackCellCount
        END DO
     END SUBROUTINE cellCount_solid

     SUBROUTINE computeNormDistance

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER            ::  nel2u1, nel2u2, nel2v1, nel2v2, nel2w1, nel2w2
        INTEGER            :: k, nel2p, g, ibxx, m
        REAL(dp) :: n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, &
                    dis, dis1, dis2, dis3, dis4, dis5, dis6, &
                    minDis, minDis1, minDis2, minDis3, minDis4, minDis5, minDis6, &
                    cent_x, cent_y,cent_z

        print*, 'computeNormDistance started'
        DO g=blk_start,nblocks

        ibxx=block(g)%ibCellCount
        print*,ibxx
        ALLOCATE(block(g)%pNormDis(ibxx), block(g)%nelp(ibxx), &
        block(g)%nelu1(ibxx), block(g)%nelu2(ibxx), block(g)%nelv1(ibxx), &
        block(g)%nelv2(ibxx), block(g)%nelw1(ibxx), block(g)%nelw2(ibxx), &
        block(g)%u1NormDis(ibxx), block(g)%u2NormDis(ibxx) , block(g)%v1NormDis(ibxx), &
        block(g)%v2NormDis(ibxx) , block(g)%w1NormDis(ibxx), block(g)%w2NormDis(ibxx))

        !$acc parallel loop gang vector default(present) &
        !$acc private (i, j, k, n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, nel2p)
        DO k = 1, block(g)%ibCellCount

           m = 0
           minDis  = 1e14_dp
           minDis1 = 1e14_dp
           minDis2 = 1e14_dp
           minDis3 = 1e14_dp
           minDis4 = 1e14_dp
           minDis5 = 1e14_dp
           minDis6 = 1e14_dp
           n1x = block(g)%xp(block(g)%interceptedIndexPtr(k, 1))
           n2x = block(g)%x1(block(g)%interceptedIndexPtr(k, 1))
           n3x = block(g)%x1(block(g)%interceptedIndexPtr(k, 1)+1)
           n1y = block(g)%yp(block(g)%interceptedIndexPtr(k, 2))
           n2y = block(g)%y1(block(g)%interceptedIndexPtr(k, 2))
           n3y = block(g)%y1(block(g)%interceptedIndexPtr(k, 2)+1)
           n1z = block(g)%zp(block(g)%interceptedIndexPtr(k, 3))
           n2z = block(g)%z1(block(g)%interceptedIndexPtr(k, 3))
           n3z = block(g)%z1(block(g)%interceptedIndexPtr(k, 3)+1)
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
        END DO
        !$acc end parallel

        ! DEALLOCATE (block(g)%minElemcell)
         END DO
         print*, 'computeNormDistance done'

     END SUBROUTINE computeNormDistance

        SUBROUTINE fine_block_cell

        INTEGER(int64) :: i, j, k, g, factor, a_blk_no, b_blk_no


        block(1)%cell_n=0
        block(1)%cell=0
        block(1)%cell_pr=0

        DO g=1, intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk
           factor=intfr(g)%b_msh/intfr(g)%a_msh

         DO k = 2, block(a_blk_no)%nz +1
         DO j = 2, block(a_blk_no)%ny +1
         DO i = 2, block(a_blk_no)%nx +1

      ! TN: Apologies for the horrible formatting, this is to please
      ! the linter, we will refactor this in the future anyway
        if(block(a_blk_no)%xp(i) >= intfr(g)%xintf_start+ block(b_blk_no)%cintp*block(a_blk_no)%dx &
     .and. block(a_blk_no)%xp(i) <= intfr(g)%xintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dx  &
     .and. block(a_blk_no)%zp(k) >= intfr(g)%zintf_start+ block(b_blk_no)%cintp*block(a_blk_no)%dz &
     .and. block(a_blk_no)%zp(k) <= intfr(g)%zintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dz  &
     .and. block(a_blk_no)%yp(j) >= intfr(g)%yintf_start + block(b_blk_no)%cintp*block(a_blk_no)%dy&
     .and. block(a_blk_no)%yp(j) <= intfr(g)%yintf_end - block(b_blk_no)%cintp*block(a_blk_no)%dy  &
          )then

                block(a_blk_no)%cell_n(i,j,k)=1
                block(a_blk_no)%cell(i,j,k)=1

        endif
        if ( block(a_blk_no)%xp(i) >= intfr(g)%xintf_start .and. &
             block(a_blk_no)%xp(i) <= intfr(g)%xintf_end   .and. &
             block(a_blk_no)%zp(k) >= intfr(g)%zintf_start .and. &
             block(a_blk_no)%zp(k) <= intfr(g)%zintf_end   .and. &
             block(a_blk_no)%yp(j) >= intfr(g)%yintf_start .and. &
             block(a_blk_no)%yp(j) <= intfr(g)%yintf_end   )then

                block(a_blk_no)%cell_pr(i,j,k)=1

        endif
        ENDDO
        ENDDO
        ENDDO
        ENDDO


        end subroutine fine_block_cell

        SUBROUTINE cellCount_solid_coarse_mv

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER(int64) ::  n, iPt, iPt1, iPt2, i, j, k
        INTEGER(int64) ::  g
        g=1

        if ( coarse_flcnt_check == 1)then

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
         DO k = 2, block(g)%nz +1
         DO j = 2, block(g)%ny +1
         DO i = 2, block(g)%nx +1
               n=    i-1  + (block(g)%nx)*(j-2)
               IF (block(g)%cell(i,j,k)==0) THEN
                  iPt1 = iPt1 + 1
                  block(g)%fluidIndexPtr(iPt1, 1) = i
                  block(g)%fluidIndexPtr(iPt1, 2) = j
                  block(g)%fluidIndexPtr(iPt1, 3) = k
               ENDIF
         END DO
         END DO
         END DO
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
         DEALLOCATE (block(g)%redCellIndexPtr ,block(g)%blackCellIndexPtr)
         ALLOCATE(block(g)%redCellIndexPtr(block(g)%redCellCount,3) , &
                  block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))
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

            print*, g, block(g)%fluidCellCount, block(g)%redCellCount, block(g)%blackCellCount
          ENDIF
          coarse_flcnt_check=0

        END SUBROUTINE cellCount_solid_coarse_mv

        SUBROUTINE block_move_check

        INTEGER(int64) :: i,j,k,g, a_blk_no, b_blk_no, factor
        REAL(dp) :: ydisp1, xdisp1, zdisp1, mg1
        REAL(dp) :: marginx, marginy, marginz, yval_up, yval_dw, xval_lt, xval_rt

        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        mg1=block(b_blk_no)%cintp*block(1)%dx
        marginx=4.5_dp*mg1
        marginy=5*mg1
        marginz=1.52_dp*mg1
        factor=intfr(g)%b_msh/intfr(g)%a_msh
        yval_up=dmax1(block(b_blk_no)%ynode1(block(b_blk_no)%mk),block(b_blk_no)%nxty_cent+0.02_dp)
        yval_dw=dmin1(block(b_blk_no)%ynode1(block(b_blk_no)%mk),block(b_blk_no)%nxty_cent-0.02_dp)
        xval_lt=(block(b_blk_no)%xnode1(block(b_blk_no)%mkx1))
        xval_rt=(block(b_blk_no)%xnode1(block(b_blk_no)%mkx2))
        ydisp1= dmin1(abs(yval_dw-intfr(g)%yintf_start),abs(yval_up - intfr(g)%yintf_end))
        xdisp1= (abs(xval_lt - intfr(g)%xintf_start))
        zdisp1= dmin1(abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start), &
                      abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end))
        print*,'blk_check_cond', marginx, marginy
        print*,'ydsip',abs(yval_dw-intfr(g)%yintf_start),abs(yval_up - intfr(g)%yintf_end)
        print*,block(b_blk_no)%xnode1(block(b_blk_no)%mkx1),block(b_blk_no)%nxtx_cent
        print*,block(b_blk_no)%ynode1(block(b_blk_no)%mk),block(b_blk_no)%nxty_cent
        block(b_blk_no)%move_amty=0
        block(b_blk_no)%move_amtx=0
        block(b_blk_no)%move_amtz=0
        if ( (abs(xval_lt - intfr(g)%xintf_start)      <= marginx )  .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start)      <= marginz )  .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end)        <= marginz )  .or. &
             ((abs(yval_dw - intfr(g)%yintf_start) <= marginy ) )   .or. &
             ((abs(yval_up - intfr(g)%yintf_end)  <= marginy ) ) )then

                block(b_blk_no)%move_check=1
             block(b_blk_no)%blk_mv_tag=1.
                 coarse_flcnt_check=1
            if( abs(yval_dw - intfr(g)%yintf_start) <= marginy    .or. &
                (abs(yval_up - intfr(g)%yintf_end)  <= marginy ) )then

                block(b_blk_no)%move_amty = &
                  floor((block(b_blk_no)%nxty_cent -block(b_blk_no)%inity_cent)/block(b_blk_no)%dy)
                if ( abs(block(b_blk_no)%move_amty) < factor)then
                        if ( block(b_blk_no)%move_amty < 0)then
                        block(b_blk_no)%move_amty=-factor
                        else
                        block(b_blk_no)%move_amty=factor
                        endif

                else
                        if ( block(b_blk_no)%move_amty < 0)then
                        block(b_blk_no)%move_amty= block(b_blk_no)%move_amty &
                                                   + mod(abs(block(b_blk_no)%move_amty),factor)
                        else
                        block(b_blk_no)%move_amty= block(b_blk_no)%move_amty &
                                                   - mod(abs(block(b_blk_no)%move_amty),factor)
                        endif

                endif
                block(b_blk_no)%inity_cent=block(b_blk_no)%nxty_cent

                endif
        if ( (abs(xval_lt- intfr(g)%xintf_start) < marginx ) )then
                block(b_blk_no)%move_amtx = &
                  floor((block(b_blk_no)%nxtx_cent - block(b_blk_no)%initx_cent) &
                  / block(b_blk_no)%dx) - factor
                print*,'move_amtx',block(b_blk_no)%move_amtx
                if ( abs(block(b_blk_no)%move_amtx) < factor)then
                        if ( block(b_blk_no)%move_amtx < 0)then
                        block(b_blk_no)%move_amtx=-factor
                        else
                        block(b_blk_no)%move_amtx=factor
                        endif

                else
                        if ( block(b_blk_no)%move_amtx < 0)then
                        block(b_blk_no)%move_amtx= block(b_blk_no)%move_amtx &
                                                   + mod(abs(block(b_blk_no)%move_amtx),factor)
                        else
                        block(b_blk_no)%move_amtx = block(b_blk_no)%move_amtx &
                                                    - mod(abs(block(b_blk_no)%move_amtx),factor)
                        endif

                endif
                block(b_blk_no)%initx_cent=block(b_blk_no)%nxtx_cent
                endif
            if( (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start) < (0.500_dp-0.410_dp-mg1)) .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end)   < (0.600_dp-0.500_dp-mg1) ))then
                block(b_blk_no)%move_amtz = &
                  floor((block(b_blk_no)%nxtz_cent -block(b_blk_no)%initz_cent)/block(b_blk_no)%dz)
                if ( abs(block(b_blk_no)%move_amtz) < factor)then
                        if ( block(b_blk_no)%move_amtz < 0)then
                        block(b_blk_no)%move_amtz=-factor
                        else
                        block(b_blk_no)%move_amtz=factor
                        endif

                else
                        if ( block(b_blk_no)%move_amtz < 0)then
                        block(b_blk_no)%move_amtz = block(b_blk_no)%move_amtz &
                                                    + mod(abs(block(b_blk_no)%move_amtz),factor)
                        else
                        block(b_blk_no)%move_amtz = block(b_blk_no)%move_amtz &
                                                    - mod(abs(block(b_blk_no)%move_amtz),factor)
                        endif

                endif
                block(b_blk_no)%initz_cent=block(b_blk_no)%nxtz_cent
                endif
                 print*,'blk_movez',block(b_blk_no)%move_amtz
                 print*,'blk_movey',block(b_blk_no)%move_amty
                 print*,'blk_movex',block(b_blk_no)%move_amtx

       intfr(g)%xintf_st_new=dmax1(intfr(g)%xintf_start,intfr(g)%xintf_start &
                             + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx))
       intfr(g)%xintf_en_new=dmin1(intfr(g)%xintf_end,intfr(g)%xintf_end &
                             + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx))
       intfr(g)%yintf_st_new=dmax1(intfr(g)%yintf_start,intfr(g)%yintf_start &
                             + (block(b_blk_no)%move_amty * block(b_blk_no)%dy))
       intfr(g)%yintf_en_new=dmin1(intfr(g)%yintf_end,intfr(g)%yintf_end &
                             + (block(b_blk_no)%move_amty * block(b_blk_no)%dy))
       intfr(g)%zintf_st_new=dmax1(intfr(g)%zintf_start,intfr(g)%zintf_start &
                             + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz))
       intfr(g)%zintf_en_new=dmin1(intfr(g)%zintf_end,intfr(g)%zintf_end &
                             + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz))
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
        endif

        ENDDO
        end subroutine block_move_check

        SUBROUTINE change_block_coords


        INTEGER(int64) :: i,j,k,g, a_blk_no, b_blk_no,countx_st,countz_st,county_st
        REAL(dp) :: change_y_f,change_x_f
        REAL(dp) :: change_z_f
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        DO g=blk_start,nblocks

        if ( block(g)% move_check == 1) then


       change_z_f= block(g)%move_amtz * block(g)%dz
       change_y_f= block(g)%move_amty * block(g)%dy
       change_x_f= block(g)%move_amtx * block(g)%dx
        print*, 'chz',block(g)%move_amtz , block(g)%dz
        print*, 'chy',block(g)%move_amty , block(g)%dy
        print*, 'chx',block(g)%move_amtx , block(g)%dx

        DO i = 1, block(g)%nx+3
           block(g)%x1(i) = block(g)%x1(i) + change_x_f
        ENDDO

        DO i = 1, block(g)%ny+3
                  ! print*,'y1p',block(g)%y1(i)
           block(g)%y1(i) = block(g)%y1(i) + change_y_f
                  ! print*,'y1a',block(g)%y1(i)
        ENDDO
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
        OPEN(UNIT=12,FILE='log.dat',STATUS='unknown',access='append')
        block(b_blk_no)%u=0
        block(b_blk_no)%v=0
        block(b_blk_no)%w=0
        block(b_blk_no)%p=0
        DO k=block(b_blk_no)%cpy_z_start_mv,block(b_blk_no)%cpy_z_end_mv
        countx_st=block(b_blk_no)%cpy_x_start
        DO i=block(b_blk_no)%cpy_x_start_mv,block(b_blk_no)%cpy_x_end_mv
        county_st=block(b_blk_no)%cpy_y_start
        DO j=block(b_blk_no)%cpy_y_start_mv,block(b_blk_no)%cpy_y_end_mv
                block(b_blk_no)%u(i,j,k)=block(b_blk_no)%u_dum(countx_st,county_st,countz_st)
                block(b_blk_no)%v(i,j,k)=block(b_blk_no)%v_dum(countx_st,county_st,countz_st)
                block(b_blk_no)%w(i,j,k)=block(b_blk_no)%w_dum(countx_st,county_st,countz_st)
                block(b_blk_no)%p(i,j,k)=block(b_blk_no)%p_dum(countx_st,county_st,countz_st)
                county_st=county_st+1
        ENDDO
                countx_st=countx_st+1
        ENDDO
                countz_st=countz_st+1
        ENDDO
        print*,'aft'
        close(12)
        ENDIF
        ENDDO



        end subroutine change_block_coords

        SUBROUTINE change_block_interface

        INTEGER(int64) :: i,j,g, a_blk_no, b_blk_no, factor
        INTEGER, PARAMETER :: rk = selected_real_kind(8)



        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
         factor=intfr(g)%b_msh/intfr(g)%a_msh
        if ( block(b_blk_no)% move_check == 1) then
        DO j=1,intfr(g)%counterxp

        intfr(g)%px_interface_det(1,j) = intfr(g)%px_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtx/factor)
        ENDDO
        DO j=1,intfr(g)%counterxu

        intfr(g)%ux_interface_det(1,j) = intfr(g)%ux_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtx/factor)

        ENDDO
        DO j=1,intfr(g)%counterxv

        intfr(g)%vx_interface_det(1,j) = intfr(g)%vx_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtx/factor)

        ENDDO
        DO j=1,intfr(g)%counterxw

        intfr(g)%wx_interface_det(1,j) = intfr(g)%wx_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtx/factor)

        ENDDO

        DO j=1,intfr(g)%counteryp

        intfr(g)%py_interface_det(1,j) = intfr(g)%py_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amty/factor)
        ENDDO
        DO j=1,intfr(g)%counteryu

        intfr(g)%uy_interface_det(1,j) = intfr(g)%uy_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counteryv

        intfr(g)%vy_interface_det(1,j) = intfr(g)%vy_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counteryw

        intfr(g)%wy_interface_det(1,j) = intfr(g)%wy_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counterzp

        intfr(g)%pz_interface_det(1,j) = intfr(g)%pz_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtz/factor)
        ENDDO
        DO j=1,intfr(g)%counterzu

        intfr(g)%uz_interface_det(1,j) = intfr(g)%uz_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtz/factor)

        ENDDO
        DO j=1,intfr(g)%counterzv

        intfr(g)%vz_interface_det(1,j) = intfr(g)%vz_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtz/factor)

        ENDDO
        DO j=1,intfr(g)%counterzw

        intfr(g)%wz_interface_det(1,j) = intfr(g)%wz_interface_det(1,j) &
                                         + (block(b_blk_no)%move_amtz/factor)

        ENDDO
        DO j=1,intflines
        print*,'**********************px***************************'
        DO i=1,intfr(j)%counterxp
        WRITE(*,33) 'px', i, &
                     intfr(j)%px_interface_det(1,i), &
                     intfr(j)%px_interface_det(2,i), &
                     intfr(j)%px_interface_det(3,i), &
                     block(a_blk_no)%xp(intfr(j)%px_interface_det(1,i)), &
                     block(b_blk_no)%xp(intfr(j)%px_interface_det(2,i)), &
                     block(b_blk_no)%xp(intfr(j)%px_interface_det(3,i))
 33       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************py***************************'
        DO i=1,intfr(j)%counteryp
        WRITE(*,133) 'py', i, &
                     intfr(j)%py_interface_det(1,i), &
                     intfr(j)%py_interface_det(2,i), &
                     intfr(j)%py_interface_det(3,i), &
                     block(a_blk_no)%yp(intfr(j)%py_interface_det(1,i)), &
                     block(b_blk_no)%yp(intfr(j)%py_interface_det(2,i)), &
                     block(b_blk_no)%yp(intfr(j)%py_interface_det(3,i))
 133       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intflines
        print*,'**********************pz***************************'
        DO i=1,intfr(j)%counterzp
        WRITE(*,933) 'pz', i, &
                     intfr(j)%pz_interface_det(1,i), &
                     intfr(j)%pz_interface_det(2,i), &
                     intfr(j)%pz_interface_det(3,i), &
                     block(a_blk_no)%zp(intfr(j)%pz_interface_det(1,i)), &
                     block(b_blk_no)%zp(intfr(j)%pz_interface_det(2,i)), &
                     block(b_blk_no)%zp(intfr(j)%pz_interface_det(3,i))
 933       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do

        DO j=1,intflines
        print*,'**********************ux***************************'
        DO i=1,intfr(j)%counterxu
        WRITE(*,331) 'ux', i, &
                     intfr(j)%ux_interface_det(1,i), &
                     intfr(j)%ux_interface_det(2,i), &
                     intfr(j)%ux_interface_det(3,i),  &
                     block(a_blk_no)%xu(intfr(j)%ux_interface_det(1,i)), &
                     block(b_blk_no)%xu(intfr(j)%ux_interface_det(2,i)), &
                     block(b_blk_no)%xu(intfr(j)%ux_interface_det(3,i))
 331       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************uy***************************'
        DO i=1,intfr(j)%counteryu
        WRITE(*,1331) 'uy', i, &
                     intfr(j)%uy_interface_det(1,i), &
                     intfr(j)%uy_interface_det(2,i), &
                     intfr(j)%uy_interface_det(3,i),  &
                     block(a_blk_no)%yu(intfr(j)%uy_interface_det(1,i)), &
                     block(b_blk_no)%yu(intfr(j)%uy_interface_det(2,i)), &
                     block(b_blk_no)%yu(intfr(j)%uy_interface_det(3,i))
 1331       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************uz***************************'
        DO i=1,intfr(j)%counterzu
        WRITE(*,91331) 'uz', i, &
                     intfr(j)%uz_interface_det(1,i), &
                     intfr(j)%uz_interface_det(2,i), &
                     intfr(j)%uz_interface_det(3,i),  &
                     block(a_blk_no)%zu(intfr(j)%uz_interface_det(1,i)), &
                     block(b_blk_no)%zu(intfr(j)%uz_interface_det(2,i)), &
                     block(b_blk_no)%zu(intfr(j)%uz_interface_det(3,i))
91331       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do

        DO j=1,intflines
        print*,'**********************vx***************************'
        DO i=1,intfr(j)%counterxv
        WRITE(*,332) 'vx', i, &
                     intfr(j)%vx_interface_det(1,i), &
                     intfr(j)%vx_interface_det(2,i), &
                     intfr(j)%vx_interface_det(3,i),  &
                     block(a_blk_no)%xv(intfr(j)%vx_interface_det(1,i)), &
                     block(b_blk_no)%xv(intfr(j)%vx_interface_det(2,i)), &
                     block(b_blk_no)%xv(intfr(j)%vx_interface_det(3,i))
 332       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************vy***************************'
        DO i=1,intfr(j)%counteryv
        WRITE(*,1332) 'vy', i, &
                     intfr(j)%vy_interface_det(1,i), &
                     intfr(j)%vy_interface_det(2,i), &
                     intfr(j)%vy_interface_det(3,i),  &
                     block(a_blk_no)%yv(intfr(j)%vy_interface_det(1,i)), &
                     block(b_blk_no)%yv(intfr(j)%vy_interface_det(2,i)), &
                     block(b_blk_no)%yv(intfr(j)%vy_interface_det(3,i))
 1332       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************vz***************************'
        DO i=1,intfr(j)%counterzv
        WRITE(*,91332) 'vz', i, &
                     intfr(j)%vz_interface_det(1,i), &
                     intfr(j)%vz_interface_det(2,i), &
                     intfr(j)%vz_interface_det(3,i),  &
                     block(a_blk_no)%zv(intfr(j)%vz_interface_det(1,i)), &
                     block(b_blk_no)%zv(intfr(j)%vz_interface_det(2,i)), &
                     block(b_blk_no)%zv(intfr(j)%vz_interface_det(3,i))
91332       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do


        DO j=1,intflines
        print*,'**********************wx***************************'
        DO i=1,intfr(j)%counterxw
        WRITE(*,3329) 'wx', i, &
                      intfr(j)%wx_interface_det(1,i), &
                      intfr(j)%wx_interface_det(2,i), &
                      intfr(j)%wx_interface_det(3,i), &
                      block(a_blk_no)%xw(intfr(j)%wx_interface_det(1,i)), &
                      block(b_blk_no)%xw(intfr(j)%wx_interface_det(2,i)), &
                      block(b_blk_no)%xw(intfr(j)%wx_interface_det(3,i))
 3329       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************wy***************************'
        DO i=1,intfr(j)%counteryw
        WRITE(*,13329) 'wy', i, &
                      intfr(j)%wy_interface_det(1,i), &
                      intfr(j)%wy_interface_det(2,i), &
                      intfr(j)%wy_interface_det(3,i), &
                      block(a_blk_no)%yw(intfr(j)%wy_interface_det(1,i)), &
                      block(b_blk_no)%yw(intfr(j)%wy_interface_det(2,i)), &
                      block(b_blk_no)%yw(intfr(j)%wy_interface_det(3,i))
13329       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************wz***************************'
        DO i=1,intfr(j)%counterzw
        WRITE(*,93329) 'wz', i, &
                      intfr(j)%wz_interface_det(1,i), &
                      intfr(j)%wz_interface_det(2,i), &
                      intfr(j)%wz_interface_det(3,i), &
                      block(a_blk_no)%zw(intfr(j)%wz_interface_det(1,i)), &
                      block(b_blk_no)%zw(intfr(j)%wz_interface_det(2,i)), &
                      block(b_blk_no)%zw(intfr(j)%wz_interface_det(3,i))
93329       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        endif
        ENDDO

        end subroutine change_block_interface

        SUBROUTINE cellCount_solid_coarse

        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (int64) ::  n, iPt, iPt1, iPt2, i, j, k
        INTEGER (int64) ::  g
        g=1

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
         DO k = 2, block(g)%nz +1
         DO j = 2, block(g)%ny +1
         DO i = 2, block(g)%nx +1
            IF (block(g)%cell(i,j,k)==0) THEN
               iPt1 = iPt1 + 1
               block(g)%fluidIndexPtr(iPt1, 1) = i
               block(g)%fluidIndexPtr(iPt1, 2) = j
               block(g)%fluidIndexPtr(iPt1, 3) = k
            ENDIF
         END DO
         END DO
         END DO
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

         ALLOCATE (block(g)%redCellIndexPtr(block(g)%redCellCount,3), &
                   block(g)%blackCellIndexPtr(block(g)%blackCellCount,3))
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

            print*, g, block(g)%fluidCellCount, block(g)%redCellCount, block(g)%blackCellCount

        END SUBROUTINE cellCount_solid_coarse
end module biocfd_search

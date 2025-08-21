module biocfd_search
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64, int32
  use global, only: block, nblocks, blk_start, xfact, totime, theta_t, &
       theta_m, piv_pt, pi, phase_angle, ita, dxmin, deltat, aoa2, aoa1, aoa, &
       ang_theta, alpha_t, alpha_m, ac_z, ac_y, ac_x, a0y, re, freq, inor, char_f, &
       intflines, coarse_flcnt_check, intfr
  use biocfd_fine_interp, only: fineUpdate_mv
  use biocfd_fine_interp_bound, only : fineUpdate_bd_mv
  use biocfd_blocks, only : Blocks
#ifdef BIOCFD_MPI
   use mpi_f08
#endif
  implicit NONE

  private

  public :: findDistnode, shiftSurfaceNodesInitial, computeSurfaceNorm
  public :: tagging_th, tagging_th_move, block_move_check, cellcount_solid
  public :: cellcount_solid_coarse, cellcount_solid_coarse_mv, change_block_coords
  public :: change_block_interface, computenormdistance, computesurfacevariables, findtscells
  public :: fine_block_cell, selectiveretagging_th

  contains
    SUBROUTINE findDistnode(blk)
        type(Blocks), intent(inout) :: blk
        REAL(dp)      ::  dist, dist1, dist2
        INTEGER(int64) ::  i

        dist=0.
        dist1=999999.
        dist2=0.
        DO i = 1, blk%ibnodes
        IF (blk%ibNodeId(i)==51) THEN
                if (abs(blk%znode(i)) > dist)then
                        dist=blk%znode(i)
                        blk%mk=i
                endif
        ENDIF
        IF (blk%ibNodeId(i)==51) THEN
                if (abs(blk%xnode(i)) < dist1)then
                        dist1=blk%xnode(i)
                        blk%mkx1=i
                endif
                if (abs(blk%xnode(i)) > dist2)then
                        dist2=blk%xnode(i)
                        blk%mkx2=i
                endif
        ENDIF
        ENDDO

        end subroutine findDistnode

     SUBROUTINE shiftSurfaceNodesInitial(blk)
        type(Blocks), intent(inout) :: blk
        INTEGER(int64) ::  i
        REAL(dp)      ::  xr1, yr1, zr1, angt
        REAL(dp)      :: bdy,bdfr

        phase_angle = phase_angle*pi/180_dp
        aoa1 = aoa*pi/180_dp
        aoa2 = -aoa1
        alpha_m = alpha_m*pi/180_dp
        theta_m = theta_m*pi/180_dp
        a0y = 0.  !a0y
        ang_theta = 0.  !2._dp*pi*freq
        alpha_t=(alpha_m*0.5_dp)*(1+cos(ang_theta*(totime+deltat)+phase_angle))
        theta_t       =  theta_m*cos(ang_theta*(totime+deltat))

        ALLOCATE (blk%xnode1(blk%ibNodes), blk%ynode1(blk%ibNodes), &
                  blk%znode1(blk%ibNodes) )
        blk%xnode1 = blk%xnode
        blk%ynode1 = blk%ynode
        blk%znode1 = blk%znode
        xfact=0.05_dp
        bdfr=15
        bdy=15*dxmin
        angt  =  2._dp*pi*bdfr

        ac_x=0.
        ac_y=0.
        ac_z=0.

        blk%xmove = 0.
        blk%ymove = 0.
        blk%zmove = 0.

        blk%a0 = blk%a0*pi/180_dp
         angt  =  2._dp*pi*blk%bfreq
         blk%xpth1=blk%xshift-(ita*dxmin*xfact)
        blk%xpth2=blk%xshift-(ita*dxmin*xfact)
        blk%ypth1=(blk%yamp)*sin(angt*blk%xshift)
        blk%ypth2=(blk%yamp)*sin(angt*blk%xshift)
        blk%piv_x = blk%xshift
        blk%piv_y = blk%yshift
        blk%piv_z = blk%zshift

        blk%thetaDot  =  0.
        blk% thetaDDot =  0._dp
        blk% thetaDot   =  0.  !ang_theta*a0*cos(2._dp*pi*freq*totime + phase_angle)
        blk%thetaDot1  = 0.
        blk% thetaDot2  = 0.
        blk% alphaDot  = 0.
        blk%alphaDDot  = 0.  !-ang_theta*ang_theta*a0*sin(2._dp*pi*freq*totime + phase_angle)
        blk% thetaDDot1 = 0.
        blk%thetaDDot2 = 0.
        blk% thetaDDot  = 0.  !-ang_theta*ang_theta*a0*sin(2._dp*pi*freq*totime + phase_angle)
       blk%yt         =  bdy*sin(2*pi*bdfr*totime )
       blk%ydot       =  angt*bdy*cos(2*pi*bdfr*totime)
       blk%yddot      =  -angt*angt*bdy*sin(2*pi*bdfr*totime)
       blk%xt         =  blk%xshift- (ita*dxmin*xfact)
       blk%xdot       =  -(dxmin*xfact)/deltat
       blk%xddot      =  0.
        blk%inity_cent=blk%yshift
        blk%nxty_cent=blk%yshift
        blk%initx_cent=blk%xshift
        blk%nxtx_cent=blk%xshift
        blk%initz_cent=blk%zshift
        blk%nxtz_cent=blk%zshift
        blk%ypos=blk%yshift
        blk%xpos=blk%xshift
        DO i = 1, blk%ibnodes
        IF (blk%ibNodeId(i)==51) THEN
            xr1 =  blk%xnode(i)
            zr1 =  blk%znode(i)*cos(aoa1) + blk%ynode(i)*sin(aoa1) + piv_pt &
                   - piv_pt*cos(aoa1)
            yr1 = -blk%znode(i)*sin(aoa1) + blk%ynode(i)*cos(aoa1) + piv_pt*sin(aoa1)
         ELSEIF (blk%ibNodeId(i)==52) THEN
            xr1 =  blk%xnode(i)
            zr1 =  blk%znode(i)*cos(aoa2) + blk%ynode(i)*sin(aoa2)  + piv_pt &
                   - piv_pt*cos(aoa2)
            yr1 = -blk%znode(i)*sin(aoa2) + blk%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
           ELSE
               xr1 = blk% xnode(i)
               zr1 = blk% znode(i)  !*cos(aoa1) + ynode(i)*sin(aoa1)
               yr1 = blk% ynode(i)  !*sin(aoa1) + ynode(i)*cos(aoa1)
           ENDIF
           blk%xnode1(i) = xr1+ blk%xShift
               blk%ynode1(i) = yr1+ blk%yShift
           blk%znode1(i) = zr1+ blk%zshift
        ENDDO

      END SUBROUTINE shiftSurfaceNodesInitial

      SUBROUTINE computeSurfaceVariables(blk, g)

        type(Blocks), intent(inout) :: blk
        integer(int64), intent(in) :: g
        INTEGER(int64) ::  i
        REAL(dp)      ::  xr1, yr1, zr1
        REAL(dp)      :: angg, angt
        REAL(dp)      :: bdy,bdfr
        CHARACTER(len=150) :: filename1

        angg=90
        aoa1       =  (blk%a0)*sin(2._dp*pi*freq*(totime+deltat) + phase_angle)
        aoa2       = -aoa1
        ang_theta  =  2._dp*pi*freq
        bdfr=blk%bfreq
        bdy=blk%yamp
        angt  =  2._dp*pi*bdfr

        blk%xpth2=blk%xshift-(ita*dxmin*xfact)
        blk%ypth2=(blk%yamp)*sin(angt*blk%xpth2)
        blk%xchg=blk%xpth2-blk%xpth1
        blk%ypth1=blk%ypth2
        blk%xpth1=blk%xpth2
        PRINT*, "angles =", aoa1*180._dp/pi, aoa2*180._dp/pi
        blk%thetaDot1 = ang_theta*blk%a0*cos(2._dp*pi*freq*(totime+deltat) + phase_angle)
        blk%thetaDDot1 = -ang_theta*ang_theta*blk%a0*sin(2._dp*pi*freq*(totime+deltat) &
                              + phase_angle)
        blk%thetaDot2  = -blk%thetaDot1
        blk%thetaDDot2 = -blk%thetaDDot1

       blk%yt         =  bdy*sin(angt*(totime) )
       blk%ydot       =  angt*bdy*cos(angt*(totime))
       blk%yddot      =  -angt*angt*bdy*sin(angt*(totime))

       blk%xt         = blk%xpth2
       blk%xdot       = -(dxmin*xfact)/deltat
       blk%xddot      = 0.

           blk%ychg=blk%yt - bdy*sin(angt*(totime-deltat) )
        blk%ymove = blk%yt
        blk%xmove = blk%xchg
        blk%zmove = 0.

        blk%ypos =  blk%ypos  + blk%ychg
        blk%xpos =  blk%xpos  + blk%xmove
        blk%piv_y = blk%piv_y + blk%ychg
        blk%piv_x = blk%piv_x + blk%xmove
        blk%piv_z = blk%piv_z + blk%zmove
        blk%nxty_cent= blk%nxty_cent + blk%ychg
        blk%nxtx_cent= blk%nxtx_cent + blk%xmove
        WRITE(filename1,1) blk%fineg,re, g
      1  FORMAT('d',I4.4,'_index.',F8.2,'.',i3.1,".dat")
        OPEN(UNIT = 17, FILE = filename1,POSITION='APPEND', STATUS = 'unknown')
        write(17,14) totime, blk%nxty_cent, blk%inity_cent, blk%ymove, &
                     blk%nxtx_cent, blk%xmove
     close(17)
     14      FORMAT(7F15.8)
        print*,'centn',blk%nxty_cent,'centi',blk%inity_cent,'mv',blk%ychg
      DO i = 1, blk%ibnodes
      IF (blk%ibNodeId(i)==51) THEN
             xr1 =  blk%xnode(i)
             zr1 =  blk%znode(i)*cos(aoa1) + blk%ynode(i)*sin(aoa1) + piv_pt &
                    - piv_pt*cos(aoa1)
             yr1 = -blk%znode(i)*sin(aoa1) + blk%ynode(i)*cos(aoa1) + piv_pt*sin(aoa1)
         ELSEIF (blk%ibNodeId(i)==52) THEN
             xr1 =  blk%xnode(i)
             zr1 =  blk%znode(i)*cos(aoa2) + blk%ynode(i)*sin(aoa2)  + piv_pt &
                    - piv_pt*cos(aoa2)
             yr1 = -blk%znode(i)*sin(aoa2) + blk%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
           ELSE
               xr1 = blk% xnode(i)
               zr1 = blk% znode(i)
               yr1 = blk% ynode(i)
           ENDIF

           blk% xnode1(i) = xr1 +blk%xpth2
           blk% ynode1(i) = yr1 +blk%yshift +blk%ymove
           blk% znode1(i) = zr1 +blk%zshift +blk%zmove
        ENDDO
           write(*,*) blk%xnode1(1),blk% ynode1(1), blk% znode1(1)

      END SUBROUTINE computeSurfaceVariables

      SUBROUTINE computeSurfaceNorm(blk)

         type(Blocks), intent(inout) :: blk
        INTEGER(int64) ::  n  !c1, c2, c3, c4
        REAL(dp)      :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL, binor
        REAL(dp)      :: var_xcent, var_ycent, var_zcent

        ALLOCATE (blk%xcent(blk%ibElems), blk%ycent(blk%ibElems), &
                  blk%zcent(blk%ibElems), &
                  blk%cosAlpha(blk%ibElems), blk%cosBeta(blk%ibElems), &
                  blk%cosGamma(blk%ibElems))

        !compute centroid and direction cosines
       !$acc parallel loop gang vector default(present) private (var_xcent, var_ycent, var_zcent,p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL)  firstprivate (inor)
        DO n = 1, blk%ibElems
           p1x = blk%xnode1(blk%ibElP1(n))                       !x coordinate element node 1
           p1y = blk%ynode1(blk%ibElP1(n))                       !y coordinate element node 1
           p1z = blk%znode1(blk%ibElP1(n))                       !z coordinate element node 1

           p2x = blk%xnode1(blk%ibElP2(n))                       !x coordinate element node 2
           p2y = blk%ynode1(blk%ibElP2(n))                       !y coordinate element node 2
           p2z = blk%znode1(blk%ibElP2(n))                       !z coordinate element node 2

           p3x = blk%xnode1(blk%ibElP3(n))                       !x coordinate element node 3
           p3y = blk%ynode1(blk%ibElP3(n))                       !y coordinate element node 3
           p3z = blk%znode1(blk%ibElP3(n))                       !z coordinate element node 3


           var_xcent =  (p2x+p1x+p3x)/3._dp                  !centroid x coordinate element
           var_ycent =  (p2y+p1y+p3y)/3._dp                  !centroid y coordinate element
           var_zcent =  (p2z+p1z+p3z)/3._dp                  !centroid z coordinate element

           blk%xcent(n) =  var_xcent                  !centroid x coordinate element
           blk%ycent(n) =  var_ycent                  !centroid y coordinate element
           blk%zcent(n) =  var_zcent                  !centroid z coordinate element

           blk%cosAlpha(n) = (p2y-p1y)*(p3z-p1z)-(p3y-p1y)*(p2z-p1z)
           blk%cosBeta(n)  = (p2z-p1z)*(p3x-p1x)-(p3z-p1z)*(p2x-p1x)
           blk%cosGamma(n) = (p2x-p1x)*(p3y-p1y)-(p3x-p1x)*(p2y-p1y)

           lenEL = dsqrt(blk%cosAlpha(n)**2 + blk%cosBeta(n)**2 + blk%cosGamma(n)**2)   !length of element

                binor=inor

           blk%cosAlpha(n) = blk%cosAlpha(n)/lenEl*binor               !direction cosine unit normal along x
           blk%cosBeta(n)  = blk%cosBeta(n)/lenEl*binor                !direction cosine unit normal along y
           blk%cosGamma(n) = blk%cosGamma(n)/lenEl*binor               !direction cosine unit normal along z
        ENDDO
       !$acc end parallel loop

        print*, 'SurfaceNorm done, inor =', inor

     END SUBROUTINE computeSurfaceNorm

     SUBROUTINE tagging_th(blk, g)

        type(Blocks), intent(inout) :: blk
        integer(int64), intent(in) :: g
        INTEGER(int64) :: n, m, i, j, k,  nel2Cen, nel2Pnt, sumNodeId
        REAL(dp)      :: minDis1, minDis, &
                         n2dotn, dis_cen, dis_pnt

        CHARACTER(LEN=120) :: filename1

        blk%ibCellCount = 0
        blk%fluidCellCount = 0
        blk%solidCellCount = 0
        blk%cell = 0
        blk%cell2 = 0
        blk%nodeIdTag = 0
        n2dotn = 0

 !$acc parallel loop collapse(3) default(present)
        ! There are some areas of the code where OpenMP and OpenACC
        ! are used together. In this particular case we don't want the
        ! OpenMP declarations if we are compiling with OpenACC.
        !
        ! WARNING: Do not indent the preprocessor macros, they must
        ! start at the first column.
        !
        ! WARNING 2: The extension of this file has been made
        ! uppercase "f90" -> "F90", this seems to be a convention
        ! across compilers, but I can't guarantee it will work for all
        ! compilers
#ifndef _OPENACC
        !$omp parallel do default(none) private(minDis, minDis1) &
        !$omp& private(dis_cen, dis_pnt, nel2Cen, nel2Pnt, n2dotn) &
        !$omp& shared(g, blk)
#endif
        DO k = blk%k_startSearch, blk%k_endSearch
        DO j = blk%j_startSearch, blk%j_endSearch
        DO i = blk%i_startSearch, blk%i_endSearch
            minDis  = 1e14_dp
            minDis1 = 1e14_dp

            !$acc loop seq
            DO m = 1, blk%ibElems
            ! I wanted to use associate here, but nvfortran doesn't
            ! like it on the GPU (although I can't find an existing
            ! bug report of this).
            !
            ! No need to take the sqrt because we are just looking for
            ! the minimum distance
               dis_cen  = (blk%xp(i)-blk%xcent(m))**2 &
                        + (blk%yp(j)-blk%ycent(m))**2 &
                        + (blk%zp(k)-blk%zcent(m))**2
               dis_pnt  = (blk%x1(i)-blk%xcent(m))**2 &
                        + (blk%y1(j)-blk%ycent(m))**2 &
                        + (blk%z1(k)-blk%zcent(m))**2
               IF (dis_cen<minDis) THEN
                  minDis    = dis_cen
                  nel2Cen   = m
               ENDIF
               IF (dis_pnt<minDis1) THEN
                  minDis1   = dis_pnt
                  nel2Pnt   = m
               ENDIF
            ENDDO
            IF((blk%x1(i)<=blk%xcent(nel2Cen) .AND. &
                blk%x1(i+1)>=blk%xcent(nel2Cen)).AND. &
               (blk%y1(j)<=blk%ycent(nel2Cen) .AND. &
                blk%y1(j+1)>=blk%ycent(nel2Cen)).AND. &
               (blk%z1(k)<=blk%zcent(nel2Cen) .AND. &
                blk%z1(k+1)>=blk%zcent(nel2Cen))) THEN
               blk%cell(i,j,k) = 2

            ENDIF

            n2dotn  = (blk%x1(i) - blk%xcent(nel2Pnt))*blk%cosAlpha(nel2Pnt) + &
                      (blk%y1(j) - blk%ycent(nel2Pnt))*blk%cosBeta(nel2Pnt)  + &
                      (blk%z1(k) - blk%zcent(nel2Pnt))*blk%cosGamma(nel2Pnt)

            IF (n2dotn>=-1e-16_dp) THEN
               blk%nodeIdTag(i,j,k) = 0
            ELSE
               blk%nodeIdTag(i,j,k) = 1
            ENDIF
         END DO
         END DO
         END DO
#ifndef _OPENACC
         !$omp end parallel do
#endif
!$acc end parallel loop

!$acc parallel loop collapse(3) default(present)
           DO k = blk%k_startSearch, blk%k_endSearch
           DO j = blk%j_startSearch, blk%j_endSearch
           DO i = blk%i_startSearch, blk%i_endSearch
               IF (blk%cell(i,j,k)/=2) THEN
                  sumNodeId = blk%nodeIdTag(i,j,k)      + blk%nodeIdTag(i+1,j,k)     &
                              + blk%nodeIdTag(i,j+1,k)    + blk%nodeIdTag(i+1,j+1,k)   &
                              + blk%nodeIdTag(i,j,k+1)    + blk%nodeIdTag(i+1,j,k+1)     &
                              + blk%nodeIdTag(i,j+1,k+1)  + blk%nodeIdTag(i+1,j+1,k+1)
                  IF (sumNodeId==8) THEN
                     blk%cell(i,j,k) = 1
                  ENDIF
               ENDIF
           END DO
           END DO
           END DO
!$acc end parallel loop


         WRITE(filename1,1) g
  1      FORMAT('butter_f.',i3.3,".dat")
          OPEN(11,FILE=filename1,status='unknown')
        DO k = 1, blk%nz+2
        DO j = 1, blk%ny+2
        DO i = 1, blk%nx+2
        WRITE(11,*) blk%cell(i,j,k), blk%nodeIdTag(i,j,k)
        END DO
        END DO
        END DO
        CLOSE(11)



         blk%ibCellCount = 0
         blk%solidCellCount = 0
         blk%fluidCellCount = 0
         DO k = 2, blk%nz+1
         DO j = 2, blk%ny+1
         DO i = 2, blk%nx+1
               IF (blk%cell(i,j,k)==1) THEN
                  blk%solidCellCount = blk%solidCellCount + 1
               ELSEIF (blk%cell(i,j,k)==0) THEN
                  blk%fluidCellCount  = blk%fluidCellCount + 1
               ELSEIF (blk%cell(i,j,k)==2) THEN
                  blk%ibCellCount = blk%ibCellCount + 1
               ENDIF
         END DO
         END DO
         END DO
      WRITE(filename1,22)char_f
 22   FORMAT(A3,'_inter_cell.dat')
      open(82,file=filename1,status='unknown')
      write(82,*)'variables = "x", "y","z", "var"'
    do k = 2,blk% nz+1
       do j = 2, blk%ny+1
       do i = 2, blk%nx+1
    n = i-1  + blk%nx*(j-2)  + blk%nx*blk%ny*(k-2)
       if(blk%cell(i,j,k)==2)then
       write(82,*) blk%xp(i),blk%yp(j), blk%zp(k), blk%cell(i,j,k)
       endif
       end do
       end do
    enddo
      close(82)
      WRITE(filename1,23)char_f
 23   FORMAT(A3,'_fluid_cell.dat')
      open(83,file=filename1,status='unknown')
      write(83,*)'variables = "x", "y","z","var"'
       do k = 2, blk%nz+1
       do j = 2,blk% ny+1
       do i = 2, blk%nx+1
    n = i-1  + blk%nx*(j-2)  + blk%nx*blk%ny*(k-2)
       if(blk%cell(i,j,k)==0)then
       write(83,*)blk%xp(i),blk%yp(j), blk%zp(k), 0
       endif
       end do
       end do
    enddo
      close(83)

      WRITE(filename1,24)char_f
 24   FORMAT(A3,'_solid_cell.dat')
      open(84,file=filename1,status='unknown')
      write(84,*)'variables = "x", "y","z","var"'
       do k = 2, blk%nz+1
       do j = 2, blk%ny+1
       do i = 2, blk%nx+1
    n = i-1  + blk%nx*(j-2)  + blk%nx*blk%ny*(k-2)
       if(blk%cell(i,j,k)==1)then
       write(84,*)blk%xp(i),blk%yp(j), blk%zp(k), 1
       endif
       end do
       end do
    enddo
      close(84)

         WRITE(filename1,2) g
 2       FORMAT('butter_cellcount_f.',i3.3,".dat")
         OPEN(12,FILE=filename1,FORM='formatted')
        WRITE(12,*) blk%solidCellCount, blk%fluidCellCount, blk%ibCellCount
        CLOSE(12)
         print*, 'search done'
         Print*, 'imms. cells=', blk%ibCellCount
        Print*, 'fluid cells=',blk%fluidCellCount
        Print*, 'solid cells=', blk%solidCellCount

     END SUBROUTINE tagging_th

     SUBROUTINE tagging_th_move

        INTEGER(int64) :: g, m, i, j, k, nel2Cen, nel2Pnt, sumNodeId
        INTEGER            :: a_blk_no, b_blk_no
        REAL(dp)      :: n1x, n1y, n1z, n2x,n2y,n2z, minDis1, minDis, &
                              n2dotn, cent_x, cent_y, cent_z, dis_cen, dis_pnt

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
!$acc end parallel loop

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
!$acc end parallel loop

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
#ifdef BIOCFD_MPI
          ! If we are using MPI then at this stage we need to make
          ! sure that block(1) is up-to-date on all ranks. We assume
          ! that all interfaces are from block(1) to another block
          call MPI_Bcast(block(1)%p, size(block(1)%p), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
          call MPI_Bcast(block(1)%u, size(block(1)%u), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
          call MPI_Bcast(block(1)%v, size(block(1)%v), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
          call MPI_Bcast(block(1)%w, size(block(1)%w), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
#endif
        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

        if ( block(b_blk_no)%move_check == 1)then
       call fineUpdate_mv(g)
       call fineUpdate_bd_mv(g)
        endif
        ENDDO
     END SUBROUTINE tagging_th_move

     SUBROUTINE findTScells(blk)


      type(Blocks), intent(inout) :: blk
        INTEGER            :: i, j, k, i1, j1, k1, iPt1, m, n, tscnt


        print*,'inside findTScells'

        !$acc parallel loop collapse(3) default(present)
                 DO k = 2, blk%nz+1
                 DO j = 2, blk%ny+1
                 DO i = 2, blk%nx+1
                    blk%cell2(i,j,k) = 0
                 END DO
                 END DO
                 END DO
        !$acc end parallel loop

        !$acc parallel loop default(present)
        DO n = 1, blk%ibCellCount
        i = blk%interceptedIndexPtr(n, 1)
        j = blk%interceptedIndexPtr(n, 2)
        k = blk%interceptedIndexPtr(n, 3)
         IF(blk%ibSurfID(blk%nelp(n))==51.OR.blk%ibSurfID(blk%nelp(n))==52) &
            blk%cell2(i, j, k) = 2

         END DO
        !$acc end parallel loop

        blk%TSCellCount = 0
        tscnt=0
        !$acc parallel loop collapse(3) default(present) reduction(+:tscnt)
                 DO k = 2, blk%nz+1
                 DO j = 2, blk%ny+1
                 DO i = 2, blk%nx+1
                     IF (blk%cell2(i,j,k)==2) THEN
                           !blk%TSCellCount = blk%TSCellCount + 1
                           tscnt = tscnt + 1
                     ENDIF
                 END DO
                 END DO
                 END DO
        !$acc end parallel loop

        blk%TSCellCount = tscnt
        print*, 'TScell count =', blk%TSCellCount

        ALLOCATE(blk%TSIndexPtr(blk%TSCellCount,3))

        iPt1 = 0
                 DO k = 0, blk%nz+3
                 DO j = 0, blk%ny+3
                 DO i = 0, blk%nx+3
                    IF (blk%cell2(i,j,k)==2) THEN
                       iPt1 = iPt1 + 1
                       blk%TSIndexPtr(iPt1, 1) = i
                           blk%TSIndexPtr(iPt1, 2) = j
                               blk%TSIndexPtr(iPt1, 3) = k
                    ENDIF
                 END DO
                 END DO
                 END DO

        ALLOCATE(&
          blk%u2_ghost(blk%TSCellCount),  &
          blk%u2t_ghost(blk%TSCellCount), &
          blk%v2_ghost(blk%TSCellCount), &
          blk%v2t_ghost(blk%TSCellCount), &
          blk%p_ghost(blk%TSCellCount), &
          blk%pt_ghost(blk%TSCellCount), &
          blk%u1_ghost(blk%TSCellCount),  &
          blk%u1t_ghost(blk%TSCellCount), &
          blk%v1_ghost(blk%TSCellCount), &
          blk%v1t_ghost(blk%TSCellCount), &
          blk%w2_ghost(blk%TSCellCount), &
          blk%w2t_ghost(blk%TSCellCount), &
          blk%w1_ghost(blk%TSCellCount),  &
          blk%w1t_ghost(blk%TSCellCount), &
          blk%index_ts(blk%TSCellCount))

        !$acc parallel loop default(present)
        DO n = 1, blk%TSCellCount
        i = blk%TSIndexPtr(n, 1)
        j = blk%TSIndexPtr(n, 2)
        k = blk%TSIndexPtr(n, 3)
           !$acc loop seq
           DO m = 1, blk%ibCellCount
           i1 = blk%interceptedIndexPtr(m, 1)
           j1 = blk%interceptedIndexPtr(m, 2)
           k1 = blk%interceptedIndexPtr(m, 3)
              IF (i1==i .AND. j1==j .AND. k1==k) THEN
                  blk%index_ts(n) = m
              ENDIF
           ENDDO
           blk%u2_ghost(n)  = 0.
           blk%u2t_ghost(n) = 0.
           blk%v2_ghost(n)  = 0.
           blk%v2t_ghost(n) = 0.
           blk%p_ghost(n)   = 0.
           blk%pt_ghost(n)  = 0.
           blk%u1_ghost(n)  = 0.
           blk%u1t_ghost(n) = 0.
           blk%v1_ghost(n)  = 0.
           blk%v1t_ghost(n) = 0.
           blk%w2_ghost(n)  = 0.
           blk%w2t_ghost(n) = 0.
           blk%w1_ghost(n)  = 0.
           blk%w1t_ghost(n) = 0.
        ENDDO
        !$acc end parallel loop

             END SUBROUTINE findTScells

     SUBROUTINE selectiveRetagging_th

        INTEGER(int64) ::  n, g, m, i, j, k, i1, j1, k1, nn, &
                               nel2Pnt, nel2Cen, sumNodeID
        INTEGER            :: flcnt, sdcnt, ibcnt
        REAL(dp)      :: minDis, minDis1, dis_cen, dis_pnt, n2dotn
        integer :: iprime, jprime, kprime

       DO g=blk_start,nblocks
        if( block(g)%blk_mv_tag ==0)then

       ! Set the intercepted indicies cell value to 0, we do this in a
       ! seperate loop so that we can nicely GPU-ise the computation
       !$acc parallel loop default(present) private(i1, j1, k1)
       DO nn = 1, block(g)%ibCellCount
         i1 = block(g)%interceptedIndexPtr(nn, 1)
         j1 = block(g)%interceptedIndexPtr(nn, 2)
         k1 = block(g)%interceptedIndexPtr(nn, 3)
         block(g)%cell(i1,j1,k1) = 0
        END DO
        !$acc end parallel loop

        !$acc parallel loop default(present) collapse(4) private(i, j, k) &
        !$acc private(minDis, minDis1, dis_cen, dis_pnt, nel2Cen, nel2Pnt) &
        !$acc private(n2dotn) firstprivate(g)
        DO nn = 1, block(g)%ibCellCount
           DO kprime = -1,+1
           DO jprime = -1,+1
           DO iprime = -1,+1

               i = block(g)%interceptedIndexPtr(nn, 1) + iprime
               j = block(g)%interceptedIndexPtr(nn, 2) + jprime
               k = block(g)%interceptedIndexPtr(nn, 3) + kprime

               minDis  = 1e14_dp
               minDis1 = 1e14_dp

               !$acc loop seq
               DO m = 1, block(g)%ibElems
                  ! No need to take sqrt because we just use for distance comparison
                  dis_cen  = (block(g)%xp(i)-block(g)%xcent(m))**2 &
                           + (block(g)%yp(j)-block(g)%ycent(m))**2 &
                           + (block(g)%zp(k)-block(g)%zcent(m))**2
                  dis_pnt  = (block(g)%x1(i)-block(g)%xcent(m))**2 &
                           + (block(g)%y1(j)-block(g)%ycent(m))**2 &
                           + (block(g)%z1(k)-block(g)%zcent(m))**2
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

               n2dotn  = (block(g)%x1(i) - block(g)%xcent(nel2Pnt)) * block(g)%cosAlpha(nel2Pnt) + &
                         (block(g)%y1(j) - block(g)%ycent(nel2Pnt)) * block(g)%cosBeta(nel2Pnt)  + &
                         (block(g)%z1(k) - block(g)%zcent(nel2Pnt)) * block(g)%cosGamma(nel2Pnt)

               IF (n2dotn>=-1e-16_dp) THEN
                  block(g)%nodeIdTag(i,j,k) = 0
               ELSE
                  block(g)%nodeIdTag(i,j,k) = 1
               ENDIF
           END DO
           END DO
           END DO
        ENDDO
!$acc end parallel loop

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
!$acc end parallel loop
  block(g)%ibCellCount = ibcnt
  block(g)%solidCellCount = sdcnt
block(g)%fluidCellCount = flcnt
       print*, 'selective retagging', block(g)%ibCellCount, block(g)%fluidCellCount, &
                block(g)%solidCellCount
     END IF
       ENDDO
     END SUBROUTINE selectiveRetagging_th

     SUBROUTINE cellCount_solid(blk, g)

        type(Blocks), intent(inout) :: blk
        integer(int64), intent(in) :: g
        INTEGER (int64) ::  n, iPt, iPt1, iPt2, i, j, k

        print*, "cellCount started"

         iPt  = 0
         iPt1 = 0
         iPt2 = 0

         ALLOCATE(blk%interceptedIndexPtr(blk%ibCellCount,3), &
                  blk%fluidIndexPtr(blk%fluidCellCount, 3), &
                  blk%solidIndexPtr(blk%solidCellCount, 3))

         DO k = 2, blk%nz+1
         DO j = 2, blk%ny+1
         DO i = 2, blk%nx+1

               IF (blk%cell(i,j,k)==0) THEN
                  iPt1 = iPt1 + 1
                  blk%fluidIndexPtr(iPt1, 1) = i
                  blk%fluidIndexPtr(iPt1, 2) = j
                  blk%fluidIndexPtr(iPt1, 3) = k
               ELSEIF (blk%cell(i,j,k)==1) THEN
                  iPt2 = iPt2 + 1
                  blk%solidIndexPtr(iPt2, 1) = i
                  blk%solidIndexPtr(iPt2, 2) = j
                  blk%solidIndexPtr(iPt2, 3) = k
               ELSEIF (blk%cell(i,j,k)==2) THEN
                  iPt = iPt + 1
                  blk%interceptedIndexPtr(iPt, 1) = i
                  blk%interceptedIndexPtr(iPt, 2) = j
                  blk%interceptedIndexPtr(iPt, 3) = k
               ENDIF
         END DO
         END DO
         END DO
         blk%redCellCount = 0
         blk%blackCellCount  = 0

         DO n = 1, blk%fluidCellCount
            i = blk%fluidIndexPtr(n, 1)
            j = blk%fluidIndexPtr(n, 2)
            k = blk%fluidIndexPtr(n, 3)

            IF (mod(i+j+k,2_int64)==1) THEN
               blk%redCellCount = blk%redCellCount + 1
            ELSE
               blk%blackCellCount = blk%blackCellCount + 1
            ENDIF
         ENDDO
         ALLOCATE (blk%redCellIndexPtr(blk%redCellCount,3), &
                   blk%blackCellIndexPtr(blk%blackCellCount,3))
         ipt1 = 0
         iPt = 0
         DO n = 1, blk%fluidCellCount
            i = blk%fluidIndexPtr(n, 1)
            j = blk%fluidIndexPtr(n, 2)
            k = blk%fluidIndexPtr(n, 3)

            IF (mod(i+j+k,2_int64)==1) THEN
               iPt = iPt + 1
               blk%redCellIndexPtr(iPt, 1) = i
               blk%redCellIndexPtr(iPt, 2) = j
               blk%redCellIndexPtr(iPt, 3) = k

            ELSE
               iPt1 = iPt1 + 1
               blk%blackCellIndexPtr(iPt1, 1) = i
               blk%blackCellIndexPtr(iPt1, 2) = j
               blk%blackCellIndexPtr(iPt1, 3) = k

            ENDIF
         ENDDO

         print*,g, blk%fluidCellCount, blk%redCellCount, blk%blackCellCount

     END SUBROUTINE cellCount_solid

     SUBROUTINE computeNormDistance(blk)


      type(Blocks), intent(inout) :: blk

        INTEGER            ::  nel2u1, nel2u2, nel2v1, nel2v2, nel2w1, nel2w2
        INTEGER            :: k, nel2p, ibxx, m
        REAL(dp) :: n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, &
                    dis, dis1, dis2, dis3, dis4, dis5, dis6, &
                    minDis, minDis1, minDis2, minDis3, minDis4, minDis5, minDis6, &
                    cent_x, cent_y,cent_z

        print*, 'computeNormDistance started'

        ibxx=blk%ibCellCount
        print*,ibxx
        ALLOCATE(blk%pNormDis(ibxx), blk%nelp(ibxx), &
        blk%nelu1(ibxx), blk%nelu2(ibxx), blk%nelv1(ibxx), &
        blk%nelv2(ibxx), blk%nelw1(ibxx), blk%nelw2(ibxx), &
        blk%u1NormDis(ibxx), blk%u2NormDis(ibxx) , blk%v1NormDis(ibxx), &
        blk%v2NormDis(ibxx) , blk%w1NormDis(ibxx), blk%w2NormDis(ibxx))

        !$acc parallel loop gang vector default(present) &
        !$acc private (k, n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, nel2p)
        DO k = 1, blk%ibCellCount

           m = 0
           minDis  = 1e14_dp
           minDis1 = 1e14_dp
           minDis2 = 1e14_dp
           minDis3 = 1e14_dp
           minDis4 = 1e14_dp
           minDis5 = 1e14_dp
           minDis6 = 1e14_dp
           n1x = blk%xp(blk%interceptedIndexPtr(k, 1))
           n2x = blk%x1(blk%interceptedIndexPtr(k, 1))
           n3x = blk%x1(blk%interceptedIndexPtr(k, 1)+1)
           n1y = blk%yp(blk%interceptedIndexPtr(k, 2))
           n2y = blk%y1(blk%interceptedIndexPtr(k, 2))
           n3y = blk%y1(blk%interceptedIndexPtr(k, 2)+1)
           n1z = blk%zp(blk%interceptedIndexPtr(k, 3))
           n2z = blk%z1(blk%interceptedIndexPtr(k, 3))
           n3z = blk%z1(blk%interceptedIndexPtr(k, 3)+1)
          !$acc loop seq
        DO m = 1, blk%ibElems
        cent_x = blk%xcent(m)
        cent_y = blk%ycent(m)
        cent_z = blk%zcent(m)
        ! No need to take sqrt because we just use for distance comparison
              dis   =  (n1x-cent_x)**2 + (n1y-cent_y)**2 + (n1z-cent_z)**2
              dis1  =  (n2x-cent_x)**2 + (n1y-cent_y)**2 + (n1z-cent_z)**2
              dis2  =  (n3x-cent_x)**2 + (n1y-cent_y)**2 + (n1z-cent_z)**2
              dis3  =  (n1x-cent_x)**2 + (n2y-cent_y)**2 + (n1z-cent_z)**2
              dis4  =  (n1x-cent_x)**2 + (n3y-cent_y)**2 + (n1z-cent_z)**2
              dis5  =  (n1x-cent_x)**2 + (n1y-cent_y)**2 + (n2z-cent_z)**2
              dis6  =  (n1x-cent_x)**2 + (n1y-cent_y)**2 + (n3z-cent_z)**2
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

           blk% nelp(k)  = nel2p
           blk% nelu1(k) = nel2u1
           blk% nelu2(k) = nel2u2
           blk% nelv1(k) = nel2v1
           blk% nelv2(k) = nel2v2
           blk% nelw1(k) = nel2w1
           blk% nelw2(k) = nel2w2
          blk% pNormDis(k)  = (n1x -blk% xcent(nel2p))*blk%cosAlpha(nel2p) + &
                   (n1y -blk% ycent(nel2p))*blk%cosBeta(nel2p)  + &
                   (n1z -blk% zcent(nel2p))*blk%cosGamma(nel2p)

          blk% u1NormDis(k) = (n2x -blk% xcent(nel2u1))*blk%cosAlpha(nel2u1) + &
                       (n1y -blk% ycent(nel2u1))*blk%cosBeta(nel2u1)  + &
                       (n1z - blk%zcent(nel2u1))*blk%cosGamma(nel2u1)

          blk% u2NormDis(k) = (n3x - blk%xcent(nel2u2))*blk%cosAlpha(nel2u2) + &
              (n1y -blk% ycent(nel2u2))*blk%cosBeta(nel2u2)  + &
              (n1z -blk% zcent(nel2u2))*blk%cosGamma(nel2u2)

          blk% v1NormDis(k) = (n1x - blk%xcent(nel2v1))*blk%cosAlpha(nel2v1) + &
              (n2y - blk%ycent(nel2v1))*blk%cosBeta(nel2v1)  + &
              (n1z - blk%zcent(nel2v1))*blk%cosGamma(nel2v1)

          blk% v2NormDis(k) = (n1x - blk%xcent(nel2v2))*blk%cosAlpha(nel2v2) + &
                   (n3y -blk% ycent(nel2v2))*blk%cosBeta(nel2v2)  + &
                   (n1z -blk% zcent(nel2v2))*blk%cosGamma(nel2v2)

          blk% w1NormDis(k) = (n1x -blk% xcent(nel2w1))*blk%cosAlpha(nel2w1) + &
                       (n1y -blk% ycent(nel2w1))*blk%cosBeta(nel2w1)  + &
                       (n2z -blk% zcent(nel2w1))*blk%cosGamma(nel2w1)

          blk% w2NormDis(k) = (n1x -blk% xcent(nel2w2))*blk%cosAlpha(nel2w2) + &
                           (n1y -blk% ycent(nel2w2))*blk%cosBeta(nel2w2)  + &
                           (n3z -blk% zcent(nel2w2))*blk%cosGamma(nel2w2)
        END DO
        !$acc end parallel loop

        ! DEALLOCATE (blk%minElemcell)
         print*, 'computeNormDistance done'

     END SUBROUTINE computeNormDistance

        SUBROUTINE fine_block_cell

        INTEGER(int64) :: i, j, k, g, factor, a_blk_no, b_blk_no


        block(1)%cell_n=0
        block(1)%cell=0
        block(1)%cell_pr=0

         ! TODO (MPI): This is an occasion where it is nice to have
         ! the interfaces available on all ranks as opposed to just
         ! having them attached to blocks - this probably leans
         ! towards having each block carrying its ID number?
        DO g=1, intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk
           factor=intfr(g)%b_msh/intfr(g)%a_msh

         DO k = 2, block(a_blk_no)%nz +1
         DO j = 2, block(a_blk_no)%ny +1
         DO i = 2, block(a_blk_no)%nx +1

      ! TN: Apologies for the horrible formatting, this is to please
      ! the linter, we will refactor this in the future anyway

      ! TODO: (MPI) block(b_blk_no)%cintp will be available on all ranks I think
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
            IF (mod(i+j+k,2_int64)==1) THEN
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
            IF (mod(i+j+k,2_int64)==1) THEN
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
#ifdef BIOCFD_MPI
         ! Here we use the trick we've used elsewhere that if one of
         ! the allocatable arrays isn't allocated we aren't on the
         ! right branch so carry on
         if (.not. allocated(block(b_blk_no)%u)) then
               ! If this array isn't allocated we aren't on the right
               ! rank to deal with this so keep going until we find
               ! one that is on this rank
               cycle
           end if
#endif
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
        endif

        ENDDO

#ifdef BIOCFD_MPI
        ! Here we set coarse_flcnt_check to 1 on every rank if it is 1
        ! on any rank (by computing the maximum). Since this is really
        ! a true or false flag, we should probably convert this to a
        ! logical. Then we'd use MPI_LOR instead of MPI_MAX for the
        ! reduction operation.
        call MPI_Allreduce(MPI_IN_PLACE, coarse_flcnt_check, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD)
#endif
        end subroutine block_move_check

        SUBROUTINE change_block_coords(blk, id)

        type(Blocks), intent(inout) :: blk
        !> TODO: The id is a bit of a hack to deal with interfaces in
        !> MPI. I think it would make a lot more sense for the
        !> interface to be a member of the Block type AND/OR id be a
        !> member of the Block type
        integer(int64), intent(in) :: id
        INTEGER(int64) :: i,j,k, g, countx_st,countz_st,county_st
        REAL(dp) :: change_y_f,change_x_f
        REAL(dp) :: change_z_f
       

        ! If move_check isn't 1 then do nothing
        if ( blk%move_check /= 1) return

        change_z_f= blk%move_amtz * blk%dz
        change_y_f= blk%move_amty * blk%dy
        change_x_f= blk%move_amtx * blk%dx
        print*, 'chz',blk%move_amtz , blk%dz
        print*, 'chy',blk%move_amty , blk%dy
        print*, 'chx',blk%move_amtx , blk%dx

        DO i = 1, blk%nx+3
           blk%x1(i) = blk%x1(i) + change_x_f
        ENDDO

        DO i = 1, blk%ny+3
                  ! print*,'y1p',blk%y1(i)
           blk%y1(i) = blk%y1(i) + change_y_f
                  ! print*,'y1a',blk%y1(i)
        ENDDO
        DO i = 1, blk%nz+3
           blk%z1(i) = blk%z1(i)+ change_z_f
        ENDDO

        DO i = 1, blk%nx+3
           blk%xu(i) = blk%x1(i)
        ENDDO

        DO i = 1, blk%ny+3
           blk%yv(i) = blk%y1(i)
        ENDDO

        DO i = 1, blk%nz+3
           blk%zw(i) = blk%z1(i)
        ENDDO

        DO i = 1, blk%ny+2
           blk%yu(i) = 0.5_dp*(blk%y1(i)+blk%y1(i+1))
           blk%yw(i) = blk%yu(i)
           blk%yp(i) = blk%yu(i)
        END DO

        DO i = 1, blk%nx+2
           blk%xv(i) = 0.5_dp*(blk%x1(i)+blk%x1(i+1))
           blk%xw(i) = blk%xv(i)
           blk%xp(i) = blk%xv(i)
        END DO

        DO i = 1, blk%nz+2
           blk%zu(i) = 0.5_dp*(blk%z1(i)+blk%z1(i+1))
           blk%zv(i) = blk%zu(i)
           blk%zp(i) = blk%zu(i)
        END DO

        DO g=1,intflines

         ! If b_blk isn't our id then try the next interface
          if (intfr(g)%b_blk /= id) cycle
              
        blk%cell_n=0
        DO k=1,blk%nz+2
        DO j=1,blk%ny+2
        DO i=1,blk%nx+2
        if( blk%xp(i) > (intfr(g)%xintf_st_new  + blk%dx) .and. &
               blk%xp(i) < (intfr(g)%xintf_en_new  - blk%dx).and. &
               blk%zp(k) > (intfr(g)%zintf_st_new  + blk%dz) .and. &
               blk%zp(k) < (intfr(g)%zintf_en_new  - blk%dz).and. &
               blk%yp(j) > (intfr(g)%yintf_st_new  + blk%dy).and. &
               blk%yp(j) < (intfr(g)%yintf_en_new  - blk%dy))then

                blk%cell_n(i,j,k)=1

           endif
        ENDDO
        ENDDO
        ENDDO

        DO i=1,blk%nx+2
           if( blk%xp(i) > intfr(g)%xintf_st_new)then
                blk%cpy_x_start_mv=i
                print*,'xpst_f',blk%xp(i),'i',i
                exit
          endif
        ENDDO

        DO i=blk%cpy_x_start_mv,blk%nx+2
           if( blk%xp(i) > intfr(g)%xintf_en_new)then
                blk%cpy_x_end_mv=i-1
                print*,'xpen_f',blk%xp(i-1),'i',i-1
                exit
          endif
        ENDDO
        DO j=1,blk%ny+2
           if( blk%yp(j) > intfr(g)%yintf_st_new)then
                blk%cpy_y_start_mv=j
                print*,'ypst_f',blk%yp(j),'j',j
                exit
          endif
        ENDDO

        DO j=blk%cpy_y_start_mv,blk%ny+2
           if( blk%yp(j) > intfr(g)%yintf_en_new)then
                blk%cpy_y_end_mv=j-1
                print*,'ypen_f',blk%yp(j-1),'j',j-1
                exit
          endif
        ENDDO
        DO j=1,blk%nz+2
           if( blk%zp(j) > intfr(g)%zintf_st_new)then
                blk%cpy_z_start_mv=j
                print*,'zpst_f',blk%zp(j),'k',j
                exit
          endif
        ENDDO

        DO j=blk%cpy_z_start_mv,blk%nz+2
           if( blk%zp(j) > intfr(g)%zintf_en_new)then
                blk%cpy_z_end_mv=j-1
                print*,'zpen_f',blk%zp(j-1),'k',j-1
                exit
          endif
        ENDDO
        print*,'bef'
        countx_st=blk%cpy_x_start
        county_st=blk%cpy_y_start
        countz_st=blk%cpy_z_start
        OPEN(UNIT=12,FILE='log.dat',STATUS='unknown',POSITION='APPEND')
        blk%u=0
        blk%v=0
        blk%w=0
        blk%p=0
        DO k=blk%cpy_z_start_mv,blk%cpy_z_end_mv
        countx_st=blk%cpy_x_start
        DO i=blk%cpy_x_start_mv,blk%cpy_x_end_mv
        county_st=blk%cpy_y_start
        DO j=blk%cpy_y_start_mv,blk%cpy_y_end_mv
                blk%u(i,j,k)=blk%u_dum(countx_st,county_st,countz_st)
                blk%v(i,j,k)=blk%v_dum(countx_st,county_st,countz_st)
                blk%w(i,j,k)=blk%w_dum(countx_st,county_st,countz_st)
                blk%p(i,j,k)=blk%p_dum(countx_st,county_st,countz_st)
                county_st=county_st+1
        ENDDO
                countx_st=countx_st+1
        ENDDO
                countz_st=countz_st+1
        ENDDO
        print*,'aft'
        close(12)
        ENDDO



        end subroutine change_block_coords

        SUBROUTINE change_block_interface(blk, id)

        type(Blocks), intent(inout) :: blk
        !> TODO: The id is a bit of a hack to deal with interfaces in
        !> MPI. I think it would make a lot more sense for the
        !> interface to be a member of the Block type AND/OR id be a
        !> member of the Block type
        integer(int64), intent(in) :: id
        INTEGER(int64) :: i,j,g, a_blk_no, b_blk_no, factor

        ! If blk%move_check isn't 1 then do nothing
        if ( blk% move_check /= 1) return

        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk

        if (intfr(g)%b_blk /= id) cycle

         factor=intfr(g)%b_msh/intfr(g)%a_msh
        DO j=1,intfr(g)%counterxp

        intfr(g)%px_interface_det(1,j) = intfr(g)%px_interface_det(1,j) &
                                         + (blk%move_amtx/factor)
        ENDDO
        DO j=1,intfr(g)%counterxu

        intfr(g)%ux_interface_det(1,j) = intfr(g)%ux_interface_det(1,j) &
                                         + (blk%move_amtx/factor)

        ENDDO
        DO j=1,intfr(g)%counterxv

        intfr(g)%vx_interface_det(1,j) = intfr(g)%vx_interface_det(1,j) &
                                         + (blk%move_amtx/factor)

        ENDDO
        DO j=1,intfr(g)%counterxw

        intfr(g)%wx_interface_det(1,j) = intfr(g)%wx_interface_det(1,j) &
                                         + (blk%move_amtx/factor)

        ENDDO

        DO j=1,intfr(g)%counteryp

        intfr(g)%py_interface_det(1,j) = intfr(g)%py_interface_det(1,j) &
                                         + (blk%move_amty/factor)
        ENDDO
        DO j=1,intfr(g)%counteryu

        intfr(g)%uy_interface_det(1,j) = intfr(g)%uy_interface_det(1,j) &
                                         + (blk%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counteryv

        intfr(g)%vy_interface_det(1,j) = intfr(g)%vy_interface_det(1,j) &
                                         + (blk%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counteryw

        intfr(g)%wy_interface_det(1,j) = intfr(g)%wy_interface_det(1,j) &
                                         + (blk%move_amty/factor)

        ENDDO
        DO j=1,intfr(g)%counterzp

        intfr(g)%pz_interface_det(1,j) = intfr(g)%pz_interface_det(1,j) &
                                         + (blk%move_amtz/factor)
        ENDDO
        DO j=1,intfr(g)%counterzu

        intfr(g)%uz_interface_det(1,j) = intfr(g)%uz_interface_det(1,j) &
                                         + (blk%move_amtz/factor)

        ENDDO
        DO j=1,intfr(g)%counterzv

        intfr(g)%vz_interface_det(1,j) = intfr(g)%vz_interface_det(1,j) &
                                         + (blk%move_amtz/factor)

        ENDDO
        DO j=1,intfr(g)%counterzw

        intfr(g)%wz_interface_det(1,j) = intfr(g)%wz_interface_det(1,j) &
                                         + (blk%move_amtz/factor)

        ENDDO
        DO j=1,intflines
        print*,'**********************px***************************'
        DO i=1,intfr(j)%counterxp
        WRITE(*,33) 'px', i, &
                     intfr(j)%px_interface_det(1,i), &
                     intfr(j)%px_interface_det(2,i), &
                     intfr(j)%px_interface_det(3,i), &
                     block(a_blk_no)%xp(intfr(j)%px_interface_det(1,i)), &
                     blk%xp(intfr(j)%px_interface_det(2,i)), &
                     blk%xp(intfr(j)%px_interface_det(3,i))
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
                     blk%yp(intfr(j)%py_interface_det(2,i)), &
                     blk%yp(intfr(j)%py_interface_det(3,i))
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
                     blk%zp(intfr(j)%pz_interface_det(2,i)), &
                     blk%zp(intfr(j)%pz_interface_det(3,i))
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
                     blk%xu(intfr(j)%ux_interface_det(2,i)), &
                     blk%xu(intfr(j)%ux_interface_det(3,i))
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
                     blk%yu(intfr(j)%uy_interface_det(2,i)), &
                     blk%yu(intfr(j)%uy_interface_det(3,i))
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
                     blk%zu(intfr(j)%uz_interface_det(2,i)), &
                     blk%zu(intfr(j)%uz_interface_det(3,i))
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
                     blk%xv(intfr(j)%vx_interface_det(2,i)), &
                     blk%xv(intfr(j)%vx_interface_det(3,i))
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
                     blk%yv(intfr(j)%vy_interface_det(2,i)), &
                     blk%yv(intfr(j)%vy_interface_det(3,i))
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
                     blk%zv(intfr(j)%vz_interface_det(2,i)), &
                     blk%zv(intfr(j)%vz_interface_det(3,i))
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
                      blk%xw(intfr(j)%wx_interface_det(2,i)), &
                      blk%xw(intfr(j)%wx_interface_det(3,i))
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
                      blk%yw(intfr(j)%wy_interface_det(2,i)), &
                      blk%yw(intfr(j)%wy_interface_det(3,i))
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
                      blk%zw(intfr(j)%wz_interface_det(2,i)), &
                      blk%zw(intfr(j)%wz_interface_det(3,i))
93329       FORMAT(A3,I5,3I5,3F8.5)
        end do
        end do
        ENDDO

        end subroutine change_block_interface

        SUBROUTINE cellCount_solid_coarse(blk)

         type(Blocks), intent(inout) :: blk

        INTEGER (int64) ::  n, iPt, iPt1, iPt2, i, j, k

        blk%fluidCellCount=0
         iPt  = 0
         iPt1 = 0
         iPt2 = 0
         DO k = 2, blk%nz +1
         DO j = 2, blk%ny +1
         DO i = 2, blk%nx +1
            IF (blk%cell(i,j,k)==0) THEN
                blk%fluidCellCount=blk%fluidCellCount +1
            ENDIF
         end do
         end do
         end do
         ALLOCATE (blk%fluidIndexPtr(blk%fluidCellCount,3))
         DO k = 2, blk%nz +1
         DO j = 2, blk%ny +1
         DO i = 2, blk%nx +1
            IF (blk%cell(i,j,k)==0) THEN
               iPt1 = iPt1 + 1
               blk%fluidIndexPtr(iPt1, 1) = i
               blk%fluidIndexPtr(iPt1, 2) = j
               blk%fluidIndexPtr(iPt1, 3) = k
            ENDIF
         END DO
         END DO
         END DO
         blk%redCellCount = 0
         blk%blackCellCount  = 0

         DO n = 1, blk%fluidCellCount
            i = blk%fluidIndexPtr(n, 1)
            j = blk%fluidIndexPtr(n, 2)
            k = blk%fluidIndexPtr(n, 3)
            IF (mod(i+j+k,2_int64)==1) THEN
               blk%redCellCount = blk%redCellCount + 1
            ELSE
               blk%blackCellCount = blk%blackCellCount + 1
            ENDIF
         ENDDO

         ALLOCATE (blk%redCellIndexPtr(blk%redCellCount,3), &
                   blk%blackCellIndexPtr(blk%blackCellCount,3))
         ipt1 = 0
         iPt = 0

         DO n = 1, blk%fluidCellCount
            i = blk%fluidIndexPtr(n, 1)
            j = blk%fluidIndexPtr(n, 2)
            k = blk%fluidIndexPtr(n, 3)
            IF (mod(i+j+k,2_int64)==1) THEN
               iPt = iPt + 1
               blk%redCellIndexPtr(iPt, 1) = i
               blk%redCellIndexPtr(iPt, 2) = j
               blk%redCellIndexPtr(iPt, 3) = k
            ELSE
               iPt1 = iPt1 + 1
               blk%blackCellIndexPtr(iPt1, 1) = i
               blk%blackCellIndexPtr(iPt1, 2) = j
               blk%blackCellIndexPtr(iPt1, 3) = k
            ENDIF
         ENDDO

            print*, 1, blk%fluidCellCount, blk%redCellCount, blk%blackCellCount

        END SUBROUTINE cellCount_solid_coarse
end module biocfd_search

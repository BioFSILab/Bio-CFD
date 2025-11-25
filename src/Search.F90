module biocfd_search
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64, int32
  use global, only: block, blk_start, totime, &
       pi, ita, dxmin, deltat, &
       re, freq, inor, &
       coarse_flcnt_check, intfr
  use biocfd_fine_interp, only: fineUpdate_mv
  use biocfd_fine_interp_bound, only : fineUpdate_bd_mv
  use biocfd_interface_detail, only : print_interface_detail
  use biocfd_block_type, only: Block_t
  use biocfd_interface_type, only: Interface_t
  implicit NONE

  private
  real(dp), parameter :: xfact = 0.05_dp

  public :: findDistnode, shiftSurfaceNodesInitial, computeSurfaceNorm
  public :: tagging_th, tagging_th_move, block_move_check, cellcount_solid
  public :: cellcount_solid_coarse, change_block_coords
  public :: change_block_interface, computenormdistance, computesurfacevariables, findtscells
  public :: fine_block_cell, selectiveretagging_th, change_block_coords_interfaces

  contains
    SUBROUTINE findDistnode(blk)
      type(Block_t), intent(inout) :: blk
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
                end if
        END IF
        IF (blk%ibNodeId(i)==51) THEN
                if (abs(blk%xnode(i)) < dist1)then
                        dist1=blk%xnode(i)
                        blk%mkx1=i
                end if
                if (abs(blk%xnode(i)) > dist2)then
                        dist2=blk%xnode(i)
                        blk%mkx2=i
                end if
        END IF
        END DO
        end subroutine findDistnode

        SUBROUTINE shiftSurfaceNodesInitial(blk,aoa1,aoa2,piv_pt)
        type(Block_t), intent(inout) :: blk
        real(dp),intent(in) :: aoa1, aoa2, piv_pt
        INTEGER(int64) ::  i
        REAL(dp)      ::  xr1, yr1, zr1, angt
        REAL(dp)      :: bdy,bdfr

        ALLOCATE (blk%xnode1(blk%ibNodes), blk%ynode1(blk%ibNodes), &
                  blk%znode1(blk%ibNodes))
        blk%xnode1 = blk%xnode
        blk%ynode1 = blk%ynode
        blk%znode1 = blk%znode
        bdfr=15
        bdy=15*dxmin
        angt  =  2._dp*pi*bdfr

        blk%xmove = 0.
        blk%ymove = 0.
        blk%zmove = 0.

        blk%a0 = blk%a0*pi/180_dp
         angt  =  2._dp*pi*blk%bfreq
        blk%xpth1=blk%xshift-(ita*dxmin*xfact)
        blk%xpth2=blk%xshift-(ita*dxmin*xfact)
        blk%ypth2=(blk%yamp)*sin(angt*blk%xshift)
        blk%piv_x = blk%xshift
        blk%piv_y = blk%yshift
        blk%piv_z = blk%zshift

        blk%thetaDot  =  0.
        blk% thetaDDot =  0._dp
        blk% thetaDot   =  0.
        blk%thetaDot1  = 0.
        blk% thetaDot2  = 0.
        blk% alphaDot  = 0.
        blk% thetaDDot1 = 0.
        blk%thetaDDot2 = 0.
        blk% thetaDDot  = 0.
        blk%yt         =  bdy*sin(2*pi*bdfr*totime)
        blk%ydot       =  angt*bdy*cos(2*pi*bdfr*totime)
        blk%yddot      =  -angt*angt*bdy*sin(2*pi*bdfr*totime)
        blk%xt         =  blk%xshift- (ita*dxmin*xfact)
        blk%xdot       =  -(dxmin*xfact)/deltat
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
         ELSE IF (blk%ibNodeId(i)==52) THEN
            xr1 =  blk%xnode(i)
            zr1 =  blk%znode(i)*cos(aoa2) + blk%ynode(i)*sin(aoa2)  + piv_pt &
                   - piv_pt*cos(aoa2)
            yr1 = -blk%znode(i)*sin(aoa2) + blk%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
           ELSE
               xr1 = blk% xnode(i)
               zr1 = blk% znode(i)  !*cos(aoa1) + ynode(i)*sin(aoa1)
               yr1 = blk% ynode(i)  !*sin(aoa1) + ynode(i)*cos(aoa1)
           END IF
           blk%xnode1(i) = xr1+ blk%xShift
           blk%ynode1(i) = yr1+ blk%yShift
           blk%znode1(i) = zr1+ blk%zshift
        END DO

      END SUBROUTINE shiftSurfaceNodesInitial

      SUBROUTINE computeSurfaceVariables(blk,g,phase_angle,piv_pt)
        type(Block_t), intent(inout) :: blk
        real(dp),intent(in) :: phase_angle,piv_pt
        real(dp) :: aoa1,aoa2
        INTEGER(int64) ::  i
        INTEGER(int64), intent(in) :: g
        REAL(dp)      ::  xr1, yr1, zr1
        REAL(dp)      :: angg, angt
        REAL(dp)      :: bdy,bdfr
        REAL(dp)      :: ang_theta
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
        blk%xpth1=blk%xpth2
        PRINT*, "angles =", aoa1*180._dp/pi, aoa2*180._dp/pi
        blk%thetaDot1 = ang_theta*blk%a0*cos(2._dp*pi*freq*(totime+deltat) + phase_angle)
        blk%thetaDDot1 = -ang_theta*ang_theta*blk%a0*sin(2._dp*pi*freq*(totime+deltat) &
                              + phase_angle)
        blk%thetaDot2  = -blk%thetaDot1
        blk%thetaDDot2 = -blk%thetaDDot1

        blk%yt         =  bdy*sin(angt*(totime))
        blk%ydot       =  angt*bdy*cos(angt*(totime))
        blk%yddot      =  -angt*angt*bdy*sin(angt*(totime))

        blk%xt         = blk%xpth2
        blk%xdot       = -(dxmin*xfact)/deltat

        blk%ychg=blk%yt - bdy*sin(angt*(totime-deltat))
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
      1  FORMAT("d",I4.4,"_index.",F8.2,".",i3.1,".dat")
        OPEN(UNIT = 17, FILE = filename1,POSITION="APPEND", STATUS = "unknown")
        write(17,14) totime, blk%nxty_cent, blk%inity_cent, blk%ymove, &
                     blk%nxtx_cent, blk%xmove
     close(17)
     14      FORMAT(7F15.8)
        print*,"centn",blk%nxty_cent,"centi",blk%inity_cent,"mv",blk%ychg
      DO i = 1, blk%ibnodes
      IF (blk%ibNodeId(i)==51) THEN
             xr1 =  blk%xnode(i)
             zr1 =  blk%znode(i)*cos(aoa1) + blk%ynode(i)*sin(aoa1) + piv_pt &
                    - piv_pt*cos(aoa1)
             yr1 = -blk%znode(i)*sin(aoa1) + blk%ynode(i)*cos(aoa1) + piv_pt*sin(aoa1)
         ELSE IF (blk%ibNodeId(i)==52) THEN
             xr1 =  blk%xnode(i)
             zr1 =  blk%znode(i)*cos(aoa2) + blk%ynode(i)*sin(aoa2)  + piv_pt &
                    - piv_pt*cos(aoa2)
             yr1 = -blk%znode(i)*sin(aoa2) + blk%ynode(i)*cos(aoa2)  + piv_pt*sin(aoa2)
           ELSE
               xr1 = blk% xnode(i)
               zr1 = blk% znode(i)
               yr1 = blk% ynode(i)
           END IF

           blk%xnode1(i) = xr1 +blk%xpth2
           blk%ynode1(i) = yr1 +blk%yshift +blk%ymove
           blk%znode1(i) = zr1 +blk%zshift +blk%zmove
        END DO
           write(*,*) blk%xnode1(1),blk%ynode1(1),blk%znode1(1)

      END SUBROUTINE computeSurfaceVariables

      SUBROUTINE computeSurfaceNorm(blk)
        type(Block_t), intent(inout) :: blk
        INTEGER(int64) ::  n  !c1, c2, c3, c4
        REAL(dp)      :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL
        REAL(dp)      :: var_xcent, var_ycent, var_zcent

        ALLOCATE (blk%xcent(blk%ibElems), blk%ycent(blk%ibElems), &
                  blk%zcent(blk%ibElems), &
                  blk%cosAlpha(blk%ibElems), blk%cosBeta(blk%ibElems), &
                  blk%cosGamma(blk%ibElems))

        !compute centroid and direction cosines
       !$acc parallel loop gang vector default(present) &
       !$acc private (var_xcent, var_ycent, var_zcent,p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL) &
       !$acc firstprivate (inor)
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

           blk%cosAlpha(n) = inor * blk%cosAlpha(n) / lenEl  ! direction cosine unit normal along x
           blk%cosBeta(n)  = inor * blk%cosBeta(n) / lenEl   ! direction cosine unit normal along y
           blk%cosGamma(n) = inor * blk%cosGamma(n) / lenEl  ! direction cosine unit normal along z
        END DO
       !$acc end parallel loop

        print*, "SurfaceNorm done, inor =", inor

     END SUBROUTINE computeSurfaceNorm

SUBROUTINE tagging_th_core(blk)
       type(Block_t), intent(inout) :: blk
       INTEGER(int64) :: m, i, j, k,  nel2Cen, nel2Pnt, sumNodeId
        REAL(dp)      :: minDis1, minDis, &
                         n2dotn, dis_cen, dis_pnt

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
        !$omp& shared(blk)
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
               END IF
               IF (dis_pnt<minDis1) THEN
                  minDis1   = dis_pnt
                  nel2Pnt   = m
               END IF
            END DO
            IF((blk%x1(i)<=blk%xcent(nel2Cen) .AND. &
                blk%x1(i+1)>=blk%xcent(nel2Cen)).AND. &
               (blk%y1(j)<=blk%ycent(nel2Cen) .AND. &
                blk%y1(j+1)>=blk%ycent(nel2Cen)).AND. &
               (blk%z1(k)<=blk%zcent(nel2Cen) .AND. &
                blk%z1(k+1)>=blk%zcent(nel2Cen))) THEN
               blk%cell(i,j,k) = 2

            END IF

            n2dotn  = (blk%x1(i) - blk%xcent(nel2Pnt))*blk%cosAlpha(nel2Pnt) + &
                      (blk%y1(j) - blk%ycent(nel2Pnt))*blk%cosBeta(nel2Pnt)  + &
                      (blk%z1(k) - blk%zcent(nel2Pnt))*blk%cosGamma(nel2Pnt)

            IF (n2dotn>=-1e-16_dp) THEN
               blk%nodeIdTag(i,j,k) = 0
            ELSE
               blk%nodeIdTag(i,j,k) = 1
            END IF
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
                  END IF
               END IF
           END DO
           END DO
           END DO
!$acc end parallel loop

         blk%ibCellCount = 0
         blk%solidCellCount = 0
         blk%fluidCellCount = 0
         DO k = 2, blk%nz+1
         DO j = 2, blk%ny+1
         DO i = 2, blk%nx+1
               IF (blk%cell(i,j,k)==1) THEN
                  blk%solidCellCount = blk%solidCellCount + 1
               ELSE IF (blk%cell(i,j,k)==0) THEN
                  blk%fluidCellCount  = blk%fluidCellCount + 1
               ELSE IF (blk%cell(i,j,k)==2) THEN
                  blk%ibCellCount = blk%ibCellCount + 1
               END IF
         END DO
         END DO
         END DO

         print*, "search done"
         Print*, "imms. cells=", blk%ibCellCount
        Print*, "fluid cells=",blk%fluidCellCount
        Print*, "solid cells=", blk%solidCellCount
     END SUBROUTINE tagging_th_core

     SUBROUTINE tagging_th(blk,blk_no)
       type(Block_t), intent(inout) :: blk
       INTEGER(int64), intent(in)  :: blk_no
       CHARACTER(LEN=120) :: filename1
       INTEGER(int64) :: i, j, k

        call tagging_th_core(blk)
               WRITE(filename1,1) blk_no
  1      FORMAT("butter_f.",i3.3,".dat")
          OPEN(11,FILE=filename1,status="unknown")
        DO k = 1, blk%nz+2
        DO j = 1, blk%ny+2
        DO i = 1, blk%nx+2
        WRITE(11,*) blk%cell(i,j,k), blk%nodeIdTag(i,j,k)
        END DO
        END DO
        END DO
        CLOSE(11)

         WRITE(filename1,2) blk_no
 2       FORMAT("butter_cellcount_f.",i3.3,".dat")
         OPEN(12,FILE=filename1,FORM="formatted")
        WRITE(12,*) blk%solidCellCount, blk%fluidCellCount, blk%ibCellCount
        CLOSE(12)

     END SUBROUTINE tagging_th

     SUBROUTINE tagging_th_move(blk, id)

        type(Block_t), intent(inout) :: blk
        INTEGER(int64), intent(in) :: id
        integer(int64) :: g

        if (blk%move_check /= 1) return
        call tagging_th_core(blk)

        DO g=1, size(intfr)
          ! Find the interface which this block corresponds to. We
          ! check if this is the right interface by skipping over
          ! interfaces where the b_blk is not this one.
          if (intfr(g)%b_blk /= id) cycle
          call fineUpdate_mv(g)
          call fineUpdate_bd_mv(g)
        END DO
     END SUBROUTINE tagging_th_move

     SUBROUTINE findTScells(blk)
       type(Block_t), intent(inout) :: blk
        INTEGER            :: i, j, k, i1, j1, k1, iPt1, m, n, tscnt
        integer :: idx
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
         IF(blk%ibSurfID(blk%nelp(n))==51.OR.blk%ibSurfID(blk%nelp(n))==52) then
           blk%cell2(i, j, k) = 2
         end if

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
                     END IF
                 END DO
                 END DO
                 END DO
        !$acc end parallel loop

        blk%TSCellCount = tscnt
        print*, "TScell count =", blk%TSCellCount

        ALLOCATE(blk%TSIndexPtr(blk%TSCellCount,3))

        iPt1 = 0
        !$acc parallel loop collapse(3) private(idx)
        DO k=1,blk%nz+3
           DO j=1,blk%ny+3
              DO i=1,blk%nx+3
                 IF (blk%cell2(i,j,k)==2) THEN
                    !$acc atomic capture
                    iPt1 = iPt1 + 1
                    idx = iPt1
                    !$acc end atomic
                    blk%TSIndexPtr(idx, :) = [i, j, k]
                 END IF
              END DO
           END DO
        END DO
        !$acc end parallel loop

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
              END IF
           END DO
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
        END DO
        !$acc end parallel loop
             END SUBROUTINE findTScells

     SUBROUTINE selectiveRetagging_th(blk)
       type(Block_t), intent(inout) :: blk
        INTEGER(int64) ::  n, m, i, j, k, i1, j1, k1, nn, &
                               nel2Pnt, nel2Cen, sumNodeID
        INTEGER            :: flcnt, sdcnt, ibcnt
        REAL(dp)      :: minDis, minDis1, dis_cen, dis_pnt, n2dotn
        integer :: iprime, jprime, kprime

        if(blk%blk_mv_tag ==0)then

       ! Set the intercepted indicies cell value to 0, we do this in a
       ! seperate loop so that we can nicely GPU-ise the computation
       !$acc parallel loop default(present) private(i1, j1, k1)
       DO nn = 1, blK%ibCellCount
         i1 = blk%interceptedIndexPtr(nn, 1)
         j1 = blk%interceptedIndexPtr(nn, 2)
         k1 = blk%interceptedIndexPtr(nn, 3)
         blk%cell(i1,j1,k1) = 0
        END DO
        !$acc end parallel loop

        !$acc parallel loop default(present) collapse(4) private(i, j, k) &
        !$acc private(minDis, minDis1, dis_cen, dis_pnt, nel2Cen, nel2Pnt) &
        !$acc private(n2dotn)
        DO nn = 1, blk%ibCellCount
           DO kprime = -1,+1
           DO jprime = -1,+1
           DO iprime = -1,+1

               i = blk%interceptedIndexPtr(nn, 1) + iprime
               j = blk%interceptedIndexPtr(nn, 2) + jprime
               k = blk%interceptedIndexPtr(nn, 3) + kprime

               minDis  = 1e14_dp
               minDis1 = 1e14_dp

               !$acc loop seq
               DO m = 1, blk%ibElems
                  ! No need to take sqrt because we just use for distance comparison
                  dis_cen  = (blk%xp(i)-blk%xcent(m))**2 &
                           + (blk%yp(j)-blk%ycent(m))**2 &
                           + (blk%zp(k)-blk%zcent(m))**2
                  dis_pnt  = (blk%x1(i)-blk%xcent(m))**2 &
                           + (blk%y1(j)-blk%ycent(m))**2 &
                           + (blk%z1(k)-blk%zcent(m))**2
                  IF (dis_cen<minDis) THEN
                     minDis    = dis_cen
                     nel2Cen   = m
                  END IF
                  IF (dis_pnt<minDis1) THEN
                     minDis1   = dis_pnt
                     nel2Pnt   = m
                  END IF
               END DO
               IF((blk%x1(i)<=blk%xcent(nel2Cen).AND. &
                   blk%x1(i+1)>=blk%xcent(nel2Cen)).AND. &
                  (blk%y1(j)<=blk%ycent(nel2Cen).AND. &
                   blk%y1(j+1)>=blk%ycent(nel2Cen)).AND. &
                  (blk%z1(k)<=blk%zcent(nel2Cen).AND. &
                   blk%z1(k+1)>=blk%zcent(nel2Cen))) THEN
                        blk%cell(i,j,k) = 2
               END IF

               n2dotn  = (blk%x1(i) - blk%xcent(nel2Pnt)) * blk%cosAlpha(nel2Pnt) + &
                         (blk%y1(j) - blk%ycent(nel2Pnt)) * blk%cosBeta(nel2Pnt)  + &
                         (blk%z1(k) - blk%zcent(nel2Pnt)) * blk%cosGamma(nel2Pnt)

               IF (n2dotn>=-1e-16_dp) THEN
                  blk%nodeIdTag(i,j,k) = 0
               ELSE
                  blk%nodeIdTag(i,j,k) = 1
               END IF
           END DO
           END DO
           END DO
        END DO
!$acc end parallel loop

!$acc parallel loop gang vector default(present)
        DO nn = 1, blk%ibCellCount
        i1 = blk%interceptedIndexPtr(nn, 1)
        j1 = blk%interceptedIndexPtr(nn, 2)
        k1 = blk%interceptedIndexPtr(nn, 3)
           !$acc loop collapse(3) seq
           DO k = k1-1, k1+1
           DO j = j1-1, j1+1
           DO i = i1-1, i1+1
               IF (blk%cell(i,j,k)/=2) THEN
                  sumNodeId = 0
                  sumNodeId = blk%nodeIdTag(i,j,k)      + blk%nodeIdTag(i+1,j,k)     &
                           + blk%nodeIdTag(i,j+1,k)    + blk%nodeIdTag(i+1,j+1,k)   &
                           + blk%nodeIdTag(i,j,k+1)    + blk%nodeIdTag(i+1,j,k+1)     &
                           + blk%nodeIdTag(i,j+1,k+1)  + blk%nodeIdTag(i+1,j+1,k+1)
                  IF (sumNodeId==8) THEN
                     blk%cell(i,j,k) = 1
                  ELSE
                     blk%cell(i,j,k) = 0
                  END IF
               END IF
           END DO
           END DO
           END DO
        END DO
blk%ibCellCount = 0
blk%solidCellCount = 0
blk%fluidCellCount = 0
sdcnt=0
flcnt=0
ibcnt=0
!$acc parallel loop gang vector collapse(3) default(present) private(i,j,k,n) reduction(+: sdcnt, flcnt, ibcnt)
         DO k = 2, blk%nz+1
         DO j = 2, blk%ny+1
         DO i = 2, blk%nx+1
            n = i-1  + blk%nx*(j-2)  + blk%nx*blk%ny*(k-2)
            IF (blk%cell(i,j,k)==1) THEN
                sdcnt = sdcnt + 1
            ELSE IF (blk%cell(i,j,k)==0) THEN
               flcnt  = flcnt + 1
            ELSE IF (blk%cell(i,j,k)==2) THEN
                ibcnt = ibcnt + 1
            END IF
         END DO
         END DO
         END DO
!$acc end parallel loop
blk%ibCellCount = ibcnt
blk%solidCellCount = sdcnt
blk%fluidCellCount = flcnt
       print*, "selective retagging", blk%ibCellCount, blk%fluidCellCount, &
                blk%solidCellCount
     END IF
     END SUBROUTINE selectiveRetagging_th

     SUBROUTINE cellCount_solid(blk,blk_no)

       type(Block_t), intent(inout) :: blk
       integer(int64), intent(in) :: blk_no
       INTEGER (int64) ::  n, iPt, iPt1, iPt2, i, j, k
       INTEGER (int64) :: cell_val
       integer :: red_count, black_count
       integer (int64) :: idx

       iPt  = 0
       iPt1 = 0
       iPt2 = 0

       ALLOCATE(blk%interceptedIndexPtr(blk%ibCellCount,3), &
            blk%fluidIndexPtr(blk%fluidCellCount, 3), &
            blk%solidIndexPtr(blk%solidCellCount, 3))

       !$acc parallel loop collapse(3) private(cell_val, idx)
       DO k = 2, blk%nz+1
          DO j = 2, blk%ny+1
             DO i = 2, blk%nx+1

                cell_val = blk%cell(i,j,k)

                IF (cell_val==0) THEN
                   !$acc atomic capture
                   iPt1 = iPt1 + 1
                   idx = iPt1
                   !$acc end atomic
                   blk%fluidIndexPtr(idx, :) = [i, j, k]
                ELSE IF (cell_val==1) THEN
                   !$acc atomic capture
                   iPt2 = iPt2 + 1
                   idx = iPt2
                   !$acc end atomic
                   blk%solidIndexPtr(idx, :) = [i, j, k]
                ELSE IF (cell_val==2) THEN
                   !$acc atomic capture
                   iPt = iPt + 1
                   idx = iPt
                   !$acc end atomic
                   blk%interceptedIndexPtr(idx, :) = [i, j, k]
                END IF

             END DO
          END DO
       END DO
       !$acc end parallel loop

       red_count = 0
       black_count = 0

       !$acc parallel loop reduction(+:red_count,black_count) private(i, j, k)
       DO n = 1, blk%fluidCellCount
          if (mod(sum(blk%fluidIndexPtr(n, :)), 2_int64) == 1) then
             red_count = red_count + 1
          ELSE
             black_count = black_count + 1
          END IF
       END DO
       !$acc end parallel loop

       blk%redCellCount = red_count
       blk%blackCellCount  = black_count

       ALLOCATE (blk%redCellIndexPtr(blk%redCellCount,3), &
            blk%blackCellIndexPtr(blk%blackCellCount,3))
       ipt1 = 0
       iPt = 0
       !$acc parallel loop private(idx)
       DO n = 1, blk%fluidCellCount
          if (mod(sum(blk%fluidIndexPtr(n, :)), 2_int64) == 1) then
             !$acc atomic capture
             iPt = iPt + 1
             idx = iPt
             !$acc end atomic
             blk%redCellIndexPtr(idx, :) = blk%fluidIndexPtr(n, :)
          ELSE
             !$acc atomic capture
             iPt1 = iPt1 + 1
             idx = iPt1
             !$acc end atomic
             blk%blackCellIndexPtr(idx, :) = blk%fluidIndexPtr(n, :)
          END IF
       END DO
       !$acc end parallel loop
       print*, "cellCount:", blk_no, block(blk_no)%fluidCellCount, &
            block(blk_no)%redCellCount, block(blk_no)%blackCellCount
     END SUBROUTINE cellCount_solid

     SUBROUTINE computeNormDistance(blk)
       type(Block_t), intent(inout) :: blk
        INTEGER            ::  nel2u1, nel2u2, nel2v1, nel2v2, nel2w1, nel2w2
        INTEGER            :: k, nel2p, ibxx, m
        REAL(dp) :: n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, &
                    dis, dis1, dis2, dis3, dis4, dis5, dis6, &
                    minDis, minDis1, minDis2, minDis3, minDis4, minDis5, minDis6, &
                    cent_x, cent_y,cent_z

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
              END IF
              IF (dis1<minDis1) THEN
                 minDis1   = dis1
                 nel2u1    = m
              END IF
              IF (dis2<minDis2) THEN
                 minDis2   = dis2
                 nel2u2    = m
              END IF
              IF (dis3<minDis3) THEN
                 minDis3   = dis3
                 nel2v1    = m
              END IF
              IF (dis4<minDis4) THEN
                 minDis4   = dis4
                 nel2v2    = m
              END IF
              IF (dis5<minDis5) THEN
                 minDis5   = dis5
                 nel2w1    = m
              END IF
              IF (dis6<minDis6) THEN
                 minDis6   = dis6
                 nel2w2    = m
              END IF
           END DO

           blk%nelp(k)  = nel2p
           blk%nelu1(k) = nel2u1
           blk%nelu2(k) = nel2u2
           blk%nelv1(k) = nel2v1
           blk%nelv2(k) = nel2v2
           blk%nelw1(k) = nel2w1
           blk% nelw2(k) = nel2w2
          blk%pNormDis(k)  = (n1x -blk%xcent(nel2p))*blk%cosAlpha(nel2p) + &
                   (n1y -blk%ycent(nel2p))*blk%cosBeta(nel2p)  + &
                   (n1z -blk%zcent(nel2p))*blk%cosGamma(nel2p)

          blk%u1NormDis(k) = (n2x -blk%xcent(nel2u1))*blk%cosAlpha(nel2u1) + &
                       (n1y -blk%ycent(nel2u1))*blk%cosBeta(nel2u1)  + &
                       (n1z - blk%zcent(nel2u1))*blk%cosGamma(nel2u1)

          blk%u2NormDis(k) = (n3x - blk%xcent(nel2u2))*blk%cosAlpha(nel2u2) + &
              (n1y -blk%ycent(nel2u2))*blk%cosBeta(nel2u2)  + &
              (n1z -blk%zcent(nel2u2))*blk%cosGamma(nel2u2)

          blk%v1NormDis(k) = (n1x - blk%xcent(nel2v1))*blk%cosAlpha(nel2v1) + &
              (n2y - blk%ycent(nel2v1))*blk%cosBeta(nel2v1)  + &
              (n1z - blk%zcent(nel2v1))*blk%cosGamma(nel2v1)

          blk%v2NormDis(k) = (n1x - blk%xcent(nel2v2))*blk%cosAlpha(nel2v2) + &
                   (n3y -blk%ycent(nel2v2))*blk%cosBeta(nel2v2)  + &
                   (n1z -blk%zcent(nel2v2))*blk%cosGamma(nel2v2)

          blk%w1NormDis(k) = (n1x -blk%xcent(nel2w1))*blk%cosAlpha(nel2w1) + &
                       (n1y -blk%ycent(nel2w1))*blk%cosBeta(nel2w1)  + &
                       (n2z -blk%zcent(nel2w1))*blk%cosGamma(nel2w1)

          blk%w2NormDis(k) = (n1x -blk%xcent(nel2w2))*blk%cosAlpha(nel2w2) + &
                           (n1y -blk%ycent(nel2w2))*blk%cosBeta(nel2w2)  + &
                           (n3z -blk% zcent(nel2w2))*blk%cosGamma(nel2w2)
        END DO
        !$acc end parallel loop

     END SUBROUTINE computeNormDistance

        SUBROUTINE fine_block_cell

        INTEGER(int64) :: i, j, k, g, factor, a_blk_no, b_blk_no

        !$acc parallel present(block)
        block(1)%cell_n = 0
        block(1)%cell = 0
        block(1)%cell_pr = 0
        !$acc end parallel

        DO g=1, size(intfr)
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk
           factor=intfr(g)%b_msh/intfr(g)%a_msh

        !$acc parallel loop collapse(3) default(present)
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

        end if
        if (block(a_blk_no)%xp(i) >= intfr(g)%xintf_start .and. &
             block(a_blk_no)%xp(i) <= intfr(g)%xintf_end   .and. &
             block(a_blk_no)%zp(k) >= intfr(g)%zintf_start .and. &
             block(a_blk_no)%zp(k) <= intfr(g)%zintf_end   .and. &
             block(a_blk_no)%yp(j) >= intfr(g)%yintf_start .and. &
             block(a_blk_no)%yp(j) <= intfr(g)%yintf_end)then

                block(a_blk_no)%cell_pr(i,j,k)=1

        end if
        END DO
        END DO
        END DO
        !$acc end parallel loop
        END DO


        end subroutine fine_block_cell

        SUBROUTINE block_move_check

        INTEGER(int64) :: i,j,k,g, a_blk_no, b_blk_no, factor
        REAL(dp) :: ydisp1, xdisp1, zdisp1, mg1
        REAL(dp) :: marginx, marginy, marginz, yval_up, yval_dw, xval_lt, xval_rt

        DO g=1,size(intfr)
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

        ! Here we check if array is allocated. This is important for
        ! MPI, as the arrays may not be allocated on this rank. When
        ! not running in MPI mode, this should always be allocated.
        if (.not. allocated(block(b_blk_no)%u)) cycle

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
        print*,"blk_check_cond", marginx, marginy
        print*,"ydsip",abs(yval_dw-intfr(g)%yintf_start),abs(yval_up - intfr(g)%yintf_end)
        print*,block(b_blk_no)%xnode1(block(b_blk_no)%mkx1),block(b_blk_no)%nxtx_cent
        print*,block(b_blk_no)%ynode1(block(b_blk_no)%mk),block(b_blk_no)%nxty_cent
        block(b_blk_no)%move_amty=0
        block(b_blk_no)%move_amtx=0
        block(b_blk_no)%move_amtz=0
        if ((abs(xval_lt - intfr(g)%xintf_start)      <= marginx)  .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start)      <= marginz)  .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end)        <= marginz)  .or. &
             ((abs(yval_dw - intfr(g)%yintf_start) <= marginy))   .or. &
             ((abs(yval_up - intfr(g)%yintf_end)  <= marginy)))then

                block(b_blk_no)%move_check=1
             block(b_blk_no)%blk_mv_tag=1.
                 coarse_flcnt_check=1
            if(abs(yval_dw - intfr(g)%yintf_start) <= marginy    .or. &
                (abs(yval_up - intfr(g)%yintf_end)  <= marginy))then

                block(b_blk_no)%move_amty = &
                  floor((block(b_blk_no)%nxty_cent -block(b_blk_no)%inity_cent)/block(b_blk_no)%dy)
                if (abs(block(b_blk_no)%move_amty) < factor)then
                        if (block(b_blk_no)%move_amty < 0)then
                        block(b_blk_no)%move_amty=-factor
                        else
                        block(b_blk_no)%move_amty=factor
                        end if

                else
                        if (block(b_blk_no)%move_amty < 0)then
                        block(b_blk_no)%move_amty= block(b_blk_no)%move_amty &
                                                   + mod(abs(block(b_blk_no)%move_amty),factor)
                        else
                        block(b_blk_no)%move_amty= block(b_blk_no)%move_amty &
                                                   - mod(abs(block(b_blk_no)%move_amty),factor)
                        end if

                end if
                block(b_blk_no)%inity_cent=block(b_blk_no)%nxty_cent

                end if
        if ((abs(xval_lt- intfr(g)%xintf_start) < marginx))then
                block(b_blk_no)%move_amtx = &
                  floor((block(b_blk_no)%nxtx_cent - block(b_blk_no)%initx_cent) &
                  / block(b_blk_no)%dx) - factor
                print*,"move_amtx",block(b_blk_no)%move_amtx
                if (abs(block(b_blk_no)%move_amtx) < factor)then
                        if (block(b_blk_no)%move_amtx < 0)then
                        block(b_blk_no)%move_amtx=-factor
                        else
                        block(b_blk_no)%move_amtx=factor
                        end if

                else
                        if (block(b_blk_no)%move_amtx < 0)then
                        block(b_blk_no)%move_amtx= block(b_blk_no)%move_amtx &
                                                   + mod(abs(block(b_blk_no)%move_amtx),factor)
                        else
                        block(b_blk_no)%move_amtx = block(b_blk_no)%move_amtx &
                                                    - mod(abs(block(b_blk_no)%move_amtx),factor)
                        end if

                end if
                block(b_blk_no)%initx_cent=block(b_blk_no)%nxtx_cent
                end if
            if((abs(block(b_blk_no)%piv_z - intfr(g)%zintf_start) < (0.500_dp-0.410_dp-mg1)) .or. &
             (abs(block(b_blk_no)%piv_z - intfr(g)%zintf_end)   < (0.600_dp-0.500_dp-mg1)))then
                block(b_blk_no)%move_amtz = &
                  floor((block(b_blk_no)%nxtz_cent -block(b_blk_no)%initz_cent)/block(b_blk_no)%dz)
                if (abs(block(b_blk_no)%move_amtz) < factor)then
                        if (block(b_blk_no)%move_amtz < 0)then
                        block(b_blk_no)%move_amtz=-factor
                        else
                        block(b_blk_no)%move_amtz=factor
                        end if

                else
                        if (block(b_blk_no)%move_amtz < 0)then
                        block(b_blk_no)%move_amtz = block(b_blk_no)%move_amtz &
                                                    + mod(abs(block(b_blk_no)%move_amtz),factor)
                        else
                        block(b_blk_no)%move_amtz = block(b_blk_no)%move_amtz &
                                                    - mod(abs(block(b_blk_no)%move_amtz),factor)
                        end if

                end if
                block(b_blk_no)%initz_cent=block(b_blk_no)%nxtz_cent
                end if
                 print*,"blk_movez",block(b_blk_no)%move_amtz
                 print*,"blk_movey",block(b_blk_no)%move_amty
                 print*,"blk_movex",block(b_blk_no)%move_amtx

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
       print*,"zstart_org",intfr(g)%zintf_start,"zend_org",intfr(g)%zintf_end
       print*,"ystart_org",intfr(g)%yintf_start,"yend_org",intfr(g)%yintf_end
       print*,"xstart_org",intfr(g)%xintf_start,"xend_org",intfr(g)%xintf_end

       intfr(g)%xintf_start=intfr(g)%xintf_start + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx)
       intfr(g)%xintf_end=intfr(g)%xintf_end + (block(b_blk_no)%move_amtx * block(b_blk_no)%dx)
       intfr(g)%yintf_start=intfr(g)%yintf_start + (block(b_blk_no)%move_amty * block(b_blk_no)%dy)
       intfr(g)%yintf_end=intfr(g)%yintf_end + (block(b_blk_no)%move_amty * block(b_blk_no)%dy)
       intfr(g)%zintf_start=intfr(g)%zintf_start + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz)
       intfr(g)%zintf_end=intfr(g)%zintf_end + (block(b_blk_no)%move_amtz * block(b_blk_no)%dz)
       print*,"zstart_mv",intfr(g)%zintf_start,"zend_mv",intfr(g)%zintf_end
       print*,"ystart_mv",intfr(g)%yintf_start,"yend_mv",intfr(g)%yintf_end
       print*,"xstart_mv",intfr(g)%xintf_start,"xend_mv",intfr(g)%xintf_end
       print*,"zstart_nw",intfr(g)%zintf_st_new,"zend_nw",intfr(g)%zintf_en_new
       print*,"ystart_nw",intfr(g)%yintf_st_new,"yend_nw",intfr(g)%yintf_en_new
       print*,"xstart_nw",intfr(g)%xintf_st_new,"xend_nw",intfr(g)%xintf_en_new

        DO i=1,block(b_blk_no)%nx+2
           if(block(b_blk_no)%xp(i) > intfr(g)%xintf_st_new)then
                block(b_blk_no)%cpy_x_start=i
                print*,"xpst",block(b_blk_no)%xp(i),"i",i
                exit
          end if
        END DO

        DO i=block(b_blk_no)%cpy_x_start,block(b_blk_no)%nx+2
           if(block(b_blk_no)%xp(i) > intfr(g)%xintf_en_new)then
                block(b_blk_no)%cpy_x_end=i-1
                print*,"xpen",block(b_blk_no)%xp(i-1),"i",i-1
                exit
          end if
        END DO
        DO j=1,block(b_blk_no)%ny+2
           if(block(b_blk_no)%yp(j) > intfr(g)%yintf_st_new)then
                block(b_blk_no)%cpy_y_start=j
                print*,"ypst",block(b_blk_no)%yp(j),"j",j
                exit
          end if
        END DO

        DO j=block(b_blk_no)%cpy_y_start,block(b_blk_no)%ny+2
           if(block(b_blk_no)%yp(j) > intfr(g)%yintf_en_new)then
                block(b_blk_no)%cpy_y_end=j-1
                print*,"ypen",block(b_blk_no)%yp(j-1),"j",j-1
                exit
          end if
        END DO

        DO k=1,block(b_blk_no)%nz+2
           if(block(b_blk_no)%zp(k) > intfr(g)%zintf_st_new)then
                block(b_blk_no)%cpy_z_start=k
                print*,"zpst",block(b_blk_no)%zp(k),"k",k
                exit
          end if
        END DO

        DO k=block(b_blk_no)%cpy_z_start,block(b_blk_no)%nz+2
           if(block(b_blk_no)%zp(k) > intfr(g)%zintf_en_new)then
                block(b_blk_no)%cpy_z_end=k-1
                print*,"zpen",block(b_blk_no)%zp(k-1),"k",k-1
                exit
          end if
        END DO

        DO k=1,block(b_blk_no)%nz+2
        DO i=1,block(b_blk_no)%nx+2
        DO j=1,block(b_blk_no)%ny+2
            block(b_blk_no)%u_dum(i,j,k)=block(b_blk_no)%u(i,j,k)
            block(b_blk_no)%v_dum(i,j,k)=block(b_blk_no)%v(i,j,k)
            block(b_blk_no)%w_dum(i,j,k)=block(b_blk_no)%w(i,j,k)
            block(b_blk_no)%p_dum(i,j,k)=block(b_blk_no)%p(i,j,k)
        END DO
        END DO
        END DO
        end if

        END DO
        end subroutine block_move_check

      ! This subroutine used to be combined with change_block_coords_interfaces
      ! and as a result this subroutine must be called just before
      ! change_block_coords
      SUBROUTINE change_block_coords(blk)
        type(Block_t), intent(inout) :: blk
        INTEGER(int64) :: i
        REAL(dp) :: change_y_f,change_x_f,change_z_f

        if (blk%move_check /= 1) return

        change_z_f= blk%move_amtz*blk%dz
        change_y_f= blk%move_amty*blk%dy
        change_x_f= blk%move_amtx * blk%dx
        print*, "chz",blk%move_amtz , blk%dz
        print*, "chy",blk%move_amty , blk%dy
        print*, "chx",blk%move_amtx , blk%dx

        DO i = 1, blk%nx+3
           blk%x1(i) = blk%x1(i) + change_x_f
        END DO

        DO i = 1, blk%ny+3
           blk%y1(i) = blk%y1(i) + change_y_f
        END DO

        DO i = 1, blk%nz+3
           blk%z1(i) = blk%z1(i)+ change_z_f
        END DO

        DO i = 1, blk%nx+3
           blk%xu(i) = blk%x1(i)
        END DO

        DO i = 1, blk%ny+3
           blk%yv(i) = blk%y1(i)
        END DO

        DO i = 1, blk%nz+3
           blk%zw(i) = blk%z1(i)
        END DO

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

      end subroutine change_block_coords

      ! This subroutine used to be combined with change_block_coords
      ! and as a result this subroutine must be called straight
      ! after change_block_coords
      subroutine change_block_coords_interfaces(local_intfr,blk)
        type(Interface_t),intent(in) :: local_intfr
        type(Block_t), intent(inout) :: blk
        INTEGER(int64) :: k,j,i,countx_st,countz_st,county_st

        ! Here we check if array is allocated. This is important for
        ! MPI, as the arrays may not be allocated on this rank. When
        ! not running in MPI mode, this should always be allocated.
         if (.not. allocated(blk%u)) return

        if (blk% move_check /= 1) return
        blk%cell_n=0
        DO k=1,blk%nz+2
        DO j=1,blk%ny+2
        DO i=1,blk%nx+2
           if(blk%xp(i) > (local_intfr%xintf_st_new  + blk%dx) .and. &
              blk%xp(i) < (local_intfr%xintf_en_new  - blk%dx).and. &
              blk%zp(k) > (local_intfr%zintf_st_new  + blk%dz) .and. &
              blk%zp(k) < (local_intfr%zintf_en_new  - blk%dz).and. &
              blk%yp(j) > (local_intfr%yintf_st_new  + blk%dy).and. &
              blk%yp(j) < (local_intfr%yintf_en_new  - blk%dy))then
                 blk%cell_n(i,j,k)=1
           end if
        END DO
        END DO
        END DO

        DO i=1,blk%nx+2
           if(blk%xp(i) > local_intfr%xintf_st_new)then
                 blk%cpy_x_start_mv=i
                 print*,"xpst_f",blk%xp(i),"i",i
                 exit
           end if
        END DO

        DO i=blk%cpy_x_start_mv,blk%nx+2
           if(blk%xp(i) > local_intfr%xintf_en_new)then
                 blk%cpy_x_end_mv=i-1
                 print*,"xpen_f",blk%xp(i-1),"i",i-1
                 exit
           end if
        END DO

        DO j=1,blk%ny+2
           if(blk%yp(j) > local_intfr%yintf_st_new)then
                  blk%cpy_y_start_mv=j
                  print*,"ypst_f",blk%yp(j),"j",j
                  exit
           end if
        END DO

        DO j=blk%cpy_y_start_mv,blk%ny+2
           if(blk%yp(j) > local_intfr%yintf_en_new)then
                 blk%cpy_y_end_mv=j-1
                 print*,"ypen_f",blk%yp(j-1),"j",j-1
                 exit
           end if
        END DO

        DO j=1,blk%nz+2
           if(blk%zp(j) > local_intfr%zintf_st_new)then
                 blk%cpy_z_start_mv=j
                 print*,"zpst_f",blk%zp(j),"k",j
                 exit
           end if
        END DO

        DO j=blk%cpy_z_start_mv,blk%nz+2
           if(blk%zp(j) > local_intfr%zintf_en_new)then
                 blk%cpy_z_end_mv=j-1
                 print*,"zpen_f",blk%zp(j-1),"k",j-1
                 exit
           end if
        END DO

        print*,"bef"
        countx_st=blk%cpy_x_start
        county_st=blk%cpy_y_start
        countz_st=blk%cpy_z_start

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
              END DO
              countx_st=countx_st+1
           END DO
           countz_st=countz_st+1
        END DO
        print*,"aft"

     end subroutine change_block_coords_interfaces

SUBROUTINE change_block_interface(local_intfr,blk_a,blk_b)

  type(Interface_t),intent(inout) :: local_intfr
  type(Block_t), intent(in) :: blk_a, blk_b
  INTEGER(int64) :: factor, increment

  ! Here we check if array is allocated. This is important for
  ! MPI, as the arrays may not be allocated on this rank. When
  ! not running in MPI mode, this should always be allocated.
  if (.not. allocated(blk_b%u)) return

  if (blk_b%move_check /= 1) return
  factor=local_intfr%b_msh/local_intfr%a_msh

  increment = blk_b%move_amtx/factor
  local_intfr%px_interface_det(1,1:local_intfr%counterxp) =&
          local_intfr%px_interface_det(1,1:local_intfr%counterxp) + increment

  local_intfr%ux_interface_det(1,1:local_intfr%counterxu) =&
          local_intfr%ux_interface_det(1,1:local_intfr%counterxu) + increment

  local_intfr%vx_interface_det(1,1:local_intfr%counterxv) =&
          local_intfr%vx_interface_det(1,1:local_intfr%counterxv) + increment

  local_intfr%wx_interface_det(1,1:local_intfr%counterxw) =&
          local_intfr%wx_interface_det(1,1:local_intfr%counterxw) + increment

  increment = blk_b%move_amty/factor
  local_intfr%py_interface_det(1,1:local_intfr%counteryp) =&
          local_intfr%py_interface_det(1,1:local_intfr%counteryp) + increment

  local_intfr%uy_interface_det(1,1:local_intfr%counteryu) =&
          local_intfr%uy_interface_det(1,1:local_intfr%counteryu) + increment

  local_intfr%vy_interface_det(1,1:local_intfr%counteryv) =&
          local_intfr%vy_interface_det(1,1:local_intfr%counteryv) + increment

  local_intfr%wy_interface_det(1,1:local_intfr%counteryw) =&
          local_intfr%wy_interface_det(1,1:local_intfr%counteryw) + increment

  increment = blk_b%move_amtz/factor
  local_intfr%pz_interface_det(1,1:local_intfr%counterzp) =&
          local_intfr%pz_interface_det(1,1:local_intfr%counterzp) + increment

  local_intfr%uz_interface_det(1,1:local_intfr%counterzu) =&
          local_intfr%uz_interface_det(1,1:local_intfr%counterzu) + increment

  local_intfr%vz_interface_det(1,1:local_intfr%counterzv) =&
          local_intfr%vz_interface_det(1,1:local_intfr%counterzv) + increment

  local_intfr%wz_interface_det(1,1:local_intfr%counterzw) =&
          local_intfr%wz_interface_det(1,1:local_intfr%counterzw) + increment

  call print_interface_detail("px", local_intfr%counterxp, local_intfr%px_interface_det, &
                              blk_a%xp, blk_b%xp)

  call print_interface_detail("py", local_intfr%counteryp, local_intfr%py_interface_det, &
                              blk_a%yp, blk_b%yp)

  call print_interface_detail("pz", local_intfr%counterzp, local_intfr%pz_interface_det, &
                              blk_a%zp, blk_b%zp)

  call print_interface_detail("ux", local_intfr%counterxu, local_intfr%ux_interface_det, &
                              blk_a%xu, blk_b%xu)

  call print_interface_detail("uy", local_intfr%counteryu, local_intfr%uy_interface_det, &
                              blk_a%yu, blk_b%yu)

  call print_interface_detail("uz", local_intfr%counterzu, local_intfr%uz_interface_det, &
                              blk_a%zu, blk_b%zu)

  call print_interface_detail("vx", local_intfr%counterxv, local_intfr%vx_interface_det, &
                              blk_a%xv, blk_b%xv)

  call print_interface_detail("vy", local_intfr%counteryv, local_intfr%vy_interface_det, &
                              blk_a%yv, blk_b%yv)

  call print_interface_detail("vz", local_intfr%counterzv, local_intfr%vz_interface_det, &
                              blk_a%zv, blk_b%zv)

  call print_interface_detail("wx", local_intfr%counterxw, local_intfr%wx_interface_det, &
                              blk_a%xw, blk_b%xw)

  call print_interface_detail("wy", local_intfr%counteryw, local_intfr%wy_interface_det, &
                              blk_a%yw, blk_b%yw)

  call print_interface_detail("wz", local_intfr%counterzw, local_intfr%wz_interface_det, &
                              blk_a%zw, blk_b%zw)

end subroutine change_block_interface

        SUBROUTINE cellCount_solid_coarse(blk)

        type(Block_t), intent(inout) :: blk
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
            END IF
         end do
         end do
         end do
         if (allocated(blk%fluidIndexPtr)) deallocate(blk%fluidIndexPtr)
         ALLOCATE (blk%fluidIndexPtr(blk%fluidCellCount,3))
         DO k = 2, blk%nz +1
         DO j = 2, blk%ny +1
         DO i = 2, blk%nx +1
            IF (blk%cell(i,j,k)==0) THEN
               iPt1 = iPt1 + 1
               blk%fluidIndexPtr(iPt1, 1) = i
               blk%fluidIndexPtr(iPt1, 2) = j
               blk%fluidIndexPtr(iPt1, 3) = k
            END IF
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
            END IF
         END DO

         if (allocated(blk%redCellIndexPtr)) deallocate(blk%redCellIndexPtr)
         if (allocated(blk%blackCellIndexPtr)) deallocate(blk%blackCellIndexPtr)

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
            END IF
         END DO

            print*, "1", blk%fluidCellCount, blk%redCellCount, blk%blackCellCount

        END SUBROUTINE cellCount_solid_coarse
end module biocfd_search

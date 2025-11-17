module biocfd_read_input
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, uc, u0, totime, &
       re, pi, omega4, &
       omega3, omega2, omega1, ita1, ita, &
       inor, freq, epsi, dxmin, dt_order, deltat, &
       blk_start, intfr
  use biocfd_interface_type, only: Interface_t
  use biocfd_block_type,only : Block_t
  implicit none

  private

  public :: readInput, readBlockInterface, readSurfaceMeshGmsh

  contains

      SUBROUTINE readInput(surGeoPoints,char_f,istart,itamax,pcItaMax,aoa,phase_angle,piv_pt)
       INTEGER (int64) :: i, g,io
       CHARACTER(len=160)  :: filename1
       INTEGER (int64),INTENT(OUT)   :: surGeoPoints, itamax,pcItaMax
       REAL(dp),intent(out) :: aoa,phase_angle,piv_pt
       CHARACTER (LEN = 3), INTENT(OUT)  :: char_f
       INTEGER,INTENT(OUT)               :: istart
       ! MB: Temporary variables added, to separate them out from type Blocks. Kept until
       !     not dependent on diff for checking code changes don't break code
       !     Variables removed from Blocks 'xstart, xend, ystart, yend, zstart, zend'
       REAL(dp),ALLOCATABLE,DIMENSION(:) :: xstart_temp,xend_temp,&
       ystart_temp,yend_temp,zstart_temp,zend_temp
       REAL(dp) :: alpha_m,theta_m,alpha_m1,theta_m1,mu_f,rho_f,l_c,u_tip,disp
       INTEGER (int64) :: intflines, nblocks
       NAMELIST /input_data/ nblocks, intflines,  &
                   itamax, epsi, pcItaMax,omega1,omega2,omega3,omega4, &
                   re,rho_f, mu_f, l_c, &
                   u0,  &
                   surGeoPoints, phase_angle, freq, aoa, piv_pt,alpha_m, theta_m, &
                   istart, dt_order,  inor, dxmin

  open(newunit=io, file="input_data.nml", status="old", action="read")
  read(io, NML=input_data)
  close(io)

        allocate(Block_t :: block(nblocks))
        allocate(Interface_t :: intfr(intflines))
        allocate(xstart_temp(size(block)),xend_temp(size(block)),&
        ystart_temp(size(block)),yend_temp(size(block)),&
        zstart_temp(size(block)),zend_temp(size(block)))
        OPEN(77, FILE = "body_search.dat", FORM = "formatted")
        DO g=1,size(block)
                READ(77,*) block(g)%i_startSearch, block(g)%i_endSearch, &
                           block(g)%j_startSearch, block(g)%j_endSearch, &
                           block(g)%k_startSearch, block(g)%k_endSearch
                print*, block(g)%i_startSearch, block(g)%i_endSearch, &
                           block(g)%j_startSearch, block(g)%j_endSearch, &
                           block(g)%k_startSearch, block(g)%k_endSearch
        END DO
        CLOSE(77)
        OPEN(77, FILE = "grid_shift.dat", FORM = "formatted")
        DO g=1,size(block)
        READ(77,*) block(g)%gx_shift, block(g)%gy_shift, block(g)%gz_shift
        block(g)%gx_shift=block(g)%gx_shift*0.001_dp
        block(g)%gy_shift=block(g)%gy_shift*0.001_dp
        block(g)%gz_shift=block(g)%gz_shift*0.001_dp
        END DO
        CLOSE(77)
        OPEN(77, FILE = "shift.dat", FORM = "formatted")

        DO g=1,size(block)
        READ(77,*) block(g)%xshift, block(g)%yshift, block(g)%zshift
        write(*,*) block(g)%xshift, block(g)%yshift, block(g)%zshift
        block(g)%xshift=block(g)%xshift*0.001_dp
        block(g)%yshift=block(g)%yshift*0.001_dp
        block(g)%zshift=block(g)%zshift*0.001_dp
        END DO
        CLOSE(77)
        if ( alpha_m /= 0 .and. theta_m /=0 ) then
                char_f = "bot"
        end if
        if ( alpha_m == 0 .and. theta_m /=0 ) then
                char_f = "amp"
        end if
        if ( alpha_m /= 0 .and. theta_m ==0 ) then
                char_f = "ang"
        end if
        blk_start=2
        alpha_m1=abs(alpha_m)
        theta_m1=abs(theta_m)
        dxmin = 0.001_dp*dxmin

        OPEN(77, FILE = "flap_amp.dat", FORM = "formatted")
        DO g=blk_start,size(block)
                READ(77,*) block(g)%a0
                write(*,*)g, block(g)%a0
        END DO
        CLOSE(77)
        mu_f = mu_f*1e-5_dp
        rho_f = rho_f
        l_c = 0.001_dp*l_c
        u0 = re*mu_f/(rho_f*l_c)
        uc = u0
        re = rho_f/mu_f
        freq = freq*u0/l_c
        deltat = 1._dp/(4._dp*freq*dt_order)  !0.00041666666666_dp! *5e-4
        disp = block(blk_start)%a0*cos(2*pi*freq*deltat)
       u_tip=4*((pi*45)/180)*l_c*freq

        Print*, "dxmin =", dxmin
        Print*, "u0 =", u0
        print*, "dt =",  deltat
        print*, "freq =", freq
        print*, "disp =", disp
        print*, "re =", re
        print*, "alpha_m =", alpha_m
        print*, "theta_m =", theta_m
        print*, "u_tip =", u_tip

        OPEN(77, FILE = "butter_move.dat", FORM = "formatted")
        DO g=blk_start,size(block)
                READ(77,*) block(g)%yamp, block(g)%bfreq
                block(g)%yamp=dxmin*block(g)%yamp
                write(*,*)g, block(g)%yamp, block(g)%bfreq
        END DO
        CLOSE(77)


        OPEN(111,FILE="init_params.dat",POSITION="APPEND",STATUS="unknown")
        WRITE(111,*) "dxmin, u0, deltat, re"
        WRITE(111,166) dxmin, u0, deltat, re
166     FORMAT(2F10.6,E15.6,F8.2)
        WRITE(111,*) "freq, disp, utip"
        WRITE(111,266) freq, disp, u_tip
266     FORMAT(3F10.6)
        WRITE(111,*) "alpha_m, theta_m"
        WRITE(111,366) alpha_m1, theta_m1
366     FORMAT(2F5.2)
        WRITE(111,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        CLOSE(111)

        ita = 0
        ita1 = 0
        totime = 0._dp
        print*, "dt =",  deltat, "ita = ", ita, "totime = ", totime

        OPEN(51, FILE = "block_details.dat", FORM = "formatted")
       DO i=1,size(block)

        read(51, *) xstart_temp(i),xend_temp(i),ystart_temp(i),&
                    yend_temp(i),zstart_temp(i),zend_temp(i),&
                    block(i)%nx, block(i)%ny, block(i)%nz, &
                    block(i)%dx, block(i)%dy, block(i)%dz
       END DO
       CLOSE(51)
        print*,"after allocation"

         DO i=1,size(block)

            ALLOCATE(block(i)%x1(block(i)%nx+3), &
                     block(i)%y1(block(i)%ny+3), &
                     block(i)%z1(block(i)%nz+3),&
                     block(i)%deltax(block(i)%nx+2), &
                     block(i)%deltay(block(i)%ny+2), &
                     block(i)%deltaz(block(i)%nz+2))
            ALLOCATE(block(i)%xu(block(i)%nx+3), &
                     block(i)%yu(block(i)%ny+2), &
                     block(i)%zu(block(i)%nz+2),&
                     block(i)%xv(block(i)%nx+2), &
                     block(i)%yv(block(i)%ny+3), &
                     block(i)%zv(block(i)%nz+2),&
                     block(i)%xw(block(i)%nx+2), &
                     block(i)%yw(block(i)%ny+2), &
                     block(i)%zw(block(i)%nz+3),&
                     block(i)%xp(block(i)%nx+2), &
                     block(i)%yp(block(i)%ny+2), &
                     block(i)%zp(block(i)%nz+2))

          block(i)%cintp=3
            END DO


         print*, "after allocation"

        do i=1, size(block)
            block(i)%fineg=block(i)%dx
            block(i)%dx=(1._dp/block(i)%dx)*0.001_dp
            block(i)%dy=(1._dp/block(i)%dy)*0.001_dp
            block(i)%dz=(1._dp/block(i)%dz)*0.001_dp

        end do

         DO i=1,size(block)
        print *,"For block blockno,xstart,xend,ystart,yend,nx,ny,dx,dy:",&
             i, xstart_temp(i)*0.001_dp,xend_temp(i)*0.001_dp,&
             ystart_temp(i)*0.001_dp,yend_temp(i)*0.001_dp,&
             zstart_temp(i)*0.001_dp,zend_temp(i)*0.001_dp,&
              block(i)%nx, block(i)%ny,block(i)%nz, block(i)%dx, block(i)%dy,block(i)%dz
         END DO


         Do g=1,size(block)

         WRITE(filename1,8282) g,block(g)%nx+1

 8282    FORMAT("xgrid_bk",i1,"_",i3.3,".txt")
         OPEN(61, FILE = filename1, FORM = "formatted")
         DO i = 2, block(g)%nx+2
            READ(61, *) block(g)%x1(i)
                block(g)%x1(i)=0.001_dp*block(g)%x1(i)
                    block(g)%x1(i)=block(g)%x1(i) + block(g)%gx_shift
         END DO
         CLOSE(61)

         DO i = 2, block(g)%nx+1
            block(g)%deltax(i) = block(g)%x1(i+1) - block(g)%x1(i)
         END DO

         block(g)%deltax(1) = block(g)%deltax(2)
         block(g)%deltax(block(g)%nx+2) = block(g)%deltax(block(g)%nx+1)
         block(g)%x1(1) = block(g)%x1(2) - block(g)%deltax(1)
         block(g)%x1(block(g)%nx+3) = block(g)%x1(block(g)%nx+2) + block(g)%deltax(block(g)%nx+2)

         WRITE(filename1,8383) g,block(g)%ny+1

 8383    FORMAT("ygrid_bk",i1,"_",i3.3,".txt")
         OPEN(62, FILE = filename1, FORM = "formatted")
         DO i = 2, block(g)%ny+2
            READ(62, *) block(g)%y1(i)
                 block(g)%y1(i)=0.001_dp*block(g)%y1(i)
                    block(g)%y1(i)=block(g)%y1(i) + block(g)%gy_shift
         END DO
         CLOSE(62)

         DO i = 2, block(g)%ny+1
            block(g)%deltay(i) = block(g)%y1(i+1) - block(g)%y1(i)
         END DO

         block(g)%deltay(1) = block(g)%deltay(2)
         block(g)%deltay(block(g)%ny+2) = block(g)%deltay(block(g)%ny+1)
         block(g)%y1(1) = block(g)%y1(2) - block(g)%deltay(1)
         block(g)%y1(block(g)%ny+3) = block(g)%y1(block(g)%ny+2) + block(g)%deltay(block(g)%ny+2)

         WRITE(filename1,8484) g,block(g)%nz+1
 8484    FORMAT("zgrid_bk",i1,"_",i3.3,".txt")
         OPEN(63, FILE =filename1, FORM = "formatted")
         DO i = 2, block(g)%nz+2
            READ(63, *) block(g)%z1(i)
                 block(g)%z1(i)=0.001_dp*block(g)%z1(i)
                    block(g)%z1(i)=block(g)%z1(i) + block(g)%gz_shift
         END DO
         CLOSE(63)

         DO i = 2, block(g)%nz+1
            block(g)%deltaz(i) = block(g)%z1(i+1) - block(g)%z1(i)
         END DO

         block(g)%deltaz(1) = block(g)%deltaz(2)
         block(g)%deltaz(block(g)%nz+2) = block(g)%deltaz(block(g)%nz+1)
         block(g)%z1(1) = block(g)%z1(2) - block(g)%deltaz(1)
         block(g)%z1(block(g)%nz+3) = block(g)%z1(block(g)%nz+2) + block(g)%deltaz(block(g)%nz+2)

        END DO



        DO g=1,size(block)
        DO i = 1, block(g)%nx+3
           block(g)%xu(i) = block(g)%x1(i)
        END DO
        END DO

        DO g=1,size(block)
        DO i = 1, block(g)%ny+3
           block(g)%yv(i) = block(g)%y1(i)
        END DO
        END DO

        DO g=1,size(block)
          DO i = 1, block(g)%nz+3
           block(g)%zw(i) = block(g)%z1(i)
          END DO
       END DO

        DO g=1,size(block)
           DO i = 1, block(g)%ny+2
             block(g)%yu(i) = 0.5_dp*(block(g)%y1(i)+block(g)%y1(i+1))
             block(g)%yw(i) = block(g)%yu(i)
             block(g)%yp(i) = block(g)%yu(i)
           END DO
        END DO

        DO g=1,size(block)
           DO i = 1, block(g)%nx+2
             block(g)%xv(i) = 0.5_dp*(block(g)%x1(i)+block(g)%x1(i+1))
             block(g)%xw(i) = block(g)%xv(i)
             block(g)%xp(i) = block(g)%xv(i)
            print*,"xp",g,i,block(g)%xp(i)
         END DO
        END DO

        DO g=1,size(block)
           DO i = 1, block(g)%nz+2
            block(g)%zu(i) = 0.5_dp*(block(g)%z1(i)+block(g)%z1(i+1))
            block(g)%zv(i) = block(g)%zu(i)
            block(g)%zp(i) = block(g)%zu(i)
          END DO
        END DO

      END SUBROUTINE readInput

      SUBROUTINE readSurfaceMeshGmsh(blk,surGeoPoints)
       type(Block_t), intent(inout) :: blk
       INTEGER(int64) :: n, i1, i2, i3, i5
       INTEGER (int64),INTENT(IN)   :: surGeoPoints

       OPEN(121, FILE ="geometries/butterflyMedium.msh", form = "formatted")               !READ SURFACE MESH FILE
        DO n = 1, 4
           READ (121,*)
        END DO
        READ (121,*) blk%ibNodes  !nsurf=total no. of points in file
        ALLOCATE(blk%ibNodeId(blk%ibNodes), blk%xnode(blk%ibNodes), &
                 blk%ynode(blk%ibNodes), blk%znode(blk%ibNodes))
        blk%ibNodeId = 50
        DO n = 1, blk%ibNodes
           READ (121,*) i1, blk%xnode(n), blk%ynode(n), blk%znode(n)
           blk%xnode(n)=blk%xnode(n)*0.001_dp
           blk%ynode(n)=blk%ynode(n)*0.001_dp
           blk%znode(n)=blk%znode(n)*0.001_dp
        END DO
        DO n = 1, 2
          READ (121,*)
        END DO
        READ (121,*) blk%ibElems   !no. of elements
        DO n = 1, surGeoPoints
          READ (121,*)
        END DO
        blk%ibElems = blk%ibElems-surGeoPoints
        ALLOCATE(blk%ibSurfId(blk%ibElems), blk%ibElP1(blk%ibElems), &
                 blk%ibElP2(blk%ibElems), blk%ibElP3(blk%ibElems))
        blk%ibElP1 = 0
        blk%ibElP2 = 0
        blk%ibElP3 = 0
        blk%ibSurfId = 0
        DO n = 1, blk%ibElems
        READ (121,*) i1, i2, i3, blk%ibSurfId(n), i5, &
        blk%ibElP1(n), blk%ibElP2(n), blk%ibElP3(n)
        END DO
       DO n = 1, blk%ibElems
       IF (blk%ibSurfId(n)==51) THEN
           blk%ibNodeId(blk%ibELP1(n)) = 51
           blk%ibNodeId(blk%ibELP2(n)) = 51
           blk%ibNodeId(blk%ibELP3(n)) = 51
           ELSE IF (blk%ibSurfId(n)==52) THEN
           blk%ibNodeId(blk%ibELP1(n)) = 52
           blk%ibNodeId(blk%ibELP2(n)) = 52
           blk%ibNodeId(blk%ibELP3(n)) = 52
          END IF
        END DO
        DO n = 1, blk%ibElems
        IF (blk%ibSurfId(n)==50) THEN
            blk%ibNodeId(blk%ibELP1(n)) = 50
            blk%ibNodeId(blk%ibELP2(n)) = 50
            blk%ibNodeId(blk%ibELP3(n)) = 50
          END IF
        END DO

       CLOSE(121)
       PRINT *, "SURFACE MESH READING COMPLETE"
       PRINT *, "ibNodes =", blk%ibNodes, "ibElems =", blk%ibElems

      END SUBROUTINE readSurfaceMeshGmsh

        SUBROUTINE readBlockInterface
        INTEGER(int64) :: i

         print*,"Inside readBlockInterface"
        OPEN(51, FILE = "interface_details.dat", FORM = "formatted")
        DO i=1, size(intfr)
        READ(51,*) intfr(i)%a_blk, intfr(i)%a_msh, intfr(i)%a_intf, &
                   intfr(i)%b_blk, intfr(i)%b_msh, intfr(i)%b_intf, &
                   intfr(i)%xintf_start, intfr(i)%xintf_end, &
                   intfr(i)%yintf_start, intfr(i)%yintf_end, &
                   intfr(i)%zintf_start, intfr(i)%zintf_end

        intfr(i)%xintf_start=intfr(i)%xintf_start *0.001_dp
        intfr(i)%xintf_end = intfr(i)%xintf_end   *0.001_dp
        intfr(i)%yintf_start=intfr(i)%yintf_start *0.001_dp
        intfr(i)%yintf_end = intfr(i)%yintf_end   *0.001_dp
        intfr(i)%zintf_start=intfr(i)%zintf_start *0.001_dp
        intfr(i)%zintf_end  =intfr(i)%zintf_end   *0.001_dp
        print*, intfr(i)%a_blk, intfr(i)%a_msh, intfr(i)%a_intf, &
                intfr(i)%b_blk, intfr(i)%b_msh, intfr(i)%b_intf, &
                intfr(i)%xintf_start, intfr(i)%xintf_end, &
                intfr(i)%yintf_start, intfr(i)%yintf_end, &
                intfr(i)%zintf_start, intfr(i)%zintf_end
        END DO
        CLOSE(51)

         end subroutine readBlockInterface
end module biocfd_read_input

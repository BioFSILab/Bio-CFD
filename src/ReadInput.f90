module biocfd_read_input
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  ! allow(use-all) - TODO: Aim to fix this in the future
  use global
  implicit none

  private

  public :: readInput, readBlockInterface, readSurfaceMeshGmsh

  contains

      SUBROUTINE readInput
       INTEGER (int64) :: i, j , k , g, nx_var, ny_var, nz_var
        CHARACTER(len=160)  :: filename1

        OPEN(60, FILE = 'inputdata', FORM = 'formatted')
        READ(60,*) nblocks, intflines,  &
                   itamax, epsi, pcItaMax,omega1,omega2,omega3,omega4, &
                   re,rho_f, mu_f, l_c, &
                   u0, v0,  w0,  &
                   surGeoPoints,a0y, phase_angle, freq, aoa, piv_pt,alpha_m, theta_m, &
                   istart, dt_order, inor, dxmin
        CLOSE(60)
        allocate(Blocks :: block(nblocks))
        allocate(Interfaces :: intfr(intflines))

        OPEN(77, FILE = 'body_search.dat', FORM = 'formatted')
        DO g=1,nblocks
                READ(77,*) block(g)%i_startSearch, block(g)%i_endSearch, &
                           block(g)%j_startSearch, block(g)%j_endSearch, &
                           block(g)%k_startSearch, block(g)%k_endSearch
                print*, block(g)%i_startSearch, block(g)%i_endSearch, &
                           block(g)%j_startSearch, block(g)%j_endSearch, &
                           block(g)%k_startSearch, block(g)%k_endSearch
        END DO
        CLOSE(77)
        OPEN(77, FILE = 'grid_shift.dat', FORM = 'formatted')
        DO g=1,nblocks
        READ(77,*) block(g)%gx_shift, block(g)%gy_shift, block(g)%gz_shift
        block(g)%gx_shift=block(g)%gx_shift*0.001_dp
        block(g)%gy_shift=block(g)%gy_shift*0.001_dp
        block(g)%gz_shift=block(g)%gz_shift*0.001_dp
        END DO
        CLOSE(77)
        OPEN(77, FILE = 'shift.dat', FORM = 'formatted')

        DO g=1,nblocks
        READ(77,*) block(g)%xshift, block(g)%yshift, block(g)%zshift
        write(*,*) block(g)%xshift, block(g)%yshift, block(g)%zshift
        block(g)%xshift=block(g)%xshift*0.001_dp
        block(g)%yshift=block(g)%yshift*0.001_dp
        block(g)%zshift=block(g)%zshift*0.001_dp
        END DO
        CLOSE(77)
        if ( alpha_m /= 0 .and. theta_m /=0 ) then
                char_f = 'bot'
        end if
        if ( alpha_m == 0 .and. theta_m /=0 ) then
                char_f = 'amp'
        end if
        if ( alpha_m /= 0 .and. theta_m ==0 ) then
                char_f = 'ang'
        end if
        pi = 4.D0*ATAN(1.D0)
        blk_start=2
        alpha_m1=abs(alpha_m)
        theta_m1=abs(theta_m)
        dxmin = 0.001_dp*dxmin

        OPEN(77, FILE = 'flap_amp.dat', FORM = 'formatted')
        DO g=blk_start,nblocks
                READ(77,*) block(g)%a0
                write(*,*)g, block(g)%a0
        END DO
        CLOSE(77)
        mu_f = mu_f*1e-5_dp
        rho_f = rho_f
        l_c = 0.001_dp*l_c
        rev = rho_f/(mu_f)
        u0 = re/(rev*l_c)
        u0 = u0*1
        uc=u0*1
        re = rev
        pi = 4.D0*ATAN(1.D0)
        a0y = 2*sin(pi/6)*l_c
        freq = freq*u0/l_c
        deltat = 1._dp/(4._dp*freq*dt_order)  !0.00041666666666_dp! *5e-4
        disp = block(blk_start)%a0*cos(2*pi*freq*deltat)
        alpha  = 1._dp
       u_tip=4*((pi*45)/180)*l_c*freq

        Print*, 'dxmin =', dxmin
        Print*, 'u0 =', u0
        print*, 'dt =',  deltat
        print*, 'freq =', freq
        print*, 'disp =', disp
        print*, 're =', re
        print*, 'alpha_m =', alpha_m
        print*, 'theta_m =', theta_m
        print*, 'u_tip =', u_tip

        OPEN(77, FILE = 'butter_move.dat', FORM = 'formatted')
        DO g=blk_start,nblocks
                READ(77,*) block(g)%yamp, block(g)%bfreq
                block(g)%yamp=dxmin*block(g)%yamp
                write(*,*)g, block(g)%yamp, block(g)%bfreq
        END DO
        CLOSE(77)


        OPEN(111,FILE='init_params.dat',ACCESS='Append',STATUS='unknown')
        WRITE(111,*) 'dxmin, u0, deltat, re'
        WRITE(111,166) dxmin, u0, deltat, re
166     FORMAT(2F10.6,E15.6,F8.2)
        WRITE(111,*) 'freq, disp, utip'
        WRITE(111,266) freq, disp, u_tip
266     FORMAT(3F10.6)
        WRITE(111,*) 'alpha_m, theta_m'
        WRITE(111,366) alpha_m1, theta_m1
366     FORMAT(2F5.2)
        WRITE(111,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        CLOSE(111)

        ita = 0
        ita1 = 0
        totime = 0._dp
        print*, 'dt =',  deltat, 'ita = ', ita, 'totime = ', totime

        OPEN(51, FILE = 'block_details.dat', FORM = 'formatted')
       DO i=1,nblocks

        read(51, *) block(i)%xstart, block(i)%xend, &
                    block(i)%ystart, block(i)%yend, &
                    block(i)%zstart, block(i)%zend, &
                    block(i)%nx, block(i)%ny, block(i)%nz, &
                    block(i)%dx, block(i)%dy, block(i)%dz
        block(i)%xstart= block(i)%xstart* 0.001_dp
        block(i)%xend=  block(i)%xend*0.001_dp
        block(i)%ystart= block(i)%ystart*0.001_dp
        block(i)%yend=block(i)%yend*0.001_dp
        block(i)%zstart=block(i)%zstart*0.001_dp
        block(i)%zend=block(i)%zend*0.001_dp
       END DO
       CLOSE(51)
        print*,'after allocation'

         DO i=1,nblocks
             nx_var=block(i)%nx
             ny_var=block(i)%ny
             nz_var=block(i)%nz

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
                     block(i)%zp(block(i)%nz+2),&
                     block(i)%xp_dum(block(i)%nx+2), &
                     block(i)%yp_dum(block(i)%ny+2), &
                     block(i)%zp_dum(block(i)%nz+2))

          ALLOCATE ( block(i)%xp1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%yp1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%zp1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%xpn1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%ypn1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%zpn1(nx_var+2,ny_var+2,nz_var+2) )
          block(i)%cintp=3
            END DO


         print*, 'after allocation'

        do i=1, nblocks
            block(i)%fineg=block(i)%dx
            block(i)%dx=(1._dp/block(i)%dx)*0.001_dp
            block(i)%dy=(1._dp/block(i)%dy)*0.001_dp
            block(i)%dz=(1._dp/block(i)%dz)*0.001_dp

        end do

         DO i=1,nblocks
        print *,'For block blockno,xstart,xend,ystart,yend,nx,ny,dx,dy:',&
                   i, block(i)%xstart, block(i)%xend, &
                   block(i)%ystart, block(i)%yend, &
                   block(i)%zstart, block(i)%zend, &
                   block(i)%nx, block(i)%ny,block(i)%nz, block(i)%dx, block(i)%dy,block(i)%dz
         END DO


         Do g=1,nblocks

         WRITE(filename1,8282) g,block(g)%nx+1

 8282    FORMAT('xgrid_bk',i1,'_',i3.3,'.txt')
         OPEN(61, FILE = filename1, FORM = 'formatted')
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

 8383    FORMAT('ygrid_bk',i1,'_',i3.3,'.txt')
         OPEN(62, FILE = filename1, FORM = 'formatted')
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
 8484    FORMAT('zgrid_bk',i1,'_',i3.3,'.txt')
         OPEN(63, FILE =filename1, FORM = 'formatted')
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



        DO g=1,nblocks
        DO i = 1, block(g)%nx+3
           block(g)%xu(i) = block(g)%x1(i)
        ENDDO
        ENDDO

        DO g=1,nblocks
        DO i = 1, block(g)%ny+3
           block(g)%yv(i) = block(g)%y1(i)
        ENDDO
        ENDDO

        DO g=1,nblocks
          DO i = 1, block(g)%nz+3
           block(g)%zw(i) = block(g)%z1(i)
          ENDDO
       ENDDO

        DO g=1,nblocks
           DO i = 1, block(g)%ny+2
             block(g)%yu(i) = 0.5_dp*(block(g)%y1(i)+block(g)%y1(i+1))
             block(g)%yw(i) = block(g)%yu(i)
             block(g)%yp(i) = block(g)%yu(i)
           END DO
        ENDDO

        DO g=1,nblocks
           DO i = 1, block(g)%nx+2
             block(g)%xv(i) = 0.5_dp*(block(g)%x1(i)+block(g)%x1(i+1))
             block(g)%xw(i) = block(g)%xv(i)
             block(g)%xp(i) = block(g)%xv(i)
            print*,'xp',g,i,block(g)%xp(i)
         END DO
        ENDDO

        DO g=1,nblocks
           DO i = 1, block(g)%nz+2
            block(g)%zu(i) = 0.5_dp*(block(g)%z1(i)+block(g)%z1(i+1))
            block(g)%zv(i) = block(g)%zu(i)
            block(g)%zp(i) = block(g)%zu(i)
          END DO
        ENDDO
        DO g=1,nblocks

        DO k=1,block(g)%nz+1
        DO j=1,block(g)%ny+1
        DO i=1,block(g)%nx+1

        block(g)%xpn1(i,j,k)=block(g)%x1(i)
        block(g)%ypn1(i,j,k)=block(g)%y1(j)
        block(g)%zpn1(i,j,k)=block(g)%z1(k)


        END DO
        END DO
        END DO
        END DO

        DO g=1,nblocks
        DO k=2,block(g)%nz+1
        DO j=2,block(g)%ny+1
        DO i=2,block(g)%nx+1

        block(g)%xp1(i,j,k)=block(g)%xp(i)
        block(g)%yp1(i,j,k)=block(g)%yp(j)
        block(g)%zp1(i,j,k)=block(g)%zp(k)


        END DO
        END DO
        END DO
        END DO

      END SUBROUTINE readInput

      SUBROUTINE readSurfaceMeshGmsh
       INTEGER(int64) :: n, i1, i2, i3, i5, g
       CHARACTER (LEN = 72) :: cLine


       DO g=blk_start, nblocks
       OPEN(121, FILE ='geometries/butterflyMedium.msh', form = 'formatted')               !READ SURFACE MESH FILE
        DO n = 1, 4
           READ (121,*) cLine
        END DO
        READ (121,*) block(g)%ibNodes  !nsurf=total no. of points in file
        ALLOCATE(block(g)%ibNodeId(block(g)%ibNodes), block(g)%xnode(block(g)%ibNodes), &
                 block(g)%ynode(block(g)%ibNodes), block(g)%znode(block(g)%ibNodes))
        block(g)%ibNodeId = 50
        DO n = 1, block(g)%ibNodes
           READ (121,*) i1, block(g)%xnode(n), block(g)%ynode(n), block(g)%znode(n)
           block(g)%xnode(n)=block(g)%xnode(n)*0.001_dp
           block(g)%ynode(n)=block(g)%ynode(n)*0.001_dp
           block(g)%znode(n)=block(g)%znode(n)*0.001_dp
        END DO
        DO n = 1, 2
          READ (121,*) line
        END DO
        READ (121,*) block(g)%ibElems   !no. of elements
        DO n = 1, surGeoPoints
          READ (121,*) cLine
        END DO
        block(g)%ibElems = block(g)%ibElems-surGeoPoints
        ALLOCATE(block(g)%ibSurfId(block(g)%ibElems), block(g)%ibElP1(block(g)%ibElems), &
                 block(g)%ibElP2(block(g)%ibElems), block(g)%ibElP3(block(g)%ibElems))
        block(g)%ibElP1 = 0
        block(g)%ibElP2 = 0
        block(g)%ibElP3 = 0
        block(g)%ibSurfId = 0
        DO n = 1, block(g)%ibElems
        READ (121,*) i1, i2, i3, block(g)%ibSurfId(n), i5, &
        block(g)%ibElP1(n), block(g)%ibElP2(n), block(g)%ibElP3(n)
        END DO
       DO n = 1, block(g)%ibElems
       IF (block(g)%ibSurfId(n)==51) THEN
           block(g)%ibNodeId(block(g)%ibELP1(n)) = 51
           block(g)%ibNodeId(block(g)%ibELP2(n)) = 51
           block(g)%ibNodeId(block(g)%ibELP3(n)) = 51
           ELSEIF (block(g)%ibSurfId(n)==52) THEN
           block(g)%ibNodeId(block(g)%ibELP1(n)) = 52
           block(g)%ibNodeId(block(g)%ibELP2(n)) = 52
           block(g)%ibNodeId(block(g)%ibELP3(n)) = 52
          ENDIF
        ENDDO
        DO n = 1, block(g)%ibElems
        IF (block(g)%ibSurfId(n)==50) THEN
            block(g)%ibNodeId(block(g)%ibELP1(n)) = 50
            block(g)%ibNodeId(block(g)%ibELP2(n)) = 50
            block(g)%ibNodeId(block(g)%ibELP3(n)) = 50
          ENDIF
        ENDDO

       CLOSE(121)
       PRINT *, 'SURFACE MESH READING COMPLETE'
       PRINT *, 'ibNodes =', block(g)%ibNodes, 'ibElems =', block(g)%ibElems
       END DO

      END SUBROUTINE readSurfaceMeshGmsh

        SUBROUTINE readBlockInterface
        INTEGER(int64) :: i

         print*,'Inside readBlockInterface'
        OPEN(51, FILE = 'interface_details.dat', FORM = 'formatted')
        DO i=1, intflines
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

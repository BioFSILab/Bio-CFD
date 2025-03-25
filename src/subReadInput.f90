module biocfd_read_input
  use global
  implicit none

  contains
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE readInput
       USE global
       IMPLICIT NONE
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER (KIND=8) :: i, j , k , g, r, nx_var, ny_var, nz_var
        CHARACTER(len=160)  :: filename1

        !OPEN(60, FILE = 'geometries/2blk/inputdata', FORM = 'formatted')
        OPEN(60, FILE = 'inputdata', FORM = 'formatted')
       !READ(60,*) nx, ny, nz,  &
       !           re,      &
       !           itamax, eps1, pcItaMax, &
       !           u0, v0,  w0, &
       !           surGeoPoints, a0, freq, aoa, piv_pt,&
       !           xShift, yShift, zShift,           &
       !           istart, dt_order, inor, &
       !           i_startSearch, i_endSearch, &
       !           j_startSearch, j_endSearch, &
       !           k_startSearch, k_endSearch
        READ(60,*) nblocks, intflines,  &
                   !re,      &
                   itamax, epsi, pcItaMax,omega1,omega2,omega3,omega4, &
                   re,rho_f, mu_f, l_c, &
                   u0, v0,  w0,  &
                   !surGeoPoints,a0, a0y, phase_angle, freq, aoa, piv_pt,alpha_m, theta_m, &
                   surGeoPoints,a0y, phase_angle, freq, aoa, piv_pt,alpha_m, theta_m, &
                   istart, dt_order, inor, dxmin
                  !i_startSearch, i_endSearch, &
                  !j_startSearch, j_endSearch, &
                  !k_startSearch, k_endSearch
        CLOSE(60)
        allocate(Blocks :: block(nblocks))
        allocate(Interfaces :: intfr(intflines))

        OPEN(77, FILE = 'body_search.dat', FORM = 'formatted')
        !OPEN(77, FILE = 'geometries/2blk/0.01/body_search.dat', FORM = 'formatted')
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
        !OPEN(77, FILE = 'geometries/2blk/0.01/body_search.dat', FORM = 'formatted')
        DO g=1,nblocks
        READ(77,*) block(g)%gx_shift, block(g)%gy_shift, block(g)%gz_shift
        block(g)%gx_shift=block(g)%gx_shift*0.001
        block(g)%gy_shift=block(g)%gy_shift*0.001
        block(g)%gz_shift=block(g)%gz_shift*0.001
        END DO
        CLOSE(77)
        OPEN(77, FILE = 'shift.dat', FORM = 'formatted')
        !OPEN(77, FILE = 'geometries/2blk/0.01/body_search.dat', FORM = 'formatted')
        DO g=1,nblocks
        READ(77,*) block(g)%xshift, block(g)%yshift, block(g)%zshift
        write(*,*) block(g)%xshift, block(g)%yshift, block(g)%zshift
        block(g)%xshift=block(g)%xshift*0.001
        block(g)%yshift=block(g)%yshift*0.001
        block(g)%zshift=block(g)%zshift*0.001
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
        dxmin = 0.001*dxmin
        !xShift = 0.001*xShift
        !yShift = 0.001*yShift
        !zShift = 0.001*zShift
        OPEN(77, FILE = 'flap_amp.dat', FORM = 'formatted')
        DO g=blk_start,nblocks
                READ(77,*) block(g)%a0
                write(*,*)g, block(g)%a0
        END DO
        CLOSE(77)
        mu_f = mu_f*1e-5
        rho_f = rho_f
        l_c = 0.001*l_c
        rev = rho_f/(mu_f)
        u0 = re/(rev*l_c)
        u0 = u0*1
        uc=u0*1
        re = rev
        pi = 4.D0*ATAN(1.D0)
        a0y = 2*sin(pi/6)*l_c
        freq = freq*u0/l_c
        deltat = 1./(4.*freq*dt_order)  !0.00041666666666_rk! *5e-4
        disp = block(blk_start)%a0*cos(2*pi*freq*deltat)
        alpha  = 1._rk
      !!freq=(mu_f*re*180)/(4*theta_m*pi*l_c*l_c)
      !!u_tip=4*((pi*theta_m)/180)*l_c*freq
      ! freq=(mu_f*re*180)/(4*45*pi*l_c*l_c)
       u_tip=4*((pi*45)/180)*l_c*freq
      !rev = rho_f/(mu_f)
      !!u0 = re/(rev*l_c)
      !u0 = 1.
      !!re = rev
      !!a0y = 2*sin(pi/6)*l_c
      !!freq = freq*u0/l_c
      !!deltat = 1./(4.*freq*dt_order)!0.00041666666666_rk! *5e-4
      !!lwing=74.
      !!lwing = 0.001*lwing

      !!deltat=(dxmin)/(u_tip*dt_order)
      !!disp = a0*cos(2*pi*freq*deltat)
      !!alpha  = 1._rk
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

        !uc = 1.
        ita = 0
        ita1 = 0
        totime = 0._rk
        print*, 'dt =',  deltat, 'ita = ', ita, 'totime = ', totime
        !print*, 'total cells =', nx*ny*nz

        !OPEN(51, FILE = 'geometries/2blk/0.01/block_details.dat', FORM = 'formatted')
        OPEN(51, FILE = 'block_details.dat', FORM = 'formatted')
       DO i=1,nblocks

        read(51, *) block(i)%xstart, block(i)%xend, block(i)%ystart, block(i)%yend,block(i)%zstart,block(i)%zend, block(i)%nx, block(i)%ny, block(i)%nz, block(i)%dx, block(i)%dy, block(i)%dz
        block(i)%xstart= block(i)%xstart* 0.001
        block(i)%xend=  block(i)%xend*0.001
        block(i)%ystart= block(i)%ystart*0.001
        block(i)%yend=block(i)%yend*0.001
        block(i)%zstart=block(i)%zstart*0.001
        block(i)%zend=block(i)%zend*0.001
       END DO
       CLOSE(51)
        print*,'after allocation'
        !call findSolidBlk

         DO i=1,nblocks
             nx_var=block(i)%nx
             ny_var=block(i)%ny
             nz_var=block(i)%nz

            ALLOCATE( block(i)%x1(block(i)%nx+3),block(i)%y1(block(i)%ny+3), block(i)%z1(block(i)%nz+3),&
                      block(i)%deltax(block(i)%nx+2), block(i)%deltay(block(i)%ny+2), block(i)%deltaz(block(i)%nz+2))
            ALLOCATE( block(i)%xu(block(i)%nx+3), block(i)%yu(block(i)%ny+2),block(i)%zu(block(i)%nz+2),&
                      block(i)%xv(block(i)%nx+2), block(i)%yv(block(i)%ny+3),block(i)%zv(block(i)%nz+2),&
                      block(i)%xw(block(i)%nx+2), block(i)%yw(block(i)%ny+2),block(i)%zw(block(i)%nz+3),&
                      block(i)%xp(block(i)%nx+2), block(i)%yp(block(i)%ny+2),block(i)%zp(block(i)%nz+2),&
                      block(i)%xp_dum(block(i)%nx+2), block(i)%yp_dum(block(i)%ny+2),block(i)%zp_dum(block(i)%nz+2))

          ALLOCATE ( block(i)%xp1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%yp1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%zp1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%xpn1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%ypn1(nx_var+2,ny_var+2,nz_var+2) )
          ALLOCATE ( block(i)%zpn1(nx_var+2,ny_var+2,nz_var+2) )
          block(i)%cintp=3
            END DO


         print*, 'after allocation'
      ! ALLOCATE( x1(nx+3), y1(ny+3), z1(nz+3), &
      !           deltax(nx+2), deltay(ny+2), deltz(nz+2))
      ! ALLOCATE( xu(nx+3), yu(ny+2), zu(nz+2), &
      !           xv(nx+2), yv(ny+3), zv(nz+2), &
      !           xw(nx+2), yw(ny+2), zw(nz+3), &
      !           xp(nx+2), yp(ny+2), zp(nz+2))

        do i=1, nblocks
            block(i)%fineg=block(i)%dx
            block(i)%dx=(1./block(i)%dx)*0.001
            block(i)%dy=(1./block(i)%dy)*0.001
            block(i)%dz=(1./block(i)%dz)*0.001

        end do

         DO i=1,nblocks
        print *,'For block blockno,xstart,xend,ystart,yend,nx,ny,dx,dy:',&
                   i,block(i)%xstart,block(i)%xend, block(i)%ystart,block(i)%yend,block(i)%zstart,block(i)%zend, &
                   block(i)%nx, block(i)%ny,block(i)%nz, block(i)%dx, block(i)%dy,block(i)%dz
         END DO

!       ! DO r = 1, nblocks-1
!       DO r = 1, nblocks
!              !j=fl_blk(r)
!              j=r
!           !DO j =1, nblocks
!           DO i = 2, block(j)%nx+1
!              block(j)%deltax(i) = block(j)%dx
!              block(j)%x1(i)=(i-2) *block(j)%dx + block(j)%xstart
!           END DO
!           block(j)%x1(block(j)%nx+2) =block(j)%x1(block(j)%nx+1) + block(j)%deltax(block(j)%nx+1)
!
!           block(j)%deltax(1)    = block(j)%deltax(2)
!           block(j)%deltax(block(j)%nx+2) =block(j)%deltax(block(j)%nx+1)
!           block(j)%x1(1)        = block(j)%x1(2) - block(j)%deltax(1)
!           block(j)%x1(block(j)%nx+3)     = block(j)%x1(block(j)%nx+2)+ block(j)%deltax(block(j)%nx+2)
!      end do


!       !DO r = 1, nblocks-1
!       DO r = 1, nblocks
!          !j=fl_blk(r)
!          j=r
!       !DO j =1, nblocks
!        DO i = 2, block(j)%ny+1
!          block(j)%deltay(i) = block(j)%dy
!          block(j)%y1(i)=(i-2) *block(j)%dy + block(j)%ystart
!       END DO
!       block(j)%y1(block(j)%ny+2) =block(j)%y1(block(j)%ny+1) + block(j)%deltay(block(j)%ny+1)

!       block(j)%deltay(1)    = block(j)%deltay(2)
!       block(j)%deltay(block(j)%ny+2) =block(j)%deltay(block(j)%ny+1)
!       block(j)%y1(1)        = block(j)%y1(2) - block(j)%deltay(1)
!       block(j)%y1(block(j)%ny+3)     = block(j)%y1(block(j)%ny+2)+ block(j)%deltay(block(j)%ny+2)
!     end do

!       !DO r = 1, nblocks-1
!       DO r = 1, nblocks
!              !j=fl_blk(r)
!              j=r
!           !DO j =1, nblocks
!            DO i = 2, block(j)%nz+1
!              block(j)%deltaz(i) = block(j)%dz
!              block(j)%z1(i)=(i-2) *block(j)%dz + block(j)%zstart
!           END DO
!           block(j)%z1(block(j)%nz+2) =block(j)%z1(block(j)%nz+1) +block(j)%deltaz(block(j)%nz+1)
!
!           block(j)%deltaz(1)    = block(j)%deltaz(2)
!           block(j)%deltaz(block(j)%nz+2) =block(j)%deltaz(block(j)%nz+1)
!           block(j)%z1(1)        = block(j)%z1(2) - block(j)%deltaz(1)
!           block(j)%z1(block(j)%nz+3)     = block(j)%z1(block(j)%nz+2)+ block(j)%deltaz(block(j)%nz+2)
!      end do



         Do g=1,nblocks
         !OPEN(61, FILE = 'xgrid_d50_498_3d_30.txt', FORM = 'formatted')
         !OPEN(61, FILE = 'xgrid_257_d20_3d_30.txt', FORM = 'formatted')

         WRITE(filename1,8282) g,block(g)%nx+1
 !8282    FORMAT('geometries/2blk/0.01/xgrid_bk',i1,'_',i3.3,'.txt')
 8282    FORMAT('xgrid_bk',i1,'_',i3.3,'.txt')
         OPEN(61, FILE = filename1, FORM = 'formatted')
         DO i = 2, block(g)%nx+2
            READ(61, *) block(g)%x1(i)
                block(g)%x1(i)=0.001*block(g)%x1(i)
                    block(g)%x1(i)=block(g)%x1(i) + block(g)%gx_shift
                !print*,'x1',g,block(g)%x1(i)
         END DO
         CLOSE(61)

         DO i = 2, block(g)%nx+1
            block(g)%deltax(i) = block(g)%x1(i+1) - block(g)%x1(i)
         END DO

         block(g)%deltax(1)    = block(g)%deltax(2)
         block(g)%deltax(block(g)%nx+2) = block(g)%deltax(block(g)%nx+1)
         block(g)%x1(1)        = block(g)%x1(2) - block(g)%deltax(1)
         block(g)%x1(block(g)%nx+3)     = block(g)%x1(block(g)%nx+2) + block(g)%deltax(block(g)%nx+2)

         !OPEN(62, FILE = 'ygrid_111_d20_3d_30.txt', FORM = 'formatted')
         WRITE(filename1,8383) g,block(g)%ny+1
! 8383    FORMAT('geometries/2blk/0.01/ygrid_bk',i1,'_',i3.3,'.txt')
 8383    FORMAT('ygrid_bk',i1,'_',i3.3,'.txt')
         OPEN(62, FILE = filename1, FORM = 'formatted')
         DO i = 2, block(g)%ny+2
            READ(62, *) block(g)%y1(i)
               ! print*,'y1',g,block(g)%y1(i)
                 block(g)%y1(i)=0.001*block(g)%y1(i)
                    block(g)%y1(i)=block(g)%y1(i) + block(g)%gy_shift
         END DO
         CLOSE(62)

         DO i = 2, block(g)%ny+1
            block(g)%deltay(i) = block(g)%y1(i+1) - block(g)%y1(i)
         END DO

         block(g)%deltay(1)    = block(g)%deltay(2)
         block(g)%deltay(block(g)%ny+2) = block(g)%deltay(block(g)%ny+1)
         block(g)%y1(1)        = block(g)%y1(2) - block(g)%deltay(1)
         block(g)%y1(block(g)%ny+3)     = block(g)%y1(block(g)%ny+2) + block(g)%deltay(block(g)%ny+2)


         !OPEN(63, FILE = 'zgrid_41_d20_3d_2.txt', FORM = 'formatted')
         WRITE(filename1,8484) g,block(g)%nz+1
 !8484    FORMAT('geometries/2blk/0.01/zgrid_bk',i1,'_',i3.3,'.txt')
 8484    FORMAT('zgrid_bk',i1,'_',i3.3,'.txt')
         OPEN(63, FILE =filename1, FORM = 'formatted')
         DO i = 2, block(g)%nz+2
            READ(63, *) block(g)%z1(i)
               ! print*,'z1',g,block(g)%z1(i)
                 block(g)%z1(i)=0.001*block(g)%z1(i)
                    block(g)%z1(i)=block(g)%z1(i) + block(g)%gz_shift
         END DO
         CLOSE(63)

         DO i = 2, block(g)%nz+1
            block(g)%deltaz(i) = block(g)%z1(i+1) - block(g)%z1(i)
         END DO

         block(g)%deltaz(1)    = block(g)%deltaz(2)
         block(g)%deltaz(block(g)%nz+2) = block(g)%deltaz(block(g)%nz+1)
         block(g)%z1(1)        = block(g)%z1(2) - block(g)%deltaz(1)
         block(g)%z1(block(g)%nz+3)     = block(g)%z1(block(g)%nz+2) + block(g)%deltaz(block(g)%nz+2)

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
           block(g)%yu(i) = 0.5_rk*(block(g)%y1(i)+block(g)%y1(i+1))
	   block(g)%yw(i) = block(g)%yu(i)
           block(g)%yp(i) = block(g)%yu(i)
        END DO
        ENDDO

        DO g=1,nblocks
        DO i = 1, block(g)%nx+2
           block(g)%xv(i) = 0.5_rk*(block(g)%x1(i)+block(g)%x1(i+1))
           block(g)%xw(i) = block(g)%xv(i)
           block(g)%xp(i) = block(g)%xv(i)
            print*,'xp',g,i,block(g)%xp(i)
        END DO
        ENDDO

        DO g=1,nblocks
	 DO i = 1, block(g)%nz+2
           block(g)%zu(i) = 0.5_rk*(block(g)%z1(i)+block(g)%z1(i+1))
	   block(g)%zv(i) = block(g)%zu(i)
           block(g)%zp(i) = block(g)%zu(i)
        END DO
        ENDDO
        DO g=1,nblocks
       !DO k=2,block(g)%nz+1
       !DO j=2,block(g)%ny+1
       !DO i=2,block(g)%nx+1
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

       !!$acc update device (deltax, deltay, deltaz
       !!$acc update device (x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp)
      END SUBROUTINE readInput

!************************************************************************************************

      SUBROUTINE readSurfaceMeshGmsh
       USE global
       IMPLICIT NONE
       INTEGER (KIND = 8) :: n, i1, i2, i3, i4, i5, i6, i7, gPoints, g
       CHARACTER (LEN = 72) :: cLine


       DO g=blk_start, nblocks
       !OPEN(121, FILE ='geometries/sphere_0.0075r_0.00025.msh', form = 'formatted')               !READ SURFACE MESH FILE
       OPEN(121, FILE ='geometries/butterflyMedium.msh', form = 'formatted')               !READ SURFACE MESH FILE
       !OPEN(121, FILE ='geometries/sphere_0.75r_0.025.msh', form = 'formatted')               !READ SURFACE MESH FILE
       !OPEN(121, FILE ='geometries/sphere_0.75r_0.03.msh', form = 'formatted')               !READ SURFACE MESH FILE
       !OPEN(121, FILE ='geometries/sphere_0.75r_0.01.msh', form = 'formatted')               !READ SURFACE MESH FILE
       !OPEN(121, FILE ='geometries/sphere_0.083r_0.004.msh', form = 'formatted')               !READ SURFACE MESH FILE
        DO n = 1, 4
           READ (121,*) cLine
        END DO
        READ (121,*) block(g)%ibNodes  !nsurf=total no. of points in file
        ALLOCATE ( block(g)%ibNodeId(block(g)%ibNodes), block(g)%xnode(block(g)%ibNodes), block(g)%ynode(block(g)%ibNodes), block(g)%znode(block(g)%ibNodes) )
        block(g)%ibNodeId = 50
        DO n = 1, block(g)%ibNodes
           READ (121,*) i1, block(g)%xnode(n), block(g)%ynode(n), block(g)%znode(n)
           !READ (121,*) i1, ynode(n), znode(n), xnode(n)
           block(g)%xnode(n)=block(g)%xnode(n)*0.001
           block(g)%ynode(n)=block(g)%ynode(n)*0.001
           block(g)%znode(n)=block(g)%znode(n)*0.001
        END DO
        DO n = 1, 2
          READ (121,*) line
        END DO
        READ (121,*) block(g)%ibElems   !no. of elements
        DO n = 1, surGeoPoints
          READ (121,*) cLine
        END DO
        block(g)%ibElems = block(g)%ibElems-surGeoPoints
            ALLOCATE ( block(g)%ibSurfId(block(g)%ibElems), block(g)%ibElP1(block(g)%ibElems), block(g)%ibElP2(block(g)%ibElems), block(g)%ibElP3(block(g)%ibElems) )
        block(g)%ibElP1 = 0
        block(g)%ibElP2 = 0
        block(g)%ibElP3 = 0
        block(g)%ibSurfId = 0
        DO n = 1, block(g)%ibElems
        READ (121,*) i1, i2, i3, block(g)%ibSurfId(n), i5, block(g)%ibElP1(n), block(g)%ibElP2(n), block(g)%ibElP3(n)
          !Print*, n, ibElP1(n), ibElP2(n)
        END DO
       !do n = 1, ibNodes
         !print*, n, xnode(n), ynode(n)
       !end do
       !  if (a4 .ne. 50)then
       !        num=num+1
       !  i1         =  a1
       !  i2         =  a2
       !  i3         =  a3
       !  ibSurfId(num)= a4
       !  i5         = a5
       !  ibElP1(num)  = a6
       !  ibElP2(num)  =  a7
       !  ibElP3(num)  = a8
       !  endif
          !Print*, n, ibsurfId(n), ibElP2(n)
       !do n = 1, ibNodes
         !print*, n, xnode(n), ynode(n)
       !end do
       ! ibElems= num
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

       !DO n = 1, ibNodes
       !   WRITE(*,*) ibNodeId(n)
       !ENDDO
       CLOSE(121)
       PRINT *, 'SURFACE MESH READING COMPLETE'
       PRINT *, 'ibNodes =', block(g)%ibNodes, 'ibElems =', block(g)%ibElems
       END DO

       !!$acc update device (xnode, ynode, znode, ibElP1, ibElP2, ibElP3)
      END SUBROUTINE readSurfaceMeshGmsh
!***********************************************************************

        SUBROUTINE readBlockInterface
        use global
        INTEGER (KIND=8) :: i, j, k, g

       !ALLOCATE(a_blk(intfLines),a_msh(intfLines),a_intf(intfLines),b_blk(intfLines),b_msh(intfLines),b_intf(intfLines),xintf_start(intfLines),xintf_end(intfLines),yintf_start(intfLines),yintf_end(i             ntfLines))

      !type(Interfaces),dimension(intflines) :: intfr
         print*,'Inside readBlockInterface'
        !OPEN(51, FILE = 'geometries/2blk/0.01/interface_details.dat', FORM = 'formatted')
        OPEN(51, FILE = 'interface_details.dat', FORM = 'formatted')
        DO i=1, intflines
        READ(51,*)intfr(i)%a_blk,intfr(i)%a_msh,intfr(i)%a_intf,intfr(i)%b_blk,intfr(i)%b_msh,intfr(i)%b_intf,intfr(i)%xintf_start,intfr(i)%xintf_end,intfr(i)%yintf_start,intfr(i)%yintf_end,intfr(i)%zintf_start,intfr(i)%zintf_end

        intfr(i)%xintf_start=intfr(i)%xintf_start *0.001
        intfr(i)%xintf_end = intfr(i)%xintf_end   *0.001
        intfr(i)%yintf_start=intfr(i)%yintf_start *0.001
        intfr(i)%yintf_end = intfr(i)%yintf_end   *0.001
        intfr(i)%zintf_start=intfr(i)%zintf_start *0.001
        intfr(i)%zintf_end  =intfr(i)%zintf_end   *0.001
        print*,intfr(i)%a_blk,intfr(i)%a_msh,intfr(i)%a_intf,intfr(i)%b_blk,intfr(i)%b_msh,intfr(i)%b_intf,intfr(i)%xintf_start,intfr(i)%xintf_end,intfr(i)%yintf_start,intfr(i)%yintf_end,intfr(i)%zintf_start, intfr(i)%zintf_end
        END DO
        CLOSE(51)

        !print*,a_blk(1),a_msh(1),a_intf(1),b_blk(1),b_msh(1),b_intf(1),xintf_start(1),xintf_end(1),yintf_start(1),yintf_end(1)
         end subroutine readBlockInterface



!***********************************************************************
!     SUBROUTINE readSurfaceMeshGambit
!      USE global
!      IMPLICIT NONE
!      INTEGER (KIND = 8) :: n, i1, i2, i3, i4
!      CHARACTER (LEN = 72) :: cLine
!
!      OPEN(121, FILE = 'geometries/bend_ellipse_ha.neu', form = 'formatted')               !READ SURFACE MESH FILE
!       DO n = 1, 6
!          READ (121,*) cLine
!       END DO
!       READ (121,*) ibNodes, ibElems, i1, i2, i3, i4
!       ALLOCATE ( xnode(ibNodes), ynode(ibNodes), znode(ibNodes) )
!       DO n = 1, 2
!          READ (121,*) cLine
!       END DO
!       DO n = 1, ibNodes
!          READ (121,*) i1, xnode(n), ynode(n), znode(n)
!       END DO
!       DO n = 1, 2
!         READ (121,*) cLine
!       END DO
!       ALLOCATE ( ibElP1(ibElems), ibElP2(ibElems), ibElP3(ibElems) )
!       ibElP1 = 0
!       ibElP2 = 0
!       ibElP3 = 0
!       DO n = 1, ibElems
!         READ (121,*) i1, i2, i3, ibElP1(n), ibElP2(n), ibElP3(n)
!       END DO
!      CLOSE(121)
!      PRINT *, 'SURFACE MESH READING COMPLETE'
!      PRINT *, 'ibNodes =', ibNodes, 'ibElems =', ibElems

!      !!$acc update device (xnode, ynode, znode, ibElP1, ibElP2, ibElP3)
!     END SUBROUTINE readSurfaceMeshGambit
!*******************************************************************
end module biocfd_read_input

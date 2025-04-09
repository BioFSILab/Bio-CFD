module biocfd_coefficient_matrix
  use global
  implicit NONE

  contains
SUBROUTINE coefficientMatrix
        USE global
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER (kind = 8) :: i, j, f,nx_var,ny_var,nz_var
        REAL (KIND = 8)    ::  rx1, rx2, rxsum, ry1, ry2, rysum, rz1, rz2, rzsum

         DO f=1,nblocks
         nx_var=block(f)%nx
         ny_var=block(f)%ny
         nz_var=block(f)%nz
         DO j = 2, ny_var+1
            ry1   = block(f)%yp(j)   - block(f)%yp(j-1)
            ry2   = block(f)%yp(j+1) - block(f)%yp(j)
            rysum = ry1 + ry2
            block(f)%Acy(j-1, 1) =   2._rk/(ry1*rysum)
            block(f)%Acy(j-1, 2) =  -2._rk/(ry1*ry2)
            block(f)%Acy(j-1, 3) =   2._rk/(ry2*rysum)
         END DO

         DO i = 2, nx_var+1
            rx1   = block(f)%xp(i)   - block(f)%xp(i-1)
            rx2   = block(f)%xp(i+1) - block(f)%xp(i)
            rxsum = rx1 + rx2
            block(f)%Acx(i-1, 1) =   2._rk/(rx1*rxsum)
            block(f)%Acx(i-1, 2) =  -2._rk/(rx1*rx2)
            block(f)%Acx(i-1, 3) =   2._rk/(rx2*rxsum)
         END DO

         DO i = 2, nz_var+1
            rz1   = block(f)%zp(i)   - block(f)%zp(i-1)
            rz2   = block(f)%zp(i+1) - block(f)%zp(i)
            rzsum = rz1 + rz2
            block(f)%Acz(i-1, 1) =   2._rk/(rz1*rzsum)
            block(f)%Acz(i-1, 2) =  -2._rk/(rz1*rz2)
            block(f)%Acz(i-1, 3) =   2._rk/(rz2*rzsum)
         END DO


         if (f ==1)then
         !inlet, i = 1
         block(f)%Acx(1, 2)   =  block(f)%Acx(1, 2)  - block(f)%Acx(1, 1)
         block(f)%Acx(1, 1)   =  0._rk
         !outlet, i = nx_var
         block(f)%Acx(nx_var, 2)  =  block(f)%Acx(nx_var, 2) - block(f)%Acx(nx_var, 3)
         block(f)%Acx(nx_var, 3)  =  0._rk

         block(f)%Acy(1, 2)   =  block(f)%Acy(1, 2)  + block(f)%Acy(1, 1)
         block(f)%Acy(1, 1)   =  0._rk
         !top, i = ny_var
         block(f)%Acy(ny_var, 2)  =  block(f)%Acy(ny_var, 2) +  block(f)%Acy(ny_var, 3)
         block(f)%Acy(ny_var, 3)  =  0._rk

        	 !front, k = 1
         block(f)%Acz(1, 2)   =  block(f)%Acz(1, 2) + block(f)%Acz(1, 1)
         block(f)%Acz(1, 1)   =  0._rk
         !back, k = nz_var
         block(f)%Acz(nz_var, 2)  =  block(f)%Acz(nz_var, 2) + block(f)%Acz(nz_var, 3)
         block(f)%Acz(nz_var, 3)  =  0._rk
       end if
        END DO
        !!$acc update device (Acx, Acy, Acz)

        print*, "Coefficient Matrix generated"
!!        DO f=1,nblocks
!!         nx_var=block(f)%nx
!!         ny_var=block(f)%ny
!!         nz_var=block(f)%nz

!!       DO k=2,nz_var+1
!!       DO j=2,ny_var+1
!!       DO i=2,nx_var+1
!!          ry1   = block(f)%yp(j)   - block(f)%yp(j-1)
!!          ry2   = block(f)%yp(j+1) - block(f)%yp(j)
!!          rysum = ry1 + ry2
!!          rx1   = block(f)%xp(i)   - block(f)%xp(i-1)
!!          rx2   = block(f)%xp(i+1) - block(f)%xp(i)
!!          rxsum = rx1 + rx2
!!          rz1   = block(f)%zp(k)   - block(f)%zp(k-1)
!!          rz2   = block(f)%zp(k+1) - block(f)%zp(k)
!!          rzsum = rz1 + rz2

!!          ip= nx_var*ny_var*(k-2)+ (j-2)*nx_var + (i-1)
!!          block(f)%Ac(ip, 1) = 2._rk/(rz1*rzsum)
!!          block(f)%Ac(ip, 2) = 2._rk/(ry1*rysum)
!!          block(f)%Ac(ip, 3) = 2._rk/(rx1*rxsum)
!!          block(f)%Ac(ip, 5) = 2._rk/(rx2*rxsum)
!!          block(f)%Ac(ip, 6) = 2._rk/(ry2*rysum)
!!          block(f)%Ac(ip, 7) = 2._rk/(rz2*rzsum)
!!          block(f)%Ac(ip, 4) = -( block(f)%Ac(ip,1) + block(f)%Ac(ip,2) + block(f)%Ac(ip,3) + block(f)%Ac(ip,5)+block(f)% Ac(ip,6)+block(f)%Ac(ip,7))
!!       END DO
!!       END DO
!!       END DO
!!      print*,'inside 1'
!!       if (f .eq.1)then
!!       !Inlet and Outlet..
!!       DO k=2,nz_var+1
!!       DO j=2,ny_var+1
!!          i=2
!!          ip= nx_var*ny_var*(k-2)+(j-2)*nx_var + (i-1)
!!          block(f)%Ac(ip, 4) = block(f)%Ac(ip, 4) - block(f)%Ac(ip, 3)
!!          block(f)%Ac(ip, 3) = 0.

!!          i=nx_var+1
!!          ip= nx_var*ny_var*(k-2)+(j-2)*nx_var + (i-1)
!!          block(f)%Ac(ip, 4) = block(f)%Ac(ip, 4) - block(f)%Ac(ip, 5)
!!          block(f)%Ac(ip, 5) = 0.
!!       END DO
!!       END DO
!!      print*,'inside 2'
!!       !Top and Bottom..


!!       DO k=2,nz_var+1
!!       DO i=2,nx_var+1

!!          j=2
!!          ip= nx_var*ny_var*(k-2)+(j-2)*nx_var + (i-1)
!!          block(f)%Ac(ip, 4) = block(f)%Ac(ip, 4) - block(f)%Ac(ip, 2)
!!          block(f)%Ac(ip, 2) = 0.

!!          ! if ( block(f)% xp(i) .lt. varx1)then
!!          j=ny_var+1
!!          ip= nx_var*ny_var*(k-2)+(j-2)*nx_var + (i-1)
!!          block(f)%Ac(ip, 4) = block(f)%Ac(ip, 4) - block(f)%Ac(ip, 6)
!!          block(f)%Ac(ip, 6) = 0.
!!       ! endif
!!       END DO
!!       END DO



!!       DO j=2,ny_var+1
!!       DO i=2,nx_var+1
!!          k=2
!!          ip= nx_var*ny_var*(k-2)+(j-2)*nx_var + (i-1)
!!          block(f)%Ac(ip, 4) = block(f)%Ac(ip, 4) - block(f)%Ac(ip, 1)
!!          block(f)%Ac(ip, 1) = 0.

!!          k=nz_var+1
!!          ip= nx_var*ny_var*(k-2)+(j-2)*nx_var + (i-1)
!!          block(f)%Ac(ip, 4) = block(f)%Ac(ip, 4) - block(f)%Ac(ip, 7)
!!          block(f)%Ac(ip, 7) = 0.
!!       END DO
!!       END DO
!!      ENDIF
!!      END DO



!
!           DO f=1,nblocks
!               nx_var=block(f)%nx
!               ny_var=block(f)%ny
!               nz_var=block(f)%nz
!
!           DO k=2,nz_var+1
!           DO j=2,ny_var+1
!           DO i=2,nx_var+1
!           counter =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!           !print*,i,j,k
!              if (block(f)%cell(i,j,k)/=0)then
!                 ! print*,counter
!                  block(f)%A(counter,1)=0
!                  block(f)%A(counter,2)=0
!                  block(f)%A(counter,3)=0
!                  block(f)%A(counter,4)=1
!                  block(f)%A(counter,5)=0
!                  block(f)%A(counter,6)=0
!                  block(f)%A(counter,7)=0


!                else

!                ! if (i == nx+1 .OR. j== ny+1 .OR. k==nz+1 .OR. i==2.OR. j==2 .OR. k==2)then
!                !  block(f)%A(counter,1)=Acz(k,1)
!                !  block(f)%A(counter,7)=Acz(k,3)
!                !  block(f)%A(counter,2)=Acy(j,1)
!                !  block(f)%A(counter,6)=Acy(j,3)
!                !  block(f)%A(counter,3)=Acx(i,1)
!                !  block(f)%A(counter,5)=Acx(i,3)
!                !  block(f)%A(counter,4)=Acx(i,2)+Acy(j,2)+Acz(k,2)
!                !  else

!                   block(f)%A(counter,1)=block(f)%Acz(k-1,1)
!                   block(f)%A(counter,7)=block(f)%Acz(k-1,3)
!                   block(f)%A(counter,2)=block(f)%Acy(j-1,1)
!                   block(f)%A(counter,6)=block(f)%Acy(j-1,3)
!                   block(f)%A(counter,3)=block(f)%Acx(i-1,1)
!                   block(f)%A(counter,5)=block(f)%Acx(i-1,3)
!                   block(f)%A(counter,4)=block(f)%Acx(i-1,2)+block(f)%Acy(j-1,2)+block(f)%Acz(k-1,2)
!             !    !print*,i,j,k,counter
!                  end if
!       end do
!       end do
!       end do
!       end do
!
!       amgx_checker=0.
!           DO f=1,nblocks
!               nx_var=block(f)%nx
!               ny_var=block(f)%ny
!               nz_var=block(f)%nz
!           DO k=2,nz_var+1
!           DO j=2,ny_var+1
!           i=2
!           counter =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!           block(f)%A(counter,4)=block(f)%A(counter,4) + block(f)%A(counter,3)
!           block(f)%A(counter,3)= 0.

!           i=nx_var
!           counter =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!           block(f)%A(counter,4)=block(f)%A(counter,4) + block(f)%A(counter,3)
!           block(f)%A(counter,3)= 0.

      end subroutine coefficientMatrix
end module biocfd_coefficient_matrix

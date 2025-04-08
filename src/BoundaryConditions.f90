module biocfd_boundary_conditions
  use global
  implicit none

  contains
!***********************************************************************
SUBROUTINE velocityBC
      USE global
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER (KIND = 8):: i, j, k, g
      REAL (KIND = 8) :: a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4, &
                         a5, b5, c5, a6, b6, c6, a7, b7, c7, a8, b8, c8, &
                         x3, uinlet, varx1
      REAL (KIND = 8) :: a00, qinlet, inletArea, outletArea

    !uc=1._rk


     !!$acc parallel loop gang vector collapse (2) present (ut, vt, wt, u, v, w, xp, yp, zp, deltax, deltay, deltaz)  &
     !!$acc firstprivate (uc, deltat, nx)


      g=1
     !$acc parallel loop gang vector collapse (2) default(present)  &
     !$acc firstprivate (uc, deltat, block(g)%nx)
      DO  k = 2, block(g)%nz+1
      DO  j = 2, block(g)%ny+1
        !uniform inlet
     !   if ( g .ne. 3)then
        !block(g)%ut(1,j,k) = block(g)%ut(block(g)%nx,j,k)
        block(g)%ut(1,j,k) = uc
       !block(g)%vt(1,j,k) = block(g)%vt(block(g)%nx+1,j,k)
       !block(g)%wt(1,j,k) = block(g)%wt(block(g)%nx+1,j,k)
        block(g)%vt(1,j,k) =-block(g)%vt( block(g)%nx+1,j,k)
        block(g)%wt(1,j,k) =-block(g)%wt( block(g)%nx+1,j,k)
        !Neumann - low Re convective flow (outlet)
        !block(g)%ut(block(g)%nx+1,j,k)  =  block(g)%ut(block(g)%nx,j,k)
        ! block(g)%vt(block(g)%nx+2,j,k)  =  block(g)%vt(2,j,k)
        ! block(g)%wt(block(g)%nx+2,j,k)  =  block(g)%wt(2,j,k)
        ! block(g)%ut(block(g)%nx+1,j,k)  =  block(g)%ut(2,j,k)
      !  block(g)%ut(1,j,k) = 0._rk
      ! !block(g)%ut(1,j,k) = block(g)%ut(2,j,k)
      ! ! block(g)%ut(1,j,k) = 1._rk
      ! block(g)%vt(1,j,k) =-block(g)%vt(2,j,k)
      ! block(g)%wt(1,j,k) =-block(g)%wt(2,j,k)
        !Neumann - low Re convective flow (outlet)
        !block(g)%ut(block(g)%nx+1,j,k)  =  block(g)%ut(block(g)%nx,j,k)
     !   endif

     !   if ( g  .ne. 2)then
        !block(g)%ut(block(g)%nx+1,j,k)  =  0._rk
      ! block(g)%ut(block(g)%nx+1,j,k)  =  block(g)%ut(block(g)%nx,j,k)
      ! block(g)%vt(block(g)%nx+2,j,k)  = - block(g)%vt(block(g)%nx+1,j,k)
      ! block(g)%wt(block(g)%nx+2,j,k)  = - block(g)%wt(block(g)%nx+1,j,k)
        !Orlanski - vortex shedding Re Convective flow (outlet)
        block(g)%ut(block(g)%nx+1,j,k)  =  block(g)%u(block(g)%nx+1,j,k)-(deltat/block(g)%deltax(block(g)%nx+1))*uc*(block(g)%u(block(g)%nx+1,j,k)-block(g)%u(block(g)%nx,j,k))
        block(g)%vt(block(g)%nx+2,j,k)  = -block(g)%vt(block(g)%nx+1,j,k) + block(g)%v(block(g)%nx+2,j,k) + block(g)%v(block(g)%nx+1,j,k) - (2._rk*deltat/block(g)%deltax(block(g)%nx+2))*uc*(block(g)%v(block(g)%nx+2,j,k)-block(g)%v(block(g)%nx+1,j,k))
        block(g)%wt(block(g)%nx+2,j,k)  = -block(g)%wt(block(g)%nx+1,j,k) + block(g)%w(block(g)%nx+2,j,k) + block(g)%w(block(g)%nx+1,j,k) - (2._rk*deltat/block(g)%deltax(block(g)%nx+2))*uc*(block(g)%w(block(g)%nx+2,j,k)-block(g)%w(block(g)%nx+1,j,k))
      !  endif
        END DO
        END DO
     !$acc end parallel


     !!$acc parallel loop gang vector collapse (2) present (ut, vt, wt) &
     !!$acc firstprivate (ny)

     !$acc parallel loop gang vector collapse (2) default(present) &
     !$acc firstprivate (block(g)%ny)
      DO k = 2, block(g)%nz+1
      DO i = 2, block(g)%nx+1

          block(g)%ut(i,1,k) = block(g)%ut(i,2,k)
          block(g)%wt(i,1,k) = block(g)%wt(i,2,k)
        !block(g)%ut(i,1,k) =  block(g)%ut(i,2,k)
        !block(g)%wt(i,1,k) =  block(g)%wt(i,2,k)
         block(g)%vt(i,2,k) =  0._rk

          block(g)%ut(i,block(g)%ny+2,k) = block(g)%ut(i,block(g)%ny+1,k)                                      !wall no slip - closed channel
          block(g)%wt(i,block(g)%ny+2,k) = block(g)%wt(i,block(g)%ny+1,k)
        !block(g)%ut(i,block(g)%ny+2,k) =  block(g)%ut(i,block(g)%ny+1,k)                                      !symmetric bc - open channel
        !block(g)%wt(i,block(g)%ny+2,k) =  block(g)%wt(i,block(g)%ny+1,k)
         block(g)%vt(i,block(g)%ny+1,k) =  0._rk



        END DO
        END DO
     !$acc end parallel

     !$acc parallel loop gang vector collapse (2) default(present) &
     !$acc firstprivate (block(g)%nz)
      DO  j = 2, block(g)%ny+1
      DO  i = 2, block(g)%nx+1
         block(g)%ut(i,j,1) =  block(g)%ut(i,j,2)
         block(g)%vt(i,j,1) =  block(g)%vt(i,j,2)
        !block(g)%ut(i,j,1) =  block(g)%ut(i,j,2)
        !block(g)%vt(i,j,1) =  block(g)%vt(i,j,2)
         block(g)%wt(i,j,1) =  0._rk

         block(g)%ut(i,j,block(g)%nz+2) =  block(g)%ut(i,j,block(g)%nz+1)                                      !wall no slip - closed channel
         block(g)%vt(i,j,block(g)%nz+2) =  block(g)%vt(i,j,block(g)%nz+1)
        !block(g)%ut(i,j,block(g)%nz+2) =  block(g)%ut(i,j,block(g)%nz+1)                                      !symmetric bc - open channel
        !block(g)%vt(i,j,block(g)%nz+2) =  block(g)%vt(i,j,block(g)%nz+1)
         block(g)%wt(i,j,block(g)%nz+1) =  0._rk

      END DO
      END DO
      !$acc end parallel
      END SUBROUTINE velocityBC


!***********************************************************************
      SUBROUTINE solidCellBC
         USE global
         INTEGER, PARAMETER :: rk = selected_real_kind(8)
         INTEGER (KIND = 8):: i, j, k, n, g
        !!$acc parallel loop gang vector &
        !!$acc present(solidIndexPtr, ut, vt, wt, p) &
        !!$acc private (i, j, k)
         !DO g=1,nblocks
         DO g=blk_start,nblocks
        !$acc parallel loop gang vector &
        !$acc default(present) &
        !$acc private (i, j, k)
         DO n = 1, block(g)%solidCellCount
            i = block(g)%solidIndexPtr(n,1)
            j = block(g)%solidIndexPtr(n,2)
	     k = block(g)%solidIndexPtr(n,3)
            block(g)%ut(i,j,k) = 0._rk
            block(g)%vt(i,j,k) = 0._rk
            block(g)%wt(i,j,k) = 0._rk
            block(g)%p(i,j,k)  = 0._rk
         END DO
        !$acc end parallel
         END DO
        ! END DO
      END SUBROUTINE solidCellBC
!***********************************************************************

!**************************************************************
      !SUBROUTINE solidCellBC_move(b_blk_no)
      SUBROUTINE solidCellBC_move(g)
         USE global
         INTEGER, PARAMETER :: rk = selected_real_kind(8)
         INTEGER (KIND = 8):: i, j, k,n, vbcOption
         !INTEGER :: g
         INTEGER (KIND = 8),INTENT(IN):: g
         ! INTEGER (KIND = 8), INTENT(IN) :: b_blk_no
        !g=b_blk_no
      !  print*,g,'inside_solid'
         if (block(g)%move_check == 1) then
         !INTEGER (KIND = 8):: g
         !DO g=blk_start, nblocks
         !$acc parallel loop gang vector &
         !$acc default(present)private (i, j)
         !!$acc present (ut, vt, p) private (i, j)
         DO n = 1, block(g)%solidCellCount
            i = block(g)%solidIndexPtr(n,1)
            j = block(g)%solidIndexPtr(n,2)
            k = block(g)%solidIndexPtr(n,3)
            block(g)%ut(i,j,k) = 0._rk
            block(g)%vt(i,j,k) = 0._rk
            block(g)%wt(i,j,k) = 0._rk
            block(g)%p(i,j,k) = 0._rk
            block(g)%u(i,j,k) = 0._rk
            block(g)%v(i,j,k) = 0._rk
            block(g)%w(i,j,k) = 0._rk
         END DO
         !$acc end parallel
        endif
        !END DO
      END SUBROUTINE solidCellBC_move

end module biocfd_boundary_conditions

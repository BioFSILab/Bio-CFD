module biocfd_boundary_conditions
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global
  implicit none
  private

  public :: velocityBC, solidCellBC, solidCellBC_move

  contains

SUBROUTINE velocityBC
      INTEGER (int64):: i, j, k, g

      g=1
     !$acc parallel loop gang vector collapse (2) default(present)  &
     !$acc firstprivate (uc, deltat, block(g)%nx)
      DO  k = 2, block(g)%nz+1
      DO  j = 2, block(g)%ny+1
        !uniform inlet
        block(g)%ut(1,j,k) = uc
        block(g)%vt(1,j,k) =-block(g)%vt( block(g)%nx+1,j,k)
        block(g)%wt(1,j,k) =-block(g)%wt( block(g)%nx+1,j,k)
        !Orlanski - vortex shedding Re Convective flow (outlet)
      block(g)%ut(block(g)%nx+1,j,k) = block(g)%u(block(g)%nx+1,j,k)- &
                                      (deltat/block(g)%deltax(block(g)%nx+1))*&
                                      uc*(block(g)%u(block(g)%nx+1,j,k)-block(g)%u(block(g)%nx,j,k))
      block(g)%vt(block(g)%nx+2,j,k) = -block(g)%vt(block(g)%nx+1,j,k)+ &
                                        block(g)%v(block(g)%nx+2,j,k)+ &
                                        block(g)%v(block(g)%nx+1,j,k)- &
                                        (2._dp*deltat/block(g)%deltax(block(g)%nx+2))*&
                                    uc*(block(g)%v(block(g)%nx+2,j,k)-block(g)%v(block(g)%nx+1,j,k))
      block(g)%wt(block(g)%nx+2,j,k) = -block(g)%wt(block(g)%nx+1,j,k)+ &
                                       block(g)%w(block(g)%nx+2,j,k)+ &
                                       block(g)%w(block(g)%nx+1,j,k)- &
                                       (2._dp*deltat/block(g)%deltax(block(g)%nx+2))*&
                                    uc*(block(g)%w(block(g)%nx+2,j,k)-block(g)%w(block(g)%nx+1,j,k))
        END DO
        END DO
     !$acc end parallel

     !$acc parallel loop gang vector collapse (2) default(present) &
     !$acc firstprivate (block(g)%ny)
      DO k = 2, block(g)%nz+1
      DO i = 2, block(g)%nx+1

          block(g)%ut(i,1,k) = block(g)%ut(i,2,k)
          block(g)%wt(i,1,k) = block(g)%wt(i,2,k)
        !block(g)%ut(i,1,k) =  block(g)%ut(i,2,k)
        !block(g)%wt(i,1,k) =  block(g)%wt(i,2,k)
         block(g)%vt(i,2,k) =  0._dp

          block(g)%ut(i,block(g)%ny+2,k) = block(g)%ut(i,block(g)%ny+1,k)                                      !wall no slip - closed channel
          block(g)%wt(i,block(g)%ny+2,k) = block(g)%wt(i,block(g)%ny+1,k)
        !block(g)%ut(i,block(g)%ny+2,k) =  block(g)%ut(i,block(g)%ny+1,k)                                      !symmetric bc - open channel
        !block(g)%wt(i,block(g)%ny+2,k) =  block(g)%wt(i,block(g)%ny+1,k)
         block(g)%vt(i,block(g)%ny+1,k) =  0._dp



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
         block(g)%wt(i,j,1) =  0._dp

         block(g)%ut(i,j,block(g)%nz+2) =  block(g)%ut(i,j,block(g)%nz+1)                                      !wall no slip - closed channel
         block(g)%vt(i,j,block(g)%nz+2) =  block(g)%vt(i,j,block(g)%nz+1)
        !block(g)%ut(i,j,block(g)%nz+2) =  block(g)%ut(i,j,block(g)%nz+1)                                      !symmetric bc - open channel
        !block(g)%vt(i,j,block(g)%nz+2) =  block(g)%vt(i,j,block(g)%nz+1)
         block(g)%wt(i,j,block(g)%nz+1) =  0._dp

      END DO
      END DO
      !$acc end parallel
      END SUBROUTINE velocityBC

      SUBROUTINE solidCellBC
         INTEGER (int64):: i, j, k, n, g
         DO g=blk_start,nblocks
        !$acc parallel loop gang vector &
        !$acc default(present) &
        !$acc private (i, j, k)
         DO n = 1, block(g)%solidCellCount
            i = block(g)%solidIndexPtr(n,1)
            j = block(g)%solidIndexPtr(n,2)
            k = block(g)%solidIndexPtr(n,3)
            block(g)%ut(i,j,k) = 0._dp
            block(g)%vt(i,j,k) = 0._dp
            block(g)%wt(i,j,k) = 0._dp
            block(g)%p(i,j,k)  = 0._dp
         END DO
        !$acc end parallel
         END DO
      END SUBROUTINE solidCellBC

      SUBROUTINE solidCellBC_move(g)
         INTEGER (int64):: i, j, k,n
         INTEGER (int64),INTENT(IN):: g

         if (block(g)%move_check == 1) then
         !$acc parallel loop gang vector &
         !$acc default(present)private (i, j)
         DO n = 1, block(g)%solidCellCount
            i = block(g)%solidIndexPtr(n,1)
            j = block(g)%solidIndexPtr(n,2)
            k = block(g)%solidIndexPtr(n,3)
            block(g)%ut(i,j,k) = 0._dp
            block(g)%vt(i,j,k) = 0._dp
            block(g)%wt(i,j,k) = 0._dp
            block(g)%p(i,j,k) = 0._dp
            block(g)%u(i,j,k) = 0._dp
            block(g)%v(i,j,k) = 0._dp
            block(g)%w(i,j,k) = 0._dp
         END DO
         !$acc end parallel
        endif
      END SUBROUTINE solidCellBC_move

end module biocfd_boundary_conditions

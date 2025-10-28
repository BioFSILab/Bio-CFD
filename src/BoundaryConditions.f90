module biocfd_boundary_conditions
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use biocfd_block_type,only : Block_t
  implicit none
  private

  public :: velocityBC, solidCellBC, solidCellBC_move

  contains

SUBROUTINE velocityBC(blk,deltat,uc)
      type(Block_t), intent(inout) :: blk
      INTEGER (int64):: i, j, k
      REAL(dp), intent(in) :: deltat,uc

     !$acc parallel loop gang vector collapse (2) default(present)  &
     !$acc firstprivate (uc, deltat)
      DO  k = 2, blk%nz+1
      DO  j = 2, blk%ny+1
        !uniform inlet
        blk%ut(1,j,k) = uc
        blk%vt(1,j,k) =-blk%vt(blk%nx+1,j,k)
        blk%wt(1,j,k) =-blk%wt(blk%nx+1,j,k)
        !Orlanski - vortex shedding Re Convective flow (outlet)
      blk%ut(blk%nx+1,j,k) = blk%u(blk%nx+1,j,k)- &
                                      (deltat/blk%deltax(blk%nx+1))*&
                                      uc*(blk%u(blk%nx+1,j,k)-blk%u(blk%nx,j,k))
      blk%vt(blk%nx+2,j,k) = -blk%vt(blk%nx+1,j,k)+ &
                                        blk%v(blk%nx+2,j,k)+ &
                                        blk%v(blk%nx+1,j,k)- &
                                        (2._dp*deltat/blk%deltax(blk%nx+2))*&
                                    uc*(blk%v(blk%nx+2,j,k)-blk%v(blk%nx+1,j,k))
      blk%wt(blk%nx+2,j,k) = -blk%wt(blk%nx+1,j,k)+ &
                                       blk%w(blk%nx+2,j,k)+ &
                                       blk%w(blk%nx+1,j,k)- &
                                       (2._dp*deltat/blk%deltax(blk%nx+2))*&
                                    uc*(blk%w(blk%nx+2,j,k)-blk%w(blk%nx+1,j,k))
        END DO
        END DO
     !$acc end parallel loop

     !$acc parallel loop gang vector collapse (2) default(present)
      DO k = 2, blk%nz+1
      DO i = 2, blk%nx+1

          blk%ut(i,1,k) = blk%ut(i,2,k)
          blk%wt(i,1,k) = blk%wt(i,2,k)
          blk%vt(i,2,k) =  0._dp

          blk%ut(i,blk%ny+2,k) = blk%ut(i,blk%ny+1,k)                                      !wall no slip - closed channel
          blk%wt(i,blk%ny+2,k) = blk%wt(i,blk%ny+1,k)
          blk%vt(i,blk%ny+1,k) =  0._dp



        END DO
        END DO
     !$acc end parallel loop

     !$acc parallel loop gang vector collapse (2) default(present)
      DO  j = 2, blk%ny+1
      DO  i = 2, blk%nx+1
         blk%ut(i,j,1) =  blk%ut(i,j,2)
         blk%vt(i,j,1) =  blk%vt(i,j,2)
         blk%wt(i,j,1) =  0._dp

         blk%ut(i,j,blk%nz+2) =  blk%ut(i,j,blk%nz+1)                                      !wall no slip - closed channel
         blk%vt(i,j,blk%nz+2) =  blk%vt(i,j,blk%nz+1)
         blk%wt(i,j,blk%nz+1) =  0._dp

      END DO
      END DO
      !$acc end parallel loop
      END SUBROUTINE velocityBC

      SUBROUTINE solidCellBC(blk)
        type(Block_t), intent(inout) :: blk
        INTEGER (int64):: i, j, k, n
        !$acc parallel loop gang vector &
        !$acc default(present) &
        !$acc private (i, j, k)
         DO n = 1, blk%solidCellCount
            i = blk%solidIndexPtr(n,1)
            j = blk%solidIndexPtr(n,2)
            k = blk%solidIndexPtr(n,3)
            blk%ut(i,j,k) = 0._dp
            blk%vt(i,j,k) = 0._dp
            blk%wt(i,j,k) = 0._dp
            blk%p(i,j,k)  = 0._dp
         END DO
        !$acc end parallel loop
      END SUBROUTINE solidCellBC

      SUBROUTINE solidCellBC_move(blk)
        type(Block_t), intent(inout) :: blk
        INTEGER (int64):: i, j, k,n

         if (blk%move_check == 1) then
         !$acc parallel loop gang vector &
         !$acc default(present)private (i, j)
         DO n = 1, blk%solidCellCount
            i = blk%solidIndexPtr(n,1)
            j = blk%solidIndexPtr(n,2)
            k = blk%solidIndexPtr(n,3)
            blk%ut(i,j,k) = 0._dp
            blk%vt(i,j,k) = 0._dp
            blk%wt(i,j,k) = 0._dp
            blk%p(i,j,k) = 0._dp
            blk%u(i,j,k) = 0._dp
            blk%v(i,j,k) = 0._dp
            blk%w(i,j,k) = 0._dp
         END DO
         !$acc end parallel loop
        endif
      END SUBROUTINE solidCellBC_move

end module biocfd_boundary_conditions

module biocfd_coefficient_matrix
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, nblocks
  implicit NONE
  private

  public :: coefficientMatrix

  contains
SUBROUTINE coefficientMatrix

        INTEGER (int64) :: i, j, f,nx_var,ny_var,nz_var
        REAL (dp)    ::  rx1, rx2, rxsum, ry1, ry2, rysum, rz1, rz2, rzsum

         DO f=1,nblocks
         nx_var=block(f)%nx
         ny_var=block(f)%ny
         nz_var=block(f)%nz
         DO j = 2, ny_var+1
            ry1   = block(f)%yp(j)   - block(f)%yp(j-1)
            ry2   = block(f)%yp(j+1) - block(f)%yp(j)
            rysum = ry1 + ry2
            block(f)%Acy(j-1, 1) =   2._dp/(ry1*rysum)
            block(f)%Acy(j-1, 2) =  -2._dp/(ry1*ry2)
            block(f)%Acy(j-1, 3) =   2._dp/(ry2*rysum)
         END DO

         DO i = 2, nx_var+1
            rx1   = block(f)%xp(i)   - block(f)%xp(i-1)
            rx2   = block(f)%xp(i+1) - block(f)%xp(i)
            rxsum = rx1 + rx2
            block(f)%Acx(i-1, 1) =   2._dp/(rx1*rxsum)
            block(f)%Acx(i-1, 2) =  -2._dp/(rx1*rx2)
            block(f)%Acx(i-1, 3) =   2._dp/(rx2*rxsum)
         END DO

         DO i = 2, nz_var+1
            rz1   = block(f)%zp(i)   - block(f)%zp(i-1)
            rz2   = block(f)%zp(i+1) - block(f)%zp(i)
            rzsum = rz1 + rz2
            block(f)%Acz(i-1, 1) =   2._dp/(rz1*rzsum)
            block(f)%Acz(i-1, 2) =  -2._dp/(rz1*rz2)
            block(f)%Acz(i-1, 3) =   2._dp/(rz2*rzsum)
         END DO


         if (f ==1)then
         !inlet, i = 1
         block(f)%Acx(1, 2)   =  block(f)%Acx(1, 2)  - block(f)%Acx(1, 1)
         block(f)%Acx(1, 1)   =  0._dp
         !outlet, i = nx_var
         block(f)%Acx(nx_var, 2)  =  block(f)%Acx(nx_var, 2) - block(f)%Acx(nx_var, 3)
         block(f)%Acx(nx_var, 3)  =  0._dp

         block(f)%Acy(1, 2)   =  block(f)%Acy(1, 2)  + block(f)%Acy(1, 1)
         block(f)%Acy(1, 1)   =  0._dp
         !top, i = ny_var
         block(f)%Acy(ny_var, 2)  =  block(f)%Acy(ny_var, 2) +  block(f)%Acy(ny_var, 3)
         block(f)%Acy(ny_var, 3)  =  0._dp

        	 !front, k = 1
         block(f)%Acz(1, 2)   =  block(f)%Acz(1, 2) + block(f)%Acz(1, 1)
         block(f)%Acz(1, 1)   =  0._dp
         !back, k = nz_var
         block(f)%Acz(nz_var, 2)  =  block(f)%Acz(nz_var, 2) + block(f)%Acz(nz_var, 3)
         block(f)%Acz(nz_var, 3)  =  0._dp
       end if
        END DO

        print*, "Coefficient Matrix generated"

      end subroutine coefficientMatrix
end module biocfd_coefficient_matrix

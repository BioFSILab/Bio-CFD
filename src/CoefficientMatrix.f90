module biocfd_coefficient_matrix
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use biocfd_blocks, only: Blocks
  implicit NONE
  private

  public :: coefficientMatrix

contains
  SUBROUTINE coefficientMatrix(blk, is_coarse)

    !> The block to construct the coefficient matrix for
    type(Blocks), intent(inout) :: blk
    logical, intent(in) :: is_coarse
    INTEGER (int64) :: i, j, nx_var, ny_var, nz_var
    REAL (dp):: rx1, rx2, rxsum, ry1, ry2, rysum, rz1, rz2, rzsum

    nx_var=blk%nx
    ny_var=blk%ny
    nz_var=blk%nz
    DO j = 2, ny_var+1
       ry1   = blk%yp(j)   - blk%yp(j-1)
       ry2   = blk%yp(j+1) - blk%yp(j)
       rysum = ry1 + ry2
       blk%Acy(j-1, 1) =   2._dp/(ry1*rysum)
       blk%Acy(j-1, 2) =  -2._dp/(ry1*ry2)
       blk%Acy(j-1, 3) =   2._dp/(ry2*rysum)
    END DO

    DO i = 2, nx_var+1
       rx1   = blk%xp(i)   - blk%xp(i-1)
       rx2   = blk%xp(i+1) - blk%xp(i)
       rxsum = rx1 + rx2
       blk%Acx(i-1, 1) =   2._dp/(rx1*rxsum)
       blk%Acx(i-1, 2) =  -2._dp/(rx1*rx2)
       blk%Acx(i-1, 3) =   2._dp/(rx2*rxsum)
    END DO

    DO i = 2, nz_var+1
       rz1   = blk%zp(i)   - blk%zp(i-1)
       rz2   = blk%zp(i+1) - blk%zp(i)
       rzsum = rz1 + rz2
       blk%Acz(i-1, 1) =   2._dp/(rz1*rzsum)
       blk%Acz(i-1, 2) =  -2._dp/(rz1*rz2)
       blk%Acz(i-1, 3) =   2._dp/(rz2*rzsum)
    END DO

    ! Some additional processing happens for the coarse block
    if (is_coarse) then
       !inlet, i = 1
       blk%Acx(1, 2)   =  blk%Acx(1, 2)  - blk%Acx(1, 1)
       blk%Acx(1, 1)   =  0._dp
       !outlet, i = nx_var
       blk%Acx(nx_var, 2)  =  blk%Acx(nx_var, 2) - blk%Acx(nx_var, 3)
       blk%Acx(nx_var, 3)  =  0._dp

       blk%Acy(1, 2)   =  blk%Acy(1, 2)  + blk%Acy(1, 1)
       blk%Acy(1, 1)   =  0._dp
       !top, i = ny_var
       blk%Acy(ny_var, 2)  =  blk%Acy(ny_var, 2) +  blk%Acy(ny_var, 3)
       blk%Acy(ny_var, 3)  =  0._dp

       !front, k = 1
       blk%Acz(1, 2)   =  blk%Acz(1, 2) + blk%Acz(1, 1)
       blk%Acz(1, 1)   =  0._dp
       !back, k = nz_var
       blk%Acz(nz_var, 2)  =  blk%Acz(nz_var, 2) + blk%Acz(nz_var, 3)
       blk%Acz(nz_var, 3)  =  0._dp
    end if

  end subroutine coefficientMatrix
end module biocfd_coefficient_matrix

module biocfd_initial_conditions
  ! allow(use-all) - TODO: Attempt to fix this in the future
  use global
  implicit none
  private

  public :: initialConditions

  contains
      SUBROUTINE initialConditions
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER::  i, j, k, n, f
        WRITE(*,*) 'Enter initialcondtitions'
       !cell variables
        DO i=1,nblocks
        block(i)%u = uc
         block(i)%ut = uc
         block(i)%v = 0._rk
         block(i)%vt = 0._rk
        block(i)%w = 0._rk
        block(i)%wt = 0._rk
        block(i)%p = 0._rk
        block(i)%u2_sum = 0._rk
        block(i)%v2_sum = 0._rk
        block(i)%w2_sum = 0._rk
        block(i)%p2_sum = 0._rk
        block(i)%u2_avg = 0._rk
        block(i)%v2_avg = 0._rk
        block(i)%w2_avg = 0._rk
        block(i)%p2_avg = 0._rk
        block(i)%u_sum = 0._rk
        block(i)%v_sum = 0._rk
        block(i)%w_sum = 0._rk
        block(i)%p_sum = 0._rk
        block(i)%uv_sum = 0._rk
        block(i)%vw_sum = 0._rk
        block(i)%uw_sum = 0._rk
        block(i)%uv_avg = 0._rk
        block(i)%vw_avg = 0._rk
        block(i)%uw_avg = 0._rk
        block(i)%u_avg = 0._rk
        block(i)%v_avg = 0._rk
        block(i)%w_avg = 0._rk
        block(i)%p_avg = 0._rk
        block(i)%ufl = 0._rk
        block(i)%vfl = 0._rk
        block(i)%wfl = 0._rk
        block(i)%resi_u = 0._rk
        block(i)%resi_v = 0._rk
        block(i)%resi_w = 0._rk
        block(i)%cell_pr=0
        block(i)%uflu_avg = 0._rk
        block(i)%vflu_avg = 0._rk
        block(i)%wflu_avg = 0._rk
        block(i)%pflu_avg = 0._rk
        block(i)%uflu_rms = 0._rk
        block(i)%vflu_rms = 0._rk
        block(i)%wflu_rms = 0._rk
        block(i)%pflu_rms = 0._rk
        block(i)%uvflu_avg = 0._rk
        block(i)%vwflu_avg = 0._rk
        block(i)%uwflu_avg = 0._rk
        end do
        SUMWSS = 0._rk            !SUMWSS global real array(nsurf)
        SIGNWSS = 0._rk           !SIGNWSS global real array(nsurf)
        ita = 0
        ita1 = 0
        totime = 0.
        DO f=1,nblocks
           DO n = 1, block(f)%fluidCellCount
              i = block(f)%fluidIndexPtr(n, 1)
              j = block(f)%fluidIndexPtr(n, 2)
              k = block(f)%fluidIndexPtr(n, 3)
              block(f)% u(i,j,k)  = uc  !396.33054782262406 !116.236233
              block(f)%v(i,j,k)  = 0.
              block(f)%w(i,j,k)  = 0._rk
              block(f)%ut(i,j,k) = uc
              block(f)%vt(i,j,k) = 0._rk
              block(f)%wt(i,j,k) = 0._rk
          END DO
        END DO

        print*, 'initial'

      END SUBROUTINE initialConditions

end module biocfd_initial_conditions

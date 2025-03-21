
      SUBROUTINE lastConditions
       USE global
       implicit none
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER::  i, j, k, g
        CHARACTER(len=150) :: filename3
        WRITE(*,*) 'Enter lastcondtitions'

        DO g=1,nblocks
        block(g)%u_sum = 0._rk
        block(g)%v_sum = 0._rk
        block(g)%w_sum = 0._rk
        block(g)% p_sum = 0._rk
        block(g)% u_avg = 0._rk
        block(g)%v_avg = 0._rk
        block(g)%w_avg = 0._rk
        block(g)% p_avg = 0._rk
        block(g)%resi_u = 0._rk
        block(g)% resi_v = 0._rk
        block(g)% resi_w = 0._rk
        END DO
        SUMWSS = 0._rk            !SUMWSS global real array(nsurf)
        SIGNWSS = 0._rk           !SIGNWSS global real array(nsurf)
        ita = 0
        ita1 = 0
        ita2 = 0


        DO g=1,nblocks
       WRITE(filename3,3) g, re
  3     FORMAT('out/aorta_chkpt.',i3.3,'.',f6.1,".dat")
       OPEN (1, FILE=filename3, FORM='formatted')
 	DO 30 k = 1, block(g)%nz+2
       DO 30 j = 1, block(g)%ny+2
       DO 30 i = 1, block(g)%nx+2
         READ(1,*) block(g)%u(i,j,k), block(g)%v(i,j,k), block(g)%w(i,j,k), block(g)%p(i,j,k), totime, ita, ita1
 30    CONTINUE
 	CLOSE(1)
        END DO


       DO g=1,nblocks
       DO k = 1, block(g)%nz+2
       DO j = 1, block(g)%ny+2
       DO i = 1, block(g)%nx+2
           block(g)%ut(i,j,k) =  block(g)%u(i,j,k)
           block(g)%vt(i,j,k) =  block(g)%v(i,j,k)
           block(g)%wt(i,j,k) =  block(g)%w(i,j,k)

        END DO
        END DO
        END DO
        END DO

      !!$acc update device(u, v, w, ut, vt, wt, p, resi_u, resi_v, resi_w)

     !  print*,'after read check point'
     !  call readComputeSumData
     !  call readComputeSumSqData
     !  print*,'after read compute'
     !  call readstressData
     !  print*,'after read stress'
       !ita=0
       !ita1=0
       !ita2=0
       !totime=0
       !ita2=ita

      END SUBROUTINE lastConditions





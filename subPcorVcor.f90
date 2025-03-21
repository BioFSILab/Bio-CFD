!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE poissonSolver
        USE global
        IMPLICIT NONE
        INTEGER(KIND=8) :: i, j,k, n, g, gg1, f, co
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        REAL (KIND = 8)    :: dalt, div, dab, dudt, dvdt, dwdt
        REAL (KIND = 8)    :: max_derr1, max_derr2, max_div, max_derrStdSt
        REAL (KIND = 8)    :: er_dudt, er_dvdt, er_dwdt, err_ds
        INTEGER(KIND=8) :: max_nIterPcor, max_nit
        CHARACTER*160 filename1                 
        


      ! DO g=1,nblocks
      !  !ALLOCATE ( block(g)%b(block(g)%nx+2, block(g)%ny+2, block(g)%nz+2) )
      !  ALLOCATE ( block(g)%b(block(g)%nx+2, block(g)%ny+2, block(g)%nz+2) )
      ! END DO
          max_derrStdst=0._rk
          max_derr1=0._rk
          max_derr2=0._rk
          max_nIterPcor=0
          max_div=0._rk
          max_nit=0._rk      
          err_ds=0._rk      
          er_dudt=0.
          er_dvdt=0.
          er_dwdt=0.
         !derrStdSt = 0._rk

         !initialize variables
         !nIterPcor = 0  
         !divmax = 0._rk
         !dalt   = 0._rk         
         !derr1  = 0._rk
        
        
        DO g=1,nblocks
        !$acc parallel loop gang vector collapse (3) default(present)
        DO k = 1, block(g)%nz+2
        DO j = 1, block(g)%ny+2
        DO i = 1, block(g)%nx+2
           block(g)%b(i,j,k)  = 0.
           block(g)%pc(i,j,k) = 0.
           block(g)%pco(i,j,k)= 0.
           
        END DO
        END DO
        END DO
        !$acc end parallel
        block(g)%nIterPcor=0
        !block(g)%divmax = 0._rk
        !block(g)%dalt   = 0._rk
        block(g)%derr2  = 0._rk
        block(g)%derrStdSt=0._rk
        block(g)%rhs=0.0
        end do
        !print*,1

!       if ( ita .eq. 258)then
!               call writeOutput1
!               print*,'before swap in ps'
!              ! pause
!       endif
       
        solverTime=0.
     CALL cpu_time(dStart)
        CALL fineUpdate_newv_bd
       !CALL writeOutput1
       !pause
        CALL coarseUpdate_newv
      ! if ( ita .eq. 258)then
      !         call writeOutput1
      !         print*,'after swap in ps'
      !         pause
      ! endif
       
!      CALL writeOutput1
!   print*,ita,'before rbsor'
!      !pause
     CALL cpu_time(dfinish)
        coupTime=coupTime + dfinish -dstart 
        DO g=1,nblocks
               CALL computeDiv(g)    !divergence vector 
        END DO
        DO g=2,nblocks
         CALL cpu_time(dStart)
         amgxita=0
         !CALL REDBLACKSOR(g)
         CALL REDBLACKSOR_linear(g)
        CALL cpu_time(dfinish)
         msTime = msTime + dfinish-dstart
        end do
        CALL coarseUpdate_pc
        g=1
        CALL cpu_time(dStart)
        !CALL AMGX_SOLVER(g)
        CALL REDBLACKSOR_linear(g)

        CALL cpu_time(dfinish)
        call fineUpdate_pc_bd
        DO g=2,nblocks
         CALL cpu_time(dStart)
         amgxita=0
         CALL REDBLACKSOR_linear(g)
        CALL cpu_time(dfinish)
         msTime = msTime + dfinish-dstart
        end do
          CALL coarseUpdate_pc

!       if ( ita .eq. 258)then
!               call writeOutput1
!               print*,'after rbsor in ps'
!               !pause
!       endif

!    call cpu_time(dstart)

       DO g=1,nblocks 
               CALL correctPressure(g)!pressure correction
               CALL correctVelocity(g) !velocity correction

        END DO
         CALL velocityBC      !correct velocity at boundaries
!       CALL writeOutput1
!   print*,ita,'after rbsor'
!       !pause
        
        !$omp parallel do private( n,i,j,k,er_dudt,er_dvdt,er_dwdt, err_ds,g) num_threads(3)
         DO g=1,nblocks
         err_ds=0.
        !$acc parallel loop gang vector firstprivate (deltat)   &
        !$acc private (i, j, k, er_dudt, er_dvdt, er_dwdt)               &
        !$acc default(present) reduction (max: err_ds)
         DO n = 1, block(g)%fluidCellCount
           i = block(g)%fluidIndexPtr(n, 1)
           j = block(g)%fluidIndexPtr(n, 2)  
           k = block(g)%fluidIndexPtr(n, 3)  
           !block(g)%dudt = dabs((block(g)%ut(i,j,k) - block(g)%u(i,j,k)))/deltat
           !block(g)%dvdt = dabs((block(g)%vt(i,j,k) - block(g)%v(i,j,k)))/deltat
           !block(g)% dwdt = dabs((block(g)%wt(i,j,k) - block(g)%w(i,j,k)))/deltat
           !block(g)%derrStdSt = dmax1(block(g)%derrStdSt, block(g)%dudt,block(g)%dvdt, block(g)%dwdt)    
           er_dudt = dabs((block(g)%ut(i,j,k) - block(g)%u(i,j,k)))/deltat
           er_dvdt = dabs((block(g)%vt(i,j,k) - block(g)%v(i,j,k)))/deltat
           er_dwdt = dabs((block(g)%wt(i,j,k) - block(g)%w(i,j,k)))/deltat
           err_ds = dmax1(err_ds, er_dudt,er_dvdt, er_dwdt)    
         ENDDO
         !$acc end parallel
           block(g)%derrStdSt = err_ds    
         
        
         ENDDO
        !$omp end parallel do



        DO i=1,nblocks
               if ( block(i)%derr2 >max_derr2)then
                  max_derr2=block(i)%derr2
               end if
              if ( block(i)%derr1 >max_derr1)then
                 max_derr1=block(i)%derr1
              end if
              if ( block(i)%derrStdSt >max_derrStdSt)then
                 max_derrStdSt=block(i)%derrStdSt
              end if
              if ( block(i)%nIterPcor >max_nIterPcor)then
                 max_nIterPcor=block(i)%nIterPcor
              end if
              totalTime=totime + totalTime
              !if ( block(i)%divmax >max_divmax)then
              !   max_divmax=block(i)%divmax
              !end if
              if ( block(i)%nit >max_nit)then
                 max_nit=block(i)%nit
              end if
          end do

            WRITE(filename1,1) 
 1          FORMAT('sphere_iter.dat')
         OPEN(111,FILE=filename1,ACCESS='Append',STATUS='unknown')
         !WRITE(111,126)   ita, block(1)%nIterPcor, block(2)%nIterPcor, block(3)%nIterPcor,max_derr2, max_derrStdSt, dmid(1),dmid(2), dmid(3)
         !WRITE(111,126)   ita, block(1)%nIterPcor, block(2)%nIterPcor, block(3)%nIterPcor, omega1, omega2,omega3, totime
         WRITE(111,126)   ita, block(1)%nIterPcor, block(2)%nIterPcor, omega1, omega2, solverTime
        !WRITE(111,126) ita, block(1)%nIterPcor,  block(2)%nIterPcor,  block(3)%nIterPcor, max_derr2, max_derrStdSt, dfinish-dstart
         WRITE(*,16) ita, max_nIterPcor, max_derr2, max_derrStdSt, totalTime
         !WRITE(*,16) ita, max_nit, max_derr2, max_derrStdSt, totalTime
!126      FORMAT(' ',I8, 3I10, 7E15.6)
 126      FORMAT(' ',I8, 2I10, 2F6.2,F14.9)
 16      FORMAT(' ',I8, I10, 4E15.6)
         CLOSE(111)	 
!         OPEN(111,FILE='iter.dat',ACCESS='Append',STATUS='unknown')
        !WRITE(111,16) ita, nIterPcor, derr2, derrStdSt, dfinish-dstart
!         WRITE(*,16) ita, max_nIterPcor, max_derr2, max_derrStdSt, totime
       ! WRITE(*,16) ita, max_nit, max_derrStdSt, totalTime
! 16      FORMAT(' ',I8, I10, 4E15.6)
!         CLOSE(111)	 
         
         DO g=1,nblocks
        !$omp parallel do collapse(3) private (i,j,k)  num_threads(48)
        !$acc parallel loop gang vector default(present) collapse (3)
         DO k = 1, block(g)%nz+2
         DO j = 1, block(g)%ny+2
         DO i = 1, block(g)%nx+2
           block(g)%u(i,j,k) = block(g)%ut(i,j,k)
           block(g)%v(i,j,k) = block(g)%vt(i,j,k)
           block(g)%w(i,j,k) = block(g)%wt(i,j,k)
         END DO
         END DO
         END DO
        !$acc end parallel
        !$omp end parallel do
         END DO
        !!$acc wait
      ! DO g=1,nblocks

      !  DEALLOCATE (block(g)%b)
      ! END DO
!        CALL updateVelocity
!       if ( ita .eq. 258)then
!               call writeOutput1
!               print*,'before last swap in ps'
!               !pause
!       endif
        CALL fineUpdate_bd
     !  CALL writeOutput1
     !  !pause
        CALL coarseUpdate
!       if ( ita .eq. 258)then
!               call writeOutput1
!               print*,'after last swap in ps'
!               !pause
!       endif
!      CALL writeOutput1
!      print*,ita,'after correct'
!      !pause
      END SUBROUTINE poissonSolver
!***********************************************************************
      SUBROUTINE computeDiv(g)
         USE global
         IMPLICIT NONE
         INTEGER :: n, i, j, k,gg, counter, nx_var, ny_var, ip
         INTEGER(KIND=8),INTENT(IN) ::g
         gg=g
         nx_var=block(g)%nx
         ny_var=block(g)%ny
        !!$acc parallel loop gang vector private (i, j, k)   &
        !!$acc present (block(g)%fluidIndexPtr, b, ut, vt, wt, block(g)%deltax, block(g)%deltay, block(g)%deltaz)         
        !$omp parallel do private (i,j,k,counter)  num_threads(48)
         !DO g=1,nblocks
        !$acc parallel loop gang vector private (i, j, k,counter)   &
        !$acc default(present)         
         DO n = 1, block(gg)%fluidCellCount
           i = block(gg)%fluidIndexPtr(n, 1)
           j = block(gg)%fluidIndexPtr(n, 2)  
           k = block(gg)%fluidIndexPtr(n, 3)  	
           counter =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)           	   
           block(gg)%b(i,j,k) =  (block(gg)%ut(i,j,k) - block(gg)%ut(i-1,j,k))/block(gg)%deltax(i) +  &
                       (block(gg)%vt(i,j,k) - block(gg)%vt(i,j-1,k))/block(gg)%deltay(j) +  &
                       (block(gg)%wt(i,j,k) - block(gg)%wt(i,j,k-1))/block(gg)%deltaz(k)

          !block(gg)%rhs(counter)=block(gg)%b(i,j,k)/deltat

         END DO
        !$acc end parallel
        !$omp end parallel do


       ! if (g .eq. 1)then
       ! DO k=2,block(g)%nz+1
       ! DO i=2,block(g)%nx+1

       !         j=block(g)%ny+1
       !         ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
       !         block(g)%rhs(ip)=block(g)%rhs(ip)-block(g)%An(ip,6)*block(g)%pc(i,j+1,k)

       ! END DO
       ! END DO
       ! end if

       ! if (g .eq. 2)then
       ! DO k=2,block(g)%nz+1
       ! DO i=2,block(g)%nx+1

       !         j=2
       !         ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
       !         block(g)%rhs(ip)=block(g)%rhs(ip)-block(g)%An(ip,2)*block(g)%pc(i,j-1,k)

       ! END DO
       ! END DO
       ! end if

       ! if (g .eq. 3)then
       ! DO k=2,block(g)%nz+1
       ! DO i=2,block(g)%nx+1

       !         j=block(g)%ny+1
       !         ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
       !         block(g)%rhs(ip)=block(g)%rhs(ip)-block(g)%An(ip,6)*block(g)%pc(i,j+1,k)

       ! END DO
       ! END DO
       ! end if




      END SUBROUTINE computeDiv

!***********************************************************************

      SUBROUTINE correctPressure(g)
         USE global
         IMPLICIT NONE
         INTEGER(KIND=8) :: n, i, j, k,gg
         REAL :: r1p, r2p
         INTEGER(KIND=8),INTENT(IN) ::g
        gg=g
        !!$acc parallel loop gang vector private (i, j, k)   &
        !!$acc present (block(g)%fluidIndexPtr, p, pc)
        !$omp parallel do private (i,j,k)  num_threads(48)
        !$acc parallel loop gang vector private (i, j, k)   &
        !$acc default(present)
         DO n = 1, block(gg)%fluidCellCount 
            i = block(gg)%fluidIndexPtr(n, 1)
            j = block(gg)%fluidIndexPtr(n, 2) 
            k = block(gg)%fluidIndexPtr(n, 3) 			
            block(gg)%p(i,j,k) = block(gg)%p(i,j,k) + block(gg)%pc(i,j,k) 
        END DO
        !$omp end parallel do
        !$acc end parallel   
      END SUBROUTINE correctPressure


!***********************************************************************

      SUBROUTINE correctVelocity(g)
         USE global
         IMPLICIT NONE
         INTEGER(KIND=8) :: n, i, j, k,gg
         INTEGER(KIND=8),INTENT(IN) ::g
         gg=g
        !!$acc parallel loop gang vector private (i, j, k) firstprivate (deltat) &
        !!$acc present (block(g)%fluidIndexPtr, ut, vt, wt, block(g)%deltax, block(g)%deltay, block(g)%deltaz, pc)
        !$omp parallel do private (i,j,k) firstprivate(deltat) num_threads(48)
        !$acc parallel loop gang vector private (i, j, k) firstprivate (deltat) &
        !$acc default(present)
         DO 30 n = 1, block(gg)%fluidCellCount 
            i = block(gg)%fluidIndexPtr(n, 1)
            j = block(gg)%fluidIndexPtr(n, 2) 
            k = block(gg)%fluidIndexPtr(n, 3) 			
			
            block(gg)%ut(i,j,k) = block(gg)%ut(i,j,k) - deltat/(0.5d0*(block(gg)%deltax(i+1)+block(gg)%deltax(i)))*(block(gg)%pc(i+1,j,k)-block(gg)%pc(i,j,k))
            block(gg)%vt(i,j,k) = block(gg)%vt(i,j,k) - deltat/(0.5d0*(block(gg)%deltay(j+1)+block(gg)%deltay(j)))*(block(gg)%pc(i,j+1,k)-block(gg)%pc(i,j,k))  
            block(gg)%wt(i,j,k) = block(gg)%wt(i,j,k) - deltat/(0.5d0*(block(gg)%deltaz(k+1)+block(gg)%deltaz(k)))*(block(gg)%pc(i,j,k+1)-block(gg)%pc(i,j,k))  			
 30      CONTINUE
         !$acc end parallel
        !$omp end parallel do
      END SUBROUTINE correctVelocity
      
!***********************************************************************

!***********************************************************************
          
      !SUBROUTINE REDBLACKSOR(epsi, isum, derr, derr2)      
      SUBROUTINE REDBLACKSOR_old(g)      
         USE global
         IMPLICIT NONE
         INTEGER, PARAMETER :: rk = selected_real_kind(8)   
         INTEGER(KIND=8) :: n, i, j, k, gg, ip, nx_var, ny_var, nz_var
         !REAL (KIND = 8) :: derr, derr2,omega, derr3, errSum,var,derr4
         REAL (KIND = 8) :: derr, derr2, derr3, errSum,var,derr4
         !REAL (KIND = 8), INTENT(IN)     :: epsi
         !REAL (KIND = 8), INTENT(OUT)    :: derr, derr2
         !INTEGER (KIND = 8), INTENT(OUT) :: isum
         INTEGER(KIND=8),INTENT(IN) ::g

         gg=g
         !isum = 0   
        ! if (gg .eq. 3)then
        !        omega = 1.67_rk 
        !else 
                !omega=1.96_rk
                omega=1.955_rk
        !endif     
         !block(g)%derr = 0._rk 

          nx_var=block(gg)%nx
          ny_var=block(gg)%ny
          nz_var=block(gg)%nz 

         block(g)%derr2 = 0._rk
         errSum = 0._rk
 3       block(g)%nIterPcor=block(g)%nIterPcor+1 

        derr4=0.
        var=0.
        !isum = isum + 1  
		 
        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k) 
        !$omp parallel do private (i,ip,j,k,n) num_threads(48)
         DO 10 n = 1, block(gg)%redCellCount 
             i = block(gg)%redCellIndexPtr(n, 1)
             j = block(gg)%redCellIndexPtr(n, 2)  
             k = block(gg)%redCellIndexPtr(n, 3)  
!             ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
!     &        -block(gg)%A(ip,3)*block(gg)%pco(i-1,j,k) -block(gg)%A(ip,5)*block(gg)%pco(i+1,j,k) &
!     &        -block(gg)%A(ip,2)*block(gg)%pco(i,j-1,k) -block(gg)%A(ip,6)*block(gg)%pco(i,j+1,k) &
!     &        -block(gg)%A(ip,1)*block(gg)%pco(i,j,k-1) -block(gg)%A(ip,7)*block(gg)%pco(i,j,k+1))/(block(gg)%A(ip,4)) 
!            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) + omega*block(gg)%pc(i,j,k)

            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
     &        -block(gg)%Acx(i-1,1)*block(gg)%pco(i-1,j,k) -block(gg)%Acx(i-1,3)*block(gg)%pco(i+1,j,k) &
     &        -block(gg)%Acy(j-1,1)*block(gg)%pco(i,j-1,k) -block(gg)%Acy(j-1,3)*block(gg)%pco(i,j+1,k) &
     &        -block(gg)%Acz(k-1,1)*block(gg)%pco(i,j,k-1) -block(gg)%Acz(k-1,3)*block(gg)%pco(i,j,k+1))/(block(gg)%Acx(i-1,2)+block(gg)%Acy(j-1,2)+block(gg)%Acz(k-1,2)) 
            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) + omega*block(gg)%pc(i,j,k)

 10      CONTINUE 
        !$omp end parallel do
        !$acc end parallel  
		 
        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k)
        !$omp parallel do private (i,ip,j,k,n) num_threads(48)
         DO 20 n = 1, block(gg)%blackCellCount   
             i = block(gg)%blackCellIndexPtr(n, 1)
             j = block(gg)%blackCellIndexPtr(n, 2)  
             k = block(gg)%blackCellIndexPtr(n, 3) 

!             ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
!     &        -block(gg)%A(ip,3)*block(gg)%pc(i-1,j,k) -block(gg)%A(ip,5)*block(gg)%pc(i+1,j,k) &
!     &        -block(gg)%A(ip,2)*block(gg)%pc(i,j-1,k) -block(gg)%A(ip,6)*block(gg)%pc(i,j+1,k) &
!     &        -block(gg)%A(ip,1)*block(gg)%pc(i,j,k-1) -block(gg)%A(ip,7)*block(gg)%pc(i,j,k+1))/(block(gg)%A(ip,4)) 
!            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) + omega*block(gg)%pc(i,j,k)


            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
     &        -block(gg)%Acx(i-1,1)*block(gg)%pc(i-1,j,k) -block(gg)%Acx(i-1,3)*block(gg)%pc(i+1,j,k) &
     &        -block(gg)%Acy(j-1,1)*block(gg)%pc(i,j-1,k) -block(gg)%Acy(j-1,3)*block(gg)%pc(i,j+1,k) &
     &        -block(gg)%Acz(k-1,1)*block(gg)%pc(i,j,k-1) -block(gg)%Acz(k-1,3)*block(gg)%pc(i,j,k+1))/(block(gg)%Acx(i-1,2)+block(gg)%Acy(j-1,2)+block(gg)%Acz(k-1,2)) 
            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) + omega*block(gg)%pc(i,j,k)

 20      CONTINUE          
        !$omp end parallel do
        !$acc end parallel  
              
        !$acc parallel loop gang vector reduction(max:derr2) default(present) private (i, j, k) 
        !$omp parallel do private (i,j,k,n,var) reduction(max:derr4) num_threads(48)
         DO 30 n = 1, block(gg)%fluidCellCount   
             i = block(gg)%fluidIndexPtr(n, 1)
             j = block(gg)%fluidIndexPtr(n, 2) 
             k = block(gg)%fluidIndexPtr(n, 3) 			 
            var = abs(block(gg)%pc(i,j,k)-block(gg)%pco(i,j,k))     
            derr4=dmax1(derr4,var)
            !errSum = errSum + (block(g)%pc(i,j,k)-block(g)%pco(i,j,k))**2       
           ! block(gg)%pco(i,j,k) = block(gg)%pc(i,j,k)  
 30      CONTINUE 
        !$omp end parallel do
        !$acc end parallel
	
        !$acc parallel loop gang vector collapse(3) default(present) private (i, j, k) 
        !$omp parallel do collapse (3) private (i,j,k) num_threads(48)
        DO k = 1, block(gg)%nz+2
        DO j = 1, block(gg)%ny+2
        DO i = 1, block(gg)%nx+2
            block(gg)%pco(i,j,k) = block(gg)%pc(i,j,k)  
           
        END DO
        END DO
        END DO
        !$omp end parallel do
        !$acc end parallel
        
         block(g)%derr2=derr4
         !derr = dsqrt(errSum/block(g)%fluidCellCount) 
         !derr3 = dmin1(derr, derr2)
         !WRITE(*,*) g,block(g)%nIterPcor, derr4
         !IF (mod(isum,2000).EQ.0)WRITE(*,*) isum, derr, derr2
         !IF (ita.LE.2.AND.isum.LT.50000) GOTO 3
         !IF (derr.gt.epsi) GOTO 3
         IF (derr4.GE.epsi) GOTO 3
         !IF (ita.lt.15.AND.isum.lt.100) GOTO 3         
      END SUBROUTINE REDBLACKSOR_old
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss   
!***********************************************************************

      !SUBROUTINE REDBLACKSOR(epsi, isum, derr, derr2)
      SUBROUTINE REDBLACKSOR(g)
         USE global
         IMPLICIT NONE
         INTEGER, PARAMETER :: rk = selected_real_kind(8)
         INTEGER(KIND=8) :: n, i, j, k, gg, ip, nx_var, ny_var,nz_var,nxy
         REAL (KIND = 8) :: derr, derr2, derr3, errSum,var,derr4
         !REAL (KIND = 8) :: derr, derr2,omega, derr3, errSum,var,derr4
         !REAL (KIND = 8), INTENT(IN)     :: epsi
         !REAL (KIND = 8), INTENT(OUT)    :: derr, derr2
         !INTEGER (KIND = 8), INTENT(OUT) :: isum
         INTEGER(KIND=8),INTENT(IN) ::g

         gg=g
                ! omega = 1.955
            if (gg .eq. 1)then
                 !omega = 1.98_rk
                 omega = omega1
           elseif (gg .eq. 2) then
                !omega=1.96_rk
                 omega=omega2
        !  elseif (gg .eq. 3) then
        !       !omega=1.96_rk
        !        omega=omega3
        !  else         
        !       omega =omega4
          endif

          nx_var=block(gg)%nx
          ny_var=block(gg)%ny
          nz_var=block(gg)%nz

        derr4=11111.
         nxy= nx_var * ny_var
         block(g)%derr2 = 0._rk
         errSum = 0._rk
 3       block(g)%nIterPcor=block(g)%nIterPcor+1

        var=0.
        !isum = isum + 1

        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k)
        !$omp parallel do private (i,ip,j,k,n) num_threads(48)
         DO 10 n = 1, block(gg)%redCellCount
            i = block(gg)%redCellIndexPtr(n, 1)
            j = block(gg)%redCellIndexPtr(n, 2)
            k = block(gg)%redCellIndexPtr(n, 3)

            ip= nx_var*ny_var*(k-2)+ (j-2)*nx_var + (i-1)
          !ip = block(gg)%redCellIndexPtr(n, 1)
          ! k = (ip -1) /(nxy) +2
          ! j = (ip -(k-2)*nxy -1)/nx_var +2
          ! i = (ip -(k-2)*nxy -(j-2)*nx_var) +1

!             ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
!     &        -block(gg)%A(ip,3)*block(gg)%pco(i-1,j,k)
!     -block(gg)%A(ip,5)*block(gg)%pco(i+1,j,k) &
!     &        -block(gg)%A(ip,2)*block(gg)%pco(i,j-1,k)
!     -block(gg)%A(ip,6)*block(gg)%pco(i,j+1,k) &
!     &        -block(gg)%A(ip,1)*block(gg)%pco(i,j,k-1)
!     -block(gg)%A(ip,7)*block(gg)%pco(i,j,k+1))/(block(gg)%A(ip,4))
!            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +
!            omega*block(gg)%pc(i,j,k)

            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
     &        -block(gg)%Ac(ip,3)*block(gg)%pco(i-1,j,k)-block(gg)%Ac(ip ,5)*block(gg)%pco(i+1,j,k) &
     &        -block(gg)%Ac(ip,2)*block(gg)%pco(i,j-1,k)-block(gg)%Ac(ip ,6)*block(gg)%pco(i,j+1,k) &
     &        -block(gg)%Ac(ip,1)*block(gg)%pco(i,j,k-1)-block(gg)%Ac(ip ,7)*block(gg)%pco(i,j,k+1))/(block(gg)%Ac(ip,4))
            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +omega *block(gg)%pc(i,j,k)

 10      CONTINUE
        !$omp end parallel do
        !$acc end parallel

        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k)
        !$omp parallel do private (i,ip,j,k,n) num_threads(48)
         DO 20 n = 1, block(gg)%blackCellCount
          !ip = block(gg)%blackCellIndexPtr(n, 1)
          ! k = (ip -1) /(nxy) +2
          ! j = (ip -(k-2)*nxy -1)/nx_var +2
          ! i = (ip -(k-2)*nxy -(j-2)*nx_var) +1
             i = block(gg)%blackCellIndexPtr(n, 1)
             j = block(gg)%blackCellIndexPtr(n, 2)
             k = block(gg)%blackCellIndexPtr(n, 3)

            ip= nx_var*ny_var*(k-2)+ (j-2)*nx_var + (i-1)
!             ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
!     &        -block(gg)%A(ip,3)*block(gg)%pc(i-1,j,k)
!     -block(gg)%A(ip,5)*block(gg)%pc(i+1,j,k) &
!     &        -block(gg)%A(ip,2)*block(gg)%pc(i,j-1,k)
!     -block(gg)%A(ip,6)*block(gg)%pc(i,j+1,k) &
!     &        -block(gg)%A(ip,1)*block(gg)%pc(i,j,k-1)
!     -block(gg)%A(ip,7)*block(gg)%pc(i,j,k+1))/(block(gg)%A(ip,4))
!            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +
!            omega*block(gg)%pc(i,j,k)


            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
     &        -block(gg)%Ac(ip ,3)*block(gg)%pc(i-1,j,k)-block(gg)%Ac(ip,5)*block(gg)%pc(i+1,j,k) &
     &        -block(gg)%Ac(ip ,2)*block(gg)%pc(i,j-1,k)-block(gg)%Ac(ip,6)*block(gg)%pc(i,j+1,k) &
     &        -block(gg)%Ac(ip ,1)*block(gg)%pc(i,j,k-1)-block(gg)%Ac(ip,7)*block(gg)%pc(i,j,k+1))/(block(gg)%Ac(ip,4))
            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +omega*block(gg)%pc(i,j,k)

 20      CONTINUE
        !$omp end parallel do
        !$acc end parallel


       if(mod(block(gg)%nIterPcor,100) .eq.0)then
       derr4=0.
        !$acc parallel loop gang vector reduction(max:derr4) default(present) private (i, j, k, var)
        !$omp parallel do private (i,j,k,n,var) reduction(max:derr4) num_threads(48)
         DO 30 n = 1, block(gg)%fluidCellCount
             i = block(gg)%fluidIndexPtr(n, 1)
             j = block(gg)%fluidIndexPtr(n, 2)
             k = block(gg)%fluidIndexPtr(n, 3)
            var = abs(block(gg)%pc(i,j,k)-block(gg)%pco(i,j,k))
            derr4=dmax1(derr4,var)
            !errSum = errSum +
            !(block(g)%pc(i,j,k)-block(g)%pco(i,j,k))**2
           ! block(gg)%pco(i,j,k) = block(gg)%pc(i,j,k)
 30      CONTINUE
        !$omp end parallel do
        !$acc end parallel
        end if
        !$acc parallel loop gang vector collapse(3) default(present) private (i, j, k)
        !$omp parallel do collapse (3) private (i,j,k) num_threads(48)
        DO k = 1, block(gg)%nz+2
        DO j = 1, block(gg)%ny+2
        DO i = 1, block(gg)%nx+2
            block(gg)%pco(i,j,k) = block(gg)%pc(i,j,k)

        END DO
        END DO
        END DO
        !$omp end parallel do
        !$acc end parallel
      ! if( gg.eq. 1)then
      !  print*,block(gg)%nIterPcor,ita,block(1)%pc(14,248,15)
      !  endif
         block(g)%derr2=derr4
         !derr = dsqrt(errSum/block(g)%fluidCellCount)
         !derr3 = dmin1(derr, derr2)
         !WRITE(*,*) g,block(g)%nIterPcor, derr4
         !IF (mod(isum,2000).EQ.0)WRITE(*,*) isum, derr, derr2
         !IF (ita.LE.2.AND.isum.LT.50000) GOTO 3
         !IF (derr.gt.epsi) GOTO 3
         IF (derr4.GE.epsi .and. block(g)%nIterPcor .le. pcItaMax) GOTO 3
         !IF (ita.lt.15.AND.isum.lt.100) GOTO 3
      END SUBROUTINE REDBLACKSOR
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!***********************************************************************

      !SUBROUTINE REDBLACKSOR(epsi, isum, derr, derr2)
      SUBROUTINE REDBLACKSOR_linear(g)
         USE global
         IMPLICIT NONE
         INTEGER, PARAMETER :: rk = selected_real_kind(8)
         INTEGER(KIND=8) :: n, i, j, k, gg, ip, nx_var, ny_var,nz_var,nxy
         REAL (KIND = 8) :: derr, derr2, derr3, errSum,var,derr4
         !REAL (KIND = 8) :: derr, derr2,omega, derr3, errSum,var,derr4
         !REAL (KIND = 8), INTENT(IN)     :: epsi
         !REAL (KIND = 8), INTENT(OUT)    :: derr, derr2
         !INTEGER (KIND = 8), INTENT(OUT) :: isum
         INTEGER(KIND=8),INTENT(IN) ::g

         gg=g
                ! omega = 1.955
            if (gg .eq. 1)then
                 !omega = 1.98_rk
                 omega = omega1
           elseif (gg .eq. 2) then
                !omega=1.96_rk
                 omega=omega2
           elseif (gg .eq. 3) then
                !omega=1.96_rk
                 omega=omega3
           else         
                omega =omega4
          endif

          nx_var=block(gg)%nx
          ny_var=block(gg)%ny
          nz_var=block(gg)%nz

        derr4=11111.
         nxy= nx_var * ny_var
         block(g)%derr2 = 0._rk
         errSum = 0._rk
 3       block(g)%nIterPcor=block(g)%nIterPcor+1

        var=0.
        !isum = isum + 1

        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k)
        !$omp parallel do private (i,ip,j,k,n) num_threads(48)
         DO 10 n = 1, block(gg)%redCellCount
            i = block(gg)%redCellIndexPtr(n, 1)
            j = block(gg)%redCellIndexPtr(n, 2)
            k = block(gg)%redCellIndexPtr(n, 3)
          !ip = block(gg)%redCellIndexPtr(n, 1)
          ! k = (ip -1) /(nxy) +2
          ! j = (ip -(k-2)*nxy -1)/nx_var +2
          ! i = (ip -(k-2)*nxy -(j-2)*nx_var) +1

!             ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
!     &        -block(gg)%A(ip,3)*block(gg)%pco(i-1,j,k)
!     -block(gg)%A(ip,5)*block(gg)%pco(i+1,j,k) &
!     &        -block(gg)%A(ip,2)*block(gg)%pco(i,j-1,k)
!     -block(gg)%A(ip,6)*block(gg)%pco(i,j+1,k) &
!     &        -block(gg)%A(ip,1)*block(gg)%pco(i,j,k-1)
!     -block(gg)%A(ip,7)*block(gg)%pco(i,j,k+1))/(block(gg)%A(ip,4))
!            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +
!            omega*block(gg)%pc(i,j,k)

            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
     &        -block(gg)%Acx(i-1,1)*block(gg)%pco(i-1,j,k)-block(gg)%Acx(i-1,3)*block(gg)%pco(i+1,j,k) &
     &        -block(gg)%Acy(j-1,1)*block(gg)%pco(i,j-1,k)-block(gg)%Acy(j-1,3)*block(gg)%pco(i,j+1,k) &
     &        -block(gg)%Acz(k-1,1)*block(gg)%pco(i,j,k-1)-block(gg)%Acz(k-1,3)*block(gg)%pco(i,j,k+1))/(block(gg)%Acx(i-1,2)+block(gg)%Acy(j-1,2)+block(gg)%Acz(k-1,2))
            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +omega*block(gg)%pc(i,j,k)

 10      CONTINUE
        !$omp end parallel do
        !$acc end parallel

        !$acc parallel loop gang vector default(present) firstprivate(deltat, omega) private (i, j, k)
        !$omp parallel do private (i,ip,j,k,n) num_threads(48)
         DO 20 n = 1, block(gg)%blackCellCount
          !ip = block(gg)%blackCellIndexPtr(n, 1)
          ! k = (ip -1) /(nxy) +2
          ! j = (ip -(k-2)*nxy -1)/nx_var +2
          ! i = (ip -(k-2)*nxy -(j-2)*nx_var) +1
             i = block(gg)%blackCellIndexPtr(n, 1)
             j = block(gg)%blackCellIndexPtr(n, 2)
             k = block(gg)%blackCellIndexPtr(n, 3)

!             ip =i-1  + nx_var*(j-2)  + nx_var*ny_var*(k-2)
!            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
!     &        -block(gg)%A(ip,3)*block(gg)%pc(i-1,j,k)
!     -block(gg)%A(ip,5)*block(gg)%pc(i+1,j,k) &
!     &        -block(gg)%A(ip,2)*block(gg)%pc(i,j-1,k)
!     -block(gg)%A(ip,6)*block(gg)%pc(i,j+1,k) &
!     &        -block(gg)%A(ip,1)*block(gg)%pc(i,j,k-1)
!     -block(gg)%A(ip,7)*block(gg)%pc(i,j,k+1))/(block(gg)%A(ip,4))
!            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +
!            omega*block(gg)%pc(i,j,k)


            block(gg)%pc(i,j,k) = (block(gg)%b(i,j,k)/deltat &
     &        -block(gg)%Acx(i-1,1)*block(gg)%pc(i-1,j,k)-block(gg)%Acx(i-1,3)*block(gg)%pc(i+1,j,k) &
     &        -block(gg)%Acy(j-1,1)*block(gg)%pc(i,j-1,k)-block(gg)%Acy(j-1,3)*block(gg)%pc(i,j+1,k) &
     &        -block(gg)%Acz(k-1,1)*block(gg)%pc(i,j,k-1)-block(gg)%Acz(k-1,3)*block(gg)%pc(i,j,k+1))/(block(gg)%Acx(i-1,2)+block(gg)%Acy(j-1,2)+block(gg)%Acz(k-1,2))
            block(gg)%pc(i,j,k) = (1._rk-omega)*block(gg)%pco(i,j,k) +omega*block(gg)%pc(i,j,k)

 20      CONTINUE
        !$omp end parallel do
        !$acc end parallel


       if(mod(block(gg)%nIterPcor,5) .eq.0)then
       derr4=0.
        !$acc parallel loop gang vector reduction(max:derr4) default(present) private (i, j, k, var)
        !$omp parallel do private (i,j,k,n,var) reduction(max:derr4) num_threads(48)
         DO 30 n = 1, block(gg)%fluidCellCount
             i = block(gg)%fluidIndexPtr(n, 1)
             j = block(gg)%fluidIndexPtr(n, 2)
             k = block(gg)%fluidIndexPtr(n, 3)
            var = abs(block(gg)%pc(i,j,k)-block(gg)%pco(i,j,k))
            derr4=dmax1(derr4,var)
            !errSum = errSum +
            !(block(g)%pc(i,j,k)-block(g)%pco(i,j,k))**2
           ! block(gg)%pco(i,j,k) = block(gg)%pc(i,j,k)
 30      CONTINUE
        !$omp end parallel do
        !$acc end parallel
        end if
        !$acc parallel loop gang vector collapse(3) default(present) private (i, j, k)
        !$omp parallel do collapse (3) private (i,j,k) num_threads(48)
        DO k = 1, block(gg)%nz+2
        DO j = 1, block(gg)%ny+2
        DO i = 1, block(gg)%nx+2
            block(gg)%pco(i,j,k) = block(gg)%pc(i,j,k)

        END DO
        END DO
        END DO 
        !$omp end parallel do
        !$acc end parallel

         block(g)%derr2=derr4
         !derr = dsqrt(errSum/block(g)%fluidCellCount)
         !derr3 = dmin1(derr, derr2)
         !WRITE(*,*) g,block(g)%nIterPcor, derr4
         !IF (mod(isum,2000).EQ.0)WRITE(*,*) isum, derr, derr2
         !IF (ita.LE.2.AND.isum.LT.50000) GOTO 3
         !IF (derr.gt.epsi) GOTO 3
         !IF (derr4.GE.epsi) GOTO 3
         IF (derr4.GE.epsi .and. block(g)%nIterPcor .le. pcItaMax) GOTO 3
         !IF (ita.lt.15.AND.isum.lt.100) GOTO 3
      END SUBROUTINE REDBLACKSOR_linear
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!*******************************************************************
      SUBROUTINE updateVelocity_newv(g)
      !SUBROUTINE updateVelocity_newv
        USE global
        IMPLICIT NONE
        INTEGER :: n, i, j, k
         INTEGER (KIND = 8), INTENT(IN) :: g
!        OPEN(UNIT=111, File='velocity.dat',STATUS='unknown')
        !!$acc parallel loop present (u, ut, v, vt)
        !DO g = blk_start, nblocks
!        print*,g,'inside_pcor'
        if (block(g)%move_check .eq. 1) then
        !$acc parallel loop gang vector collapse(2) default(present)
        DO 70 k = 1, block(g)%nz+2
        DO 70 j = 1, block(g)%ny+2
        DO 70 i = 1, block(g)%nx+2
       !!do 70 j=1,block(g)%jtc_en
       !!do 70 i=1,block(g)%itc_en
           block(g)%ut(i,j,k) = block(g)%u(i,j,k)
           block(g)%vt(i,j,k) = block(g)%v(i,j,k)
           block(g)%wt(i,j,k) = block(g)%w(i,j,k)
!         WRITE(111,*)g,i,j, block(g)%u(i,j), block(g)%v(i,j)
!111          FORMAT('',3I4,2F12.4)
 70    CONTINUE 
       !$acc end parallel
        endif
       ! END DO
     
      END SUBROUTINE updateVelocity_newv

!********************************************************************
              







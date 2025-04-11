module biocfd_write_output_corner1
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  ! allow(use-all) - TODO: Aim to fix this in the future
  use global
  implicit none

  private

  public :: writeOutput1, body_plot, writeresult

contains
      SUBROUTINE writeOutput1
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       CHARACTER(len=150)  :: filename1
       INTEGER  :: k, i, j, g
       REAL (dp) :: u1, v1, w1

         IF((mod(ita,200) ==0 .or. ita <= 2 ))then

           Do g=1,nblocks
       WRITE(filename1,1)char_f,ita,g,re,block(2)%dx,nblocks
1     FORMAT('out/',A3,'_butter_fielddata.',i9.9,'.',i3.3,'.',f7.1,'.',f8.6,'.',i3.3,".dat")
           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables="x","y","z","u","v","w","p","totime","cellid","cell_n","cell_pr"'
            WRITE(786,*) 'zone, ', 'i = ', block(g)%nx,' j = ', block(g)%ny, ' k = ', block(g)%nz

            DO 30 k = 2, block(g)%nz+1
            DO 30 j = 2, block(g)%ny+1
            DO 30 i = 2, block(g)%nx+1
               u1 = 0.5_dp*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))
               v1 = 0.5_dp*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))
               w1 = 0.5_dp*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))
               WRITE(786,*) block(g)%xp(i), block(g)%yp(j), block(g)%zp(k), u1, v1, w1, &
                            block(g)%p(i,j,k), totime, block(g)%cell(i,j,k) , &
                            block(g)%cell_n(i,j,k) , block(g)%cell_pr(i,j,k)
 30         CONTINUE
            CLOSE(786)
        end do
        !$acc wait
         ENDIF
      END SUBROUTINE writeOutput1

      SUBROUTINE writeResult
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k,g
       CHARACTER(len=70)  :: filename1
        IF(mod(ita,500)==0)THEN
           Do g=1,nblocks
           WRITE(filename1,22)char_f,g,re,block(2)%dx
 22          FORMAT('out/Chkpt/',A3,'_butter_chkpt.',i3.3,'.',f6.1,'.',f8.6,".dat")
        OPEN (1,FILE=filename1,FORM='formatted')
        DO 30 k = 1, block(g)%nz+2
        DO 30 j = 1, block(g)%ny+2
        DO 30 i = 1, block(g)%nx+2
          WRITE(1,*) block(g)%u(i,j,k), block(g)%v(i,j,k), block(g)%w(i,j,k), &
        block(g)%p(i,j,k), totime, ita, ita1
 30     CONTINUE
        CLOSE(22)
        END DO
        END IF
      END SUBROUTINE writeResult

         SUBROUTINE body_plot
         INTEGER(int64) :: inode, ielem, g
         CHARACTER(len=150) :: filename1

          DO g=1,nblocks
          if (ita == 1 )then
          !DO g=1,nblocks
          WRITE(filename1,108)
  108     FORMAT("out/butterfly.dat")
          OPEN(UNIT=857,FILE=filename1,STATUS='unknown')
          WRITE(857,*) 'TITLE = "FEstressplot"'
          WRITE(857,*) 'VARIABLES= "x", "y", "z","bd_n"'
          WRITE(857,*) 'ZONE NODES= ',block(g)%ibNodes,',ELEMENTS=',block(g)%ibElems,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
          DO inode = 1, block(g)%ibNodes
            WRITE(857,*)block(g)%xnode1(inode),block(g)%ynode1(inode),block(g)%znode1(inode), 99
          END DO
          DO ielem = 1, block(g)%ibElems
            WRITE(857,*) block(g)%ibElP1(ielem), block(g)%ibElP2(ielem), block(g)%ibElP3(ielem)
          END DO
          CLOSE(857)
         END IF
          END DO
        end subroutine body_plot

end module biocfd_write_output_corner1

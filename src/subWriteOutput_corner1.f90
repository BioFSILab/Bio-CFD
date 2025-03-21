!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!     SUBROUTINE writeOutput
!      USE global
!      INTEGER, PARAMETER :: rk = selected_real_kind(8)
!      CHARACTER*120  filename1
!      CHARACTER*70  filename2
!      INTEGER  :: k, i, j  ,g
!      REAL (KIND = 8) :: u1, v1, w1

!       IF(mod(ita,100).EQ.0)THEN
!       DO g=1,nblocks
!           WRITE(filename1,1) ita,g
!1          FORMAT('out/fielddata.',i9.9,'.',i3.3,".dat")
!           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
!           WRITE(786,*)'variables= "x","y","z","u","v","w","p","totime","cellid"'
!           WRITE(786,*) 'zone, ', 'i = ', block(g)%nx+1,' j = ', block(g)%ny+1, ' k = ', block(g)%nz+1
!           !k = block(g)%nz/2
!           DO 30 k = 1, block(g)%nz+1
!           DO 30 j = 1, block(g)%ny+1
!           DO 30 i = 1, block(g)%nx+1

!            bl_intp_y1= block(g)%yp(j)
!            bl_intp_y2= block(g)%yp(j+1)
!            bl_intp_x1= block(g)%xp(i)
!            bl_intp_x2= block(g)%xp(i+1)
!            bl_intp_f1= block(g)%p(i,j,k)
!            bl_intp_f2= block(g)%p(i+1,j,k)
!            bl_intp_f3= block(g)%p(i+1,j+1,k)
!            bl_intp_f4= block(g)%p(i,j+1,k)
!            bl_intp_valy=block(g)%y1(j+1)
!            bl_intp_valx=block(g)%x1(i+1)
!          call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)

!            p_new1=bl_interp_ans


!            bl_intp_y1= block(g)%yp(j)
!            bl_intp_y2= block(g)%yp(j+1)
!            bl_intp_x1= block(g)%xp(i)
!            bl_intp_x2= block(g)%xp(i+1)
!            bl_intp_f1= block(g)%p(i,j,k+1)
!            bl_intp_f2= block(g)%p(i+1,j,k+1)
!            bl_intp_f3= block(g)%p(i+1,j+1,k+1)
!            bl_intp_f4= block(g)%p(i,j+1,k+1)
!            bl_intp_valy=block(g)%y1(j+1)
!            bl_intp_valx=block(g)%x1(i+1)
!          call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)

!            p_new2=bl_interp_ans


!
!            l_intp_valx=block(g)%z1(k+1)
!            l_intp_x1= 0.5*(block(g)%z1(k)+block(g)%z1(k+1))
!            l_intp_x2= 0.5*(block(g)%z1(k+1)+block(g)%z1(k+2))
!            l_intp_y1=p_new1
!            l_intp_y2=p_new2
!            call linearInterp(l_intp_valx,l_intp_x1,l_intp_x2,l_intp_y1,l_intp_y2,l_intp_valy)
!            p_final=l_intp_valy


!            bl_intp_y1= 0.5*(block(g)%z1(k)+block(g)%z1(k+1))
!            bl_intp_y2= 0.5*(block(g)%z1(k+1)+block(g)%z1(k+2))
!            bl_intp_x1= 0.5*(block(g)%x1(i)+block(g)%x1(i+1))
!            bl_intp_x2= 0.5*(block(g)%x1(i+1)+block(g)%x1(i+2))
!            bl_intp_f1= block(g)%v(i,j,k)
!            bl_intp_f2= block(g)%v(i+1,j,k)
!            bl_intp_f3= block(g)%v(i+1,j,k+1)
!            bl_intp_f4= block(g)%v(i,j,k+1)
!            bl_intp_valy=block(g)%z1(k+1)
!            bl_intp_valx=block(g)%x1(i+1)
!          call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)

!            v_new=bl_interp_ans

!            bl_intp_y1= 0.5*(block(g)%z1(k)+block(g)%z1(k+1))
!            bl_intp_y2= 0.5*(block(g)%z1(k+1)+block(g)%z1(k+2))
!            bl_intp_x1= 0.5*(block(g)%y1(j)+block(g)%y1(j+1))
!            bl_intp_x2= 0.5*(block(g)%y1(j+1)+block(g)%y1(j+2))
!            bl_intp_f1= block(g)%u(i,j,k)
!            bl_intp_f2= block(g)%u(i,j,k+1)
!            bl_intp_f3= block(g)%u(i,j+1,k+1)
!            bl_intp_f4= block(g)%u(i,j+1,k)
!            bl_intp_valy=block(g)%z1(k+1)
!            bl_intp_valx=block(g)%y1(j+1)
!          call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)

!            u_new=bl_interp_ans

!            bl_intp_y1= 0.5*(block(g)%y1(j)+block(g)%y1(j+1))
!            bl_intp_y2= 0.5*(block(g)%y1(j+1)+block(g)%y1(j+2))
!            bl_intp_x1= 0.5*(block(g)%x1(i)+block(g)%x1(i+1))
!            bl_intp_x2= 0.5*(block(g)%x1(i+1)+block(g)%x1(i+2))
!            bl_intp_f1= block(g)%w(i,j,k)
!            bl_intp_f2= block(g)%w(i+1,j,k)
!            bl_intp_f3= block(g)%w(i+1,j+1,k)
!            bl_intp_f4= block(g)%w(i,j+1,k)
!            bl_intp_valy=block(g)%y1(j+1)
!            bl_intp_valx=block(g)%x1(i+1)
!          call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)

!            w_new=bl_interp_ans



!              WRITE(786,*) block(g)%x1(i+1), block(g)%y1(j+1), block(g)%z1(k+1), u_new, v_new, w_new, p_final, totime, block(g)%cell(i,j,k)
!30         CONTINUE
!           CLOSE(786)
!       END DO
!        ENDIF
!     END SUBROUTINE writeOutput
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE writeOutput1
       USE global
!!       USE mod_tecio
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       CHARACTER(len=150)  :: filename1
       CHARACTER(len=150)  :: filename2
       INTEGER  :: k, i, j, g
       REAL (KIND = 8) :: u1, v1, w1

!!         CHARACTER*1 NULLCHR
!!       INTEGER*4 iStat, FileFormat, FileType, Debug, VIsDouble
!!       integer*4 ZoneType, IMax, JMax, KMax, ICellMax, JCellMax,KCellMax
!!       integer*4 StrandID, unused, IsBlock, NFConns, FNMode, ShrConn
!!       integer,  dimension(:), pointer :: null
!!       integer*4 nTotalPts
!!       integer*4,dimension(11) :: valueLoc
!!       real*8    SolTime
!!       REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:,:) :: un1, vn1, wn1



         !IF(mod(ita,10000) .eq.0 )THEN
         !IF(mod(ita,1000) .eq.0 .or. mod(ita,3) .eq. 0)THEN
         IF((mod(ita,200) ==0 .or. ita <= 2 ))then  ! .or. &
         !   (ita .gt. 250 .and. mod(ita,1) .eq. 0)) THEN
         !IF(mod(ita,100) .eq.0 )THEN

!!      Do g=1,nblocks
!!
!!       nx_var=block(g)%nx
!!       ny_var=block(g)%ny
!!       nz_var=block(g)%nz
!!      ALLOCATE (un1(nx_var+1,ny_var+1,nz_var+1), vn1(nx_var+1,ny_var+1,nz_var+1),wn1(nx_var+1,ny_var+1,nz_var+1))
!!      WRITE(filename1,1)char_f,ita,g,re,block(2)%dx,nblocks
!!1     FORMAT('out/',A3,'_butter_fielddata.',i9.9,'.',i3.3,'.',f7.1,'.',f8.6,'.',i3.3,".plt")
!!     NULLCHR = CHAR(0)
!!     FileFormat = 0  !0= PLT;  1=SZPLT
!!     FileType = 0    !1= Grid; 2=Field; 0=Grid+Field
!!     Debug = 1
!!     VIsDouble = 0
!!
!!     ZoneType = 0
!!     IMax = block(g)%nx
!!     JMax = block(g)%ny
!!     KMax = block(g)%nz
!!     ICellMax = 0
!!     JCellMax = 0
!!     KCellMax = 0
!!
!!     StrandID = 1
!!     unused = 0
!!     isBlock = 1
!!     NFConns = 0
!!     FNMode = 0
!!     ShrConn = 0
!!     valueLoc(:) = 1
!!     SolTime = totime
!!
!!    !write(filename1,101) ita
!!    !101  format('out/field.',i10.10,'.plt')
!!
!!     iStat = tecIni142('SweptAngleStudy'//NULLCHR,  &
!!                       'x y z p  cell cell_n cell_pr u v w b'//NULLCHR,       &
!!                       filename1//NULLCHR,          &
!!                       './out/'//NULLCHR,    &   ! Scratch directory
!!                       FileFormat,                 &
!!                       FileType,                   &
!!                       Debug,                      &
!!                       VIsDouble)
!!
!!     iStat = tecZne142('Instantaneous Plots'//NULLCHR, &
!!                       ZoneType,               &
!!                       IMax,                   &
!!                       JMax,                   &
!!                       KMax,                   &
!!                       ICellMax,               &
!!                       JCellMax,               &
!!                       KCellMax,               &
!!                       SolTime,                &
!!                       StrandID,               &
!!                       unused,                 &
!!                       IsBlock,                &
!!                       NFConns,                &
!!                       FNMode,                 &
!!                       0,                      &
!!                       0,                      &
!!                       0,                      &
!!                       Null,                   &
!!                       valueLoc,               &
!!                       Null,                   &
!!                       ShrConn)
!!
!!     un1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%u(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%u(1:nx_var-1,2:ny_var,2:nz_var)) *0.5
!!     vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%v(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%v(2:nx_var,1:ny_var-1,2:nz_var)) *0.5
!!     wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%w(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%w(2:nx_var,2:ny_var,1:nz_var-1)) *0.5
!!    !un1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%ut(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%ut(1:nx_var-1,2:ny_var,2:nz_var)) *0.5
!!    !vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%vt(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%vt(2:nx_var,1:ny_var-1,2:nz_var)) *0.5
!!    !wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%wt(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%wt(2:nx_var,2:ny_var,1:nz_var-1)) *0.5
!!
!!     nTotalPts = (nx_var)*(ny_var)*(nz_var)
!!    !nTotalPts = (nx_var)*(ny_var)*(1)
!!      iStat = tecDat142(nTotalPts,      block(g)%xp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!      iStat = tecDat142(nTotalPts,      block(g)%yp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!      iStat = tecDat142(nTotalPts,      block(g)%zp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!     !iStat = tecDat142(nTotalPts,     xp(2:nx_var,4)    , 0)
!!     !iStat = tecDat142(nTotalPts,     yp(2:ny_var,4)    , 0)
!!     !iStat = tecDat142(nTotalPts,     zp(2:nz_var,4)    , 0)
!!     iStat = tecDat142(nTotalPts, real(  block(g)%p(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!     !iStat = tecDat142(nTotalPts, real( dmiu(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),     0)
!!     iStat = tecDat142(nTotalPts, real( block(g)%cell(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!     iStat = tecDat142(nTotalPts, real( block(g)%cell_n(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!     iStat = tecDat142(nTotalPts, real( block(g)%cell_pr(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!     iStat = tecDat142(nTotalPts, real(un1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!     iStat = tecDat142(nTotalPts, real(vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!     iStat = tecDat142(nTotalPts, real(wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!     iStat = tecDat142(nTotalPts, real(  block(g)%b(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!     iStat = tecEnd142()
!!
!!       DEALLOCATE(un1,vn1,wn1)
!!      end do

!             WRITE(filename1,1) ita


           Do g=1,nblocks
!!            WRITE(filename1,1) msh_t,geo_num,ita,g
!! 1          FORMAT('out/d',i3.3,'_',i6,'_fielddata.',i9.9,'.',i3.3,".dat")
       WRITE(filename1,1)char_f,ita,g,re,block(2)%dx,nblocks
1     FORMAT('out/',A3,'_butter_fielddata.',i9.9,'.',i3.3,'.',f7.1,'.',f8.6,'.',i3.3,".dat")
           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables="x","y","z","u","v","w","p","totime","cellid","cell_n","cell_pr"'
            WRITE(786,*) 'zone, ', 'i = ', block(g)%nx,' j = ', block(g)%ny, ' k = ', block(g)%nz
            !k = block(g)%nz/2
            DO 30 k = 2, block(g)%nz+1
            DO 30 j = 2, block(g)%ny+1
            DO 30 i = 2, block(g)%nx+1
               u1 = 0.5*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))
               v1 = 0.5*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))
               w1 = 0.5*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))
               WRITE(786,*) block(g)%xp(i), block(g)%yp(j), block(g)%zp(k), u1, v1, w1, block(g)%p(i,j,k), totime, block(g)%cell(i,j,k) , block(g)%cell_n(i,j,k) , block(g)%cell_pr(i,j,k)
 30         CONTINUE
            CLOSE(786)
        end do
        !$acc wait
         ENDIF
      END SUBROUTINE writeOutput1
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

!***********************************************************************

!***********************************************************************
      SUBROUTINE writeResult
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k,g
       CHARACTER(len=70)  :: filename1
        IF(mod(ita,500)==0)THEN
           Do g=1,nblocks
           WRITE(filename1,22)char_f,g,re,block(2)%dx
 !22          FORMAT('out/4blk/aorta_chkpt.',i3.3,'.',f6.1,".dat")
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
!***********************************************************************

!***********************************************************************
     SUBROUTINE writeTagging
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       CHARACTER(len=70)  :: filename1
       CHARACTER(len=70)  :: filename2
       INTEGER  :: k, i, j,g
        g=2
            WRITE(filename1,1) geo_num,ita
 1          FORMAT('out/',i6.6,'/tag/tagdata.',i2.2,".dat")
            OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables= "x","y","z","cellid"'
          WRITE(786,*) 'zone, ', 'i = ', block(g)%nx,' j = ', block(g)%ny, ' k = ', block(g)%nz
            !k = block(g)%nz/2
            DO 30 k = 2, block(g)%nz+1
            DO 30 j = 2, block(g)%ny+1
            DO 30 i = 2, block(g)%nx+1
               WRITE(786,*) block(g)%xp(i), block(g)%yp(j), block(g)%zp(k), block(g)%cell(i,j,k)
 30         CONTINUE
            CLOSE(786)
      END SUBROUTINE writeTagging
!***********************************************************************

!***********************************************************************
      SUBROUTINE computeSumData
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER         :: i, j, k,g
             DO g=1,nblocks
           !$acc parallel loop gang vector collapse (3) default(present)
	     DO 30 k = 1, block(g)%nz+2
             DO 30 j = 1, block(g)%ny+2
             DO 30 i = 1, block(g)%nx+2
               block(g)%u_sum(i,j,k) = block(g)%u_sum(i,j,k) + block(g)%u(i,j,k)
               block(g)%v_sum(i,j,k) = block(g)%v_sum(i,j,k) + block(g)%v(i,j,k)
               block(g)%w_sum(i,j,k) = block(g)%w_sum(i,j,k) + block(g)%w(i,j,k)
               block(g)%p_sum(i,j,k) = block(g)%p_sum(i,j,k) + block(g)%p(i,j,k)
 30         CONTINUE
            !$acc end parallel

            END DO
      END SUBROUTINE computeSumData
!***********************************************************************

!***********************************************************************
      SUBROUTINE writeComputeSumData
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k, g
       CHARACTER(len=70)  :: filename1


        IF(mod(ita1,1000)==0)THEN
        DO g=1,nblocks
        WRITE(filename1,22)g,re
 22     FORMAT('out/4blk/aorta_sumdata.',i3.3,'.',f6.1,".dat")
         OPEN (1,FILE=filename1,FORM='formatted')
 	 DO 30 k = 1, block(g)%nz+2
        DO 30 j = 1, block(g)%ny+2
        DO 30 i = 1, block(g)%nx+2
          WRITE(1,*) block(g)%u_sum(i,j,k), block(g)%v_sum(i,j,k), block(g)%w_sum(i,j,k),block(g)%p_sum(i,j,k), ita, ita1
 30     CONTINUE
 	 CLOSE(1)
        END DO
 	 END IF
      END SUBROUTINE  writeComputeSumData
!***********************************************************************

      SUBROUTINE writeComputeSumSqData
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k, g
       CHARACTER(len=70)  :: filename1


        IF(mod(ita1,1000)==0)THEN
        DO g=1,nblocks
        WRITE(filename1,22)g,re
 22     FORMAT('out/4blk/aorta_sumsqdata.',i3.3,'.',f6.1,".dat")
         OPEN (1,FILE=filename1,FORM='formatted')
 	 DO 30 k = 1, block(g)%nz+2
        DO 30 j = 1, block(g)%ny+2
        DO 30 i = 1, block(g)%nx+2
          WRITE(1,*) block(g)%u2_sum(i,j,k), block(g)%v2_sum(i,j,k), block(g)%w2_sum(i,j,k), block(g)%p2_sum(i,j,k), block(g)%uv_sum(i,j,k), block(g)%vw_sum(i,j,k), block(g)%uw_sum(i,j,k)
 30     CONTINUE
 	 CLOSE(1)
        END DO
 	 END IF
      END SUBROUTINE  writeComputeSumSqData
!***********************************************************************
      SUBROUTINE readComputeSumData
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k, g
       CHARACTER(len=70)  :: filename1

        DO g=1,nblocks
        WRITE(filename1,22)g,re
 22     FORMAT('out/aorta_sumdata.',i3.3,'.',f6.1,".dat")
!        WRITE(filename1,22) g
! 22     FORMAT('sumdata.',i3.3,".dat")
        OPEN (1,FILE=filename1,FORM='formatted')
        DO 30 k = 1, block(g)%nz+2
        DO 30 j = 1, block(g)%ny+2
        DO 30 i = 1, block(g)%nx+2
          READ(1,*) block(g)%u_sum(i,j,k), block(g)%v_sum(i,j,k), block(g)%w_sum(i,j,k), block(g)%p_sum(i,j,k), ita, ita1
 30     CONTINUE
        CLOSE(1)
        END DO
      END SUBROUTINE  readComputeSumData
!***********************************************************************

!***********************************************************************
      SUBROUTINE readComputeSumSqData
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k, g
       CHARACTER(len=70)  :: filename1

        DO g=1,nblocks
        WRITE(filename1,22)g,re
 22     FORMAT('out/aorta_sumsqdata.',i3.3,'.',f6.1,".dat")
!        WRITE(filename1,22) g
! 22     FORMAT('sumdata.',i3.3,".dat")
        OPEN (1,FILE=filename1,FORM='formatted')
        DO 30 k = 1, block(g)%nz+2
        DO 30 j = 1, block(g)%ny+2
        DO 30 i = 1, block(g)%nx+2
        READ(1,*) block(g)%u2_sum(i,j,k), block(g)%v2_sum(i,j,k), block(g)%w2_sum(i,j,k), block(g)%p2_sum(i,j,k), block(g)%uv_sum(i,j,k), block(g)%vw_sum(i,j,k), block(g)%uw_sum(i,j,k)
 30     CONTINUE
        CLOSE(1)
        END DO
      END SUBROUTINE  readComputeSumSqData
!***********************************************************************
!***********************************************************************
      SUBROUTINE readstressData
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  inode, g
        DO g=1,nblocks
        OPEN (1,FILE='stressdata_fl',FORM='formatted')
        DO inode = 1, block(g)%ibNodes
        READ(1,*) SQSUMWSS(inode),SUMWSS(inode), SIGNWSS(inode,1),SIGNWSS(inode,2),&
         SIGNWSS(inode,3), ita1, ita

        END DO
        END DO
        CLOSE(1)
      END SUBROUTINE readstressData
!***********************************************************************

!***********************************************************************
      SUBROUTINE computeAvgData
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER  :: n, k, i, j, g
             DO 30 g=1,nblocks
	     DO 30 k = 1, block(g)%nz+2
             DO 30 j = 1, block(g)%ny+2
             DO 30 i = 1, block(g)%nx+2
               block(g)%u_avg(i,j,k) = block(g)%u_sum(i,j,k)/ita1
               block(g)%v_avg(i,j,k) = block(g)%v_sum(i,j,k)/ita1
               block(g)%w_avg(i,j,k) = block(g)%w_sum(i,j,k)/ita1
               block(g)%p_avg(i,j,k) = block(g)%p_sum(i,j,k)/ita1
 30         CONTINUE
      END SUBROUTINE computeAvgData
!***********************************************************************

!***********************************************************************
      SUBROUTINE writeAvgoutput
       USE global
!!tec       USE mod_tecio
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       CHARACTER(len=70)  :: filename1
       CHARACTER(len=70)  :: filename2
       INTEGER  :: k, i, j, g
       REAL (KIND = 8) :: u_avg1, v_avg1, w_avg1
!!tec         CHARACTER*1 NULLCHR
!!tec       INTEGER*4 iStat, FileFormat, FileType, Debug, VIsDouble
!!tec       integer*4 ZoneType, IMax, JMax, KMax, ICellMax, JCellMax,KCellMax
!!tec       integer*4 StrandID, unused, IsBlock, NFConns, FNMode, ShrConn
!!tec       integer,  dimension(:), pointer :: null
!!tec       integer*4 nTotalPts
!!tec       integer*4,dimension(8) :: valueLoc
!!tec       real*8    SolTime
!!tec       REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:,:) :: un1, vn1, wn1

!!tec      Do g=1,nblocks
!!tec
!!tec       nx_var=block(g)%nx
!!tec       ny_var=block(g)%ny
!!tec       nz_var=block(g)%nz
!!tec      ALLOCATE (un1(nx_var+1,ny_var+1,nz_var+1), vn1(nx_var+1,ny_var+1,nz_var+1), wn1(nx_var+1,ny_var+1,nz_var+1))
!!tec      WRITE(filename1,1) ita,g,re
!!tec1     FORMAT('avgout/avg_fielddata.',i9.9,'.',i3.3,'.',f6.1,".plt")
!!tec     NULLCHR = CHAR(0)
!!tec     FileFormat = 0  !0= PLT;  1=SZPLT
!!tec     FileType = 0    !1= Grid; 2=Field; 0=Grid+Field
!!tec     Debug = 1
!!tec     VIsDouble = 0
!!tec
!!tec     ZoneType = 0
!!tec     IMax = block(g)%nx
!!tec     JMax = block(g)%ny
!!tec     KMax = block(g)%nz
!!tec     ICellMax = 0
!!tec     JCellMax = 0
!!tec     KCellMax = 0
!!tec
!!tec     StrandID = 1
!!tec     unused = 0
!!tec     isBlock = 1
!!tec     NFConns = 0
!!tec     FNMode = 0
!!tec     ShrConn = 0
!!tec     valueLoc(:) = 1
!!tec     SolTime = totime
!!tec
!!tec    !write(filename1,101) ita
!!tec    !101  format('out/field.',i10.10,'.plt')
!!tec
!!tec     iStat = tecIni142('SweptAngleStudy'//NULLCHR,  &
!!tec                       'x y z p  cell u v w'//NULLCHR,       &
!!tec                       filename1//NULLCHR,          &
!!tec                       './out/'//NULLCHR,    &   ! Scratch directory
!!tec                       FileFormat,                 &
!!tec                       FileType,                   &
!!tec                       Debug,                      &
!!tec                       VIsDouble)
!!tec
!!tec     iStat = tecZne142('Instantaneous Plots'//NULLCHR, &
!!tec                       ZoneType,               &
!!tec                       IMax,                   &
!!tec                       JMax,                   &
!!tec                       KMax,                   &
!!tec                       ICellMax,               &
!!tec                       JCellMax,               &
!!tec                       KCellMax,               &
!!tec                       SolTime,                &
!!tec                       StrandID,               &
!!tec                       unused,                 &
!!tec                       IsBlock,                &
!!tec                       NFConns,                &
!!tec                       FNMode,                 &
!!tec                       0,                      &
!!tec                       0,                      &
!!tec                       0,                      &
!!tec                       Null,                   &
!!tec                       valueLoc,               &
!!tec                       Null,                   &
!!tec                       ShrConn)
!!tec
!!tec     un1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%u_avg(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%u_avg(1:nx_var-1,2:ny_var,2:nz_var)) *0.5
!!tec     vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%v_avg(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%v_avg(2:nx_var,1:ny_var-1,2:nz_var)) *0.5
!!tec     wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%w_avg(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%w_avg(2:nx_var,2:ny_var,1:nz_var-1)) *0.5
!!tec
!!tec     nTotalPts = (nx_var)*(ny_var)*(nz_var)
!!tec    !nTotalPts = (nx_var)*(ny_var)*(1)
!!tec      iStat = tecDat142(nTotalPts,      block(g)%xp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!tec      iStat = tecDat142(nTotalPts,      block(g)%yp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!tec      iStat = tecDat142(nTotalPts,      block(g)%zp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!tec     !iStat = tecDat142(nTotalPts,     xp(2:nx_var,4)    , 0)
!!tec     !iStat = tecDat142(nTotalPts,     yp(2:ny_var,4)    , 0)
!!tec     !iStat = tecDat142(nTotalPts,     zp(2:nz_var,4)    , 0)
!!tec     iStat = tecDat142(nTotalPts, real(  block(g)%p_avg(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!tec     !iStat = tecDat142(nTotalPts, real( dmiu(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),     0)
!!tec     iStat = tecDat142(nTotalPts, real( block(g)%cell(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!tec     iStat = tecDat142(nTotalPts, real(un1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!tec     iStat = tecDat142(nTotalPts, real(vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!tec     iStat = tecDat142(nTotalPts, real(wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!tec     iStat = tecEnd142()
!!tec
!!tec       DEALLOCATE(un1,vn1,wn1)
!!tec      end do

           DO g=1,nblocks
           WRITE(filename1,1) ita1, g
1          FORMAT('avgout/fielddata.',i9.9,'.',i3.3,".dat")
           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
           WRITE(786,*)'variables= "x","y","z","u_avg","v_avg","w_avg","p_avg","totime","cellid"'
           WRITE(786,*) 'zone, ', 'i = ', block(g)%nx,' j = ', block(g)%ny, ' k = ', block(g)%nz
           !k = block(g)%nz/2
           DO 30 k = 2, block(g)%nz+1
           DO 30 j = 2, block(g)%ny+1
           DO 30 i = 2, block(g)%nx+1
              u_avg1 = 0.5*(block(g)%u_avg(i,j,k)+block(g)%u_avg(i-1,j,k))
              v_avg1 = 0.5*(block(g)%v_avg(i,j,k)+block(g)%v_avg(i,j-1,k))
              w_avg1 = 0.5*(block(g)%w_avg(i,j,k)+block(g)%w_avg(i,j,k-1))
              WRITE(786,*) block(g)%xp(i), block(g)%yp(j), block(g)%zp(k), u_avg1, v_avg1, w_avg1, block(g)%p_avg(i,j,k), totime, block(g)%cell(i,j,k)
30         CONTINUE
           CLOSE(786)
           END DO
      END SUBROUTINE writeAvgoutput
!***********************************************************************
         SUBROUTINE body_plot
         USE global
         IMPLICIT NONE
         INTEGER(KIND=8) :: inode, ielem, g
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

!****************************************************************************


      SUBROUTINE computeSumSqData
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER  :: n, k, i, j, g
	    DO  g=1,nblocks

           !!$acc parallel loop gang vector collapse (3) default(present)
            DO 30 k = 2, block(g)%nz+2
            DO 30 j = 2, block(g)%ny+2
            DO 30 i = 2, block(g)%nx+2
             !block(g)% u_sum(i,j,k,ph) =block(g)% u_sum(i,j,k,ph) + 0.5*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))
             !block(g)% v_sum(i,j,k,ph) =block(g)% v_sum(i,j,k,ph) + 0.5*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))
             !block(g)% w_sum(i,j,k,ph) =block(g)% w_sum(i,j,k,ph) + 0.5*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))
             !block(g)% p_sum(i,j,k,ph) =block(g)% p_sum(i,j,k,ph) + block(g)%p(i,j,k)
             !block(g)%uv_sum(i,j,k,ph) = block(g)%uv_sum(i,j,k,ph) + 0.5*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))*0.5*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))
             !block(g)%vw_sum(i,j,k,ph) = block(g)%vw_sum(i,j,k,ph) + 0.5*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))*0.5*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))
             !block(g)%uw_sum(i,j,k,ph) =block(g)% uw_sum(i,j,k,ph) + 0.5*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))*0.5*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))

             !block(g)%u2_sum(i,j,k,ph) = block(g)%u2_sum(i,j,k,ph) + (0.5*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k)))**2
             !block(g)%v2_sum(i,j,k,ph) = block(g)%v2_sum(i,j,k,ph) + (0.5*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k)))**2
             !block(g)%w2_sum(i,j,k,ph) =block(g)% w2_sum(i,j,k,ph) + (0.5*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1)))**2
             !block(g)%p2_sum(i,j,k,ph) = block(g)%p2_sum(i,j,k,ph) + block(g)%p(i,j,k)*block(g)%p(i,j,k)


               block(g)%uv_sum(i,j,k) = block(g)%uv_sum(i,j,k) + (block(g)%u(i,j,k) *block(g)%v(i,j,k))
               block(g)%vw_sum(i,j,k) = block(g)%vw_sum(i,j,k) + (block(g)%v(i,j,k) *block(g)%w(i,j,k))
               block(g)%uw_sum(i,j,k) =block(g)% uw_sum(i,j,k) + (block(g)%u(i,j,k) *block(g)%w(i,j,k))

              block(g)%u2_sum(i,j,k) = block(g)%u2_sum(i,j,k) + ((block(g)%u(i,j,k))**2)
              block(g)%v2_sum(i,j,k) = block(g)%v2_sum(i,j,k) + ((block(g)%v(i,j,k))**2)
              block(g)%w2_sum(i,j,k) = block(g)%w2_sum(i,j,k) + ((block(g)%w(i,j,k))**2)
              block(g)%p2_sum(i,j,k) = block(g)%p2_sum(i,j,k) + (block(g)%p(i,j,k)*block(g)%p(i,j,k))
 30         CONTINUE
            !!$acc end parallel

           END DO
      END SUBROUTINE computeSumSqData
!***********************************************************************
      SUBROUTINE computePhaseAvgData
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER  :: n, k, i, j,g

       !print*, nCycle, ph, phdif(ph)

           !!$acc update host(u_sum, v_sum, w_sum, p_sum, u2_sum, v2_sum, w2_sum, p2_sum, uv_sum, vw_sum,  uw_sum)

            DO g=1, nblocks

           !!$acc parallel loop gang vector collapse (3) default(present)
            DO 30 k = 1, block(g)%nz+2
            DO 30 j = 1, block(g)%ny+2
            DO 30 i = 1, block(g)%nx+2
              !block(g)%u_avg(i,j,k,ph) = block(g)%u_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%v_avg(i,j,k,ph) = block(g)%v_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%w_avg(i,j,k,ph) = block(g)%w_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%p_avg(i,j,k,ph) = block(g)%p_sum(i,j,k,ph)/(nCycle*phdif(ph))

              !block(g)%u2_avg(i,j,k,ph) = block(g)%u2_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%v2_avg(i,j,k,ph) = block(g)%v2_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%w2_avg(i,j,k,ph) = block(g)%w2_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%p2_avg(i,j,k,ph) = block(g)%p2_sum(i,j,k,ph)/(nCycle*phdif(ph))

              !block(g)%uv_avg(i,j,k,ph) = block(g)%uv_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%vw_avg(i,j,k,ph) = block(g)%vw_sum(i,j,k,ph)/(nCycle*phdif(ph))
              !block(g)%uw_avg(i,j,k,ph) = block(g)%uw_sum(i,j,k,ph)/(nCycle*phdif(ph))


              !block(g)%uflu_avg(i,j,k,ph) =   block(g)%u2_avg(i,j,k,ph) -  block(g)%u_avg(i,j,k,ph)**2
              !block(g)%vflu_avg(i,j,k,ph) =   block(g)%v2_avg(i,j,k,ph) -  block(g)%v_avg(i,j,k,ph)**2
              !block(g)%wflu_avg(i,j,k,ph) =   block(g)%w2_avg(i,j,k,ph) -  block(g)%w_avg(i,j,k,ph)**2
              !block(g)%pflu_avg(i,j,k,ph) =   block(g)%p2_avg(i,j,k,ph) -  block(g)%p_avg(i,j,k,ph)**2

              !block(g)%uvflu_avg(i,j,k,ph)=   block(g)%uv_avg(i,j,k,ph) -    block(g)%u_avg(i,j,k,ph)*block(g)%v_avg(i,j,k,ph)
              !block(g)%vwflu_avg(i,j,k,ph)=   block(g)%vw_avg(i,j,k,ph) -    block(g)%v_avg(i,j,k,ph)*block(g)%w_avg(i,j,k,ph)
              !block(g)%uwflu_avg(i,j,k,ph)=   block(g)%uw_avg(i,j,k,ph) -    block(g)%u_avg(i,j,k,ph)*block(g)%w_avg(i,j,k,ph)

              !block(g)%uflu_rms(i,j,k,ph) = sqrt(abs(block(g)%uflu_avg(i,j,k,ph)))
              !block(g)%vflu_rms(i,j,k,ph) = sqrt(abs(block(g)%vflu_avg(i,j,k,ph)))
              !block(g)%wflu_rms(i,j,k,ph) = sqrt(abs(block(g)%wflu_avg(i,j,k,ph)))
              !block(g)%pflu_rms(i,j,k,ph) = sqrt(abs(block(g)%pflu_avg(i,j,k,ph)))

               block(g)%u_avg(i,j,k) = block(g)%u_sum(i,j,k)/(ita1)
               block(g)%v_avg(i,j,k) = block(g)%v_sum(i,j,k)/(ita1)
               block(g)%w_avg(i,j,k) = block(g)%w_sum(i,j,k)/(ita1)
               block(g)%p_avg(i,j,k) = block(g)%p_sum(i,j,k)/(ita1)

               block(g)%u2_avg(i,j,k) = block(g)%u2_sum(i,j,k)/(ita1)
               block(g)%v2_avg(i,j,k) = block(g)%v2_sum(i,j,k)/(ita1)
               block(g)%w2_avg(i,j,k) = block(g)%w2_sum(i,j,k)/(ita1)
               block(g)%p2_avg(i,j,k) = block(g)%p2_sum(i,j,k)/(ita1)

               block(g)%uv_avg(i,j,k) = block(g)%uv_sum(i,j,k)/(ita1)
               block(g)%vw_avg(i,j,k) = block(g)%vw_sum(i,j,k)/(ita1)
               block(g)%uw_avg(i,j,k) = block(g)%uw_sum(i,j,k)/(ita1)


               block(g)%uflu_avg(i,j,k) =   block(g)%u2_avg(i,j,k) -  block(g)%u_avg(i,j,k)**2
               block(g)%vflu_avg(i,j,k) =   block(g)%v2_avg(i,j,k) -  block(g)%v_avg(i,j,k)**2
               block(g)%wflu_avg(i,j,k) =   block(g)%w2_avg(i,j,k) -  block(g)%w_avg(i,j,k)**2
               block(g)%pflu_avg(i,j,k) =   block(g)%p2_avg(i,j,k) -  block(g)%p_avg(i,j,k)**2
  block(g)%uvflu_avg(i,j,k)=block(g)%uv_avg(i,j,k)-block(g)%u_avg(i,j,k)*block(g)%v_avg(i,j,k)
  block(g)%vwflu_avg(i,j,k)=block(g)%vw_avg(i,j,k)-block(g)%v_avg(i,j,k)*block(g)%w_avg(i,j,k)
  block(g)%uwflu_avg(i,j,k)=block(g)%uw_avg(i,j,k)-block(g)%u_avg(i,j,k)*block(g)%w_avg(i,j,k)

               block(g)%uflu_rms(i,j,k) = sqrt(abs(block(g)%uflu_avg(i,j,k)))
               block(g)%vflu_rms(i,j,k) = sqrt(abs(block(g)%vflu_avg(i,j,k)))
               block(g)%wflu_rms(i,j,k) = sqrt(abs(block(g)%wflu_avg(i,j,k)))
               block(g)%pflu_rms(i,j,k) = sqrt(abs(block(g)%pflu_avg(i,j,k)))
 30         CONTINUE
            !!$acc end parallel

        END DO

      END SUBROUTINE computePhaseAvgData
!!***********************************************************************
    ! SUBROUTINE writeComputeSqSumData
    !   USE global
    !   implicit none
    !   INTEGER, PARAMETER :: rk = selected_real_kind(8)
    !   INTEGER::  i, j, k
    !  CHARACTER*70  filename1

    !  !IF(mod(ita1,5000).EQ.0.or.mod(ita1,86000).EQ.0)THEN! &
    !   IF(mod(ita1,4000).EQ.0)THEN! &
!   !        .or. ita1.eq.ph5e .or. ita1.eq.ph6e .or. ita1.eq.ph7e)THEN
    !   !$acc update host(u_sum, v_sum, w_sum, p_sum, u2_sum, v2_sum, w2_sum, p2_sum, uv_sum, vw_sum, uw_sum)
!   !       WRITE(filename1,1) nCycle, ita1
!1  !   FORMAT('out_gg/Results/result.',i2,2,'.',i9.9,".io")
!1  !   FORMAT('out_gg/Results/sumdata.',i2.2,'.',i9.9,".io")
!   !       OPEN (1,FILE=filename1,FORM='formatted')
    !  !OPEN (1,FILE=filename1,FORM='unformatted',access='stream')
    !   OPEN (1,FILE='out_gg/Results/sumdata4k.io',FORM='unformatted',access='stream')

!   !    DO 30 k = 1, nz+2
!   !    DO 30 j = 1, ny+2
!   !    DO 30 i = 1, nx+2
    !     !WRITE(1,*) u_sum(i,j,k,ph), v_sum(i,j,k,ph), w_sum(i,j,k,ph), p_sum(i,j,k,ph),  u2_sum(i,j,k,ph), v2_sum(i,j,k,ph), w2_sum(i,j,k,ph), p2_sum(i,j,k,ph), uv_sum(i,j,k,ph), vw_sum(i,j,k,ph),  uw_sum(i,j,k,ph), ph
    !     WRITE(1) u_sum, v_sum, w_sum, p_sum,  u2_sum, v2_sum, w2_sum, p2_sum, uv_sum, vw_sum, uw_sum
! 30!    CONTINUE
    !    CLOSE(1)
    !    END IF
    !   IF(mod(ita1,5000).EQ.0)THEN! &
!   !        .or. ita1.eq.ph5e .or. ita1.eq.ph6e .or. ita1.eq.ph7e)THEN
    !   !$acc update host(u_sum, v_sum, w_sum, p_sum, u2_sum, v2_sum, w2_sum, p2_sum, uv_sum, vw_sum, uw_sum)
!   !       WRITE(filename1,1) nCycle, ita1
!1  !   FORMAT('out_gg/Results/result.',i2,2,'.',i9.9,".io")
!1  !   FORMAT('out_gg/Results/sumdata.',i2.2,'.',i9.9,".io")
!   !       OPEN (1,FILE=filename1,FORM='formatted')
    !  !OPEN (1,FILE=filename1,FORM='unformatted',access='stream')
    !   OPEN (1,FILE='out_gg/Results/sumdata5k.io',FORM='unformatted',access='stream')

!   !    DO 30 k = 1, nz+2
!   !    DO 30 j = 1, ny+2
!   !    DO 30 i = 1, nx+2
    !     !WRITE(1,*) u_sum(i,j,k,ph), v_sum(i,j,k,ph), w_sum(i,j,k,ph), p_sum(i,j,k,ph),  u2_sum(i,j,k,ph), v2_sum(i,j,k,ph), w2_sum(i,j,k,ph), p2_sum(i,j,k,ph), uv_sum(i,j,k,ph), vw_sum(i,j,k,ph),  uw_sum(i,j,k,ph), ph
    !     WRITE(1) u_sum, v_sum, w_sum, p_sum,  u2_sum, v2_sum, w2_sum, p2_sum, uv_sum, vw_sum, uw_sum
! 30!    CONTINUE
    !    CLOSE(1)
    !    END IF
    !   IF(mod(ita1,86000).EQ.0)THEN! &
!   !        .or. ita1.eq.ph5e .or. ita1.eq.ph6e .or. ita1.eq.ph7e)THEN
    !   !$acc update host(u_sum, v_sum, w_sum, p_sum, u2_sum, v2_sum, w2_sum, p2_sum, uv_sum, vw_sum, uw_sum)
!   !       WRITE(filename1,1) nCycle, ita1
!1  !   FORMAT('out_gg/Results/result.',i2,2,'.',i9.9,".io")
!1  !   FORMAT('out_gg/Results/sumdata.',i2.2,'.',i9.9,".io")
!   !       OPEN (1,FILE=filename1,FORM='formatted')
    !  !OPEN (1,FILE=filename1,FORM='unformatted',access='stream')
    !   OPEN (1,FILE='out_gg/Results/sumdata86k.io',FORM='unformatted',access='stream')

!   !    DO 30 k = 1, nz+2
!   !    DO 30 j = 1, ny+2
!   !    DO 30 i = 1, nx+2
    !     !WRITE(1,*) u_sum(i,j,k,ph), v_sum(i,j,k,ph), w_sum(i,j,k,ph), p_sum(i,j,k,ph),  u2_sum(i,j,k,ph), v2_sum(i,j,k,ph), w2_sum(i,j,k,ph), p2_sum(i,j,k,ph), uv_sum(i,j,k,ph), vw_sum(i,j,k,ph),  uw_sum(i,j,k,ph), ph
    !     WRITE(1) u_sum, v_sum, w_sum, p_sum,  u2_sum, v2_sum, w2_sum, p2_sum, uv_sum, vw_sum, uw_sum
! 30!    CONTINUE
    !    CLOSE(1)
    !    END IF
    ! END SUBROUTINE  writeComputeSqSumData
!***********************************************************************


      SUBROUTINE writeAvgoutput_fl
       USE global
!!tec       USE mod_tecio
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       CHARACTER(len=70)  :: filename1
       CHARACTER(len=70)  :: filename2
       INTEGER  :: k, i, j, g
       REAL (KIND = 8) :: u_avg1, v_avg1, w_avg1
!!tec         CHARACTER*1 NULLCHR
!!tec       INTEGER*4 iStat, FileFormat, FileType, Debug, VIsDouble
!!tec       integer*4 ZoneType, IMax, JMax, KMax, ICellMax, JCellMax,KCellMax
!!tec       integer*4 StrandID, unused, IsBlock, NFConns, FNMode, ShrConn
!!tec       integer,  dimension(:), pointer :: null
!!tec       integer*4 nTotalPts
!!tec       integer*4,dimension(8) :: valueLoc
!!tec       real*8    SolTime
!!tec       REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:,:) :: un1, vn1, wn1
!!tec
!!tec      Do g=1,nblocks
!!tec
!!tec       nx_var=block(g)%nx
!!tec       ny_var=block(g)%ny
!!tec       nz_var=block(g)%nz
!!tec      ALLOCATE (un1(nx_var+1,ny_var+1,nz_var+1), vn1(nx_var+1,ny_var+1,nz_var+1),wn1(nx_var+1,ny_var+1,nz_var+1))
!!tec      WRITE(filename1,1) ita,g,re
!!tec1     FORMAT('avgout/flavg_fielddata.',i9.9,'.',i3.3,'.',f6.1,".plt")
!!tec     NULLCHR = CHAR(0)
!!tec     FileFormat = 0  !0= PLT;  1=SZPLT
!!tec     FileType = 0    !1= Grid; 2=Field; 0=Grid+Field
!!tec     Debug = 1
!!tec     VIsDouble = 0
!!tec
!!tec     ZoneType = 0
!!tec     IMax = block(g)%nx
!!tec     JMax = block(g)%ny
!!tec     KMax = block(g)%nz
!!tec     ICellMax = 0
!!tec     JCellMax = 0
!!tec     KCellMax = 0
!!tec
!!tec     StrandID = 1
!!tec     unused = 0
!!tec     isBlock = 1
!!tec     NFConns = 0
!!tec     FNMode = 0
!!tec     ShrConn = 0
!!tec     valueLoc(:) = 1
!!tec     SolTime = totime
!!tec
!!tec    !write(filename1,101) ita
!!tec    !101  format('out/field.',i10.10,'.plt')
!!tec
!!tec     iStat = tecIni142('SweptAngleStudy'//NULLCHR,  &
!!tec                       'x y z p  cell u v w'//NULLCHR,       &
!!tec                       filename1//NULLCHR,          &
!!tec                       './out/'//NULLCHR,    &   ! Scratch directory
!!tec                       FileFormat,                 &
!!tec                       FileType,                   &
!!tec                       Debug,                      &
!!tec                       VIsDouble)
!!tec
!!tec     iStat = tecZne142('Instantaneous Plots'//NULLCHR, &
!!tec                       ZoneType,               &
!!tec                       IMax,                   &
!!tec                       JMax,                   &
!!tec                       KMax,                   &
!!tec                       ICellMax,               &
!!tec                       JCellMax,               &
!!tec                       KCellMax,               &
!!tec                       SolTime,                &
!!tec                       StrandID,               &
!!tec                       unused,                 &
!!tec                       IsBlock,                &
!!tec                       NFConns,                &
!!tec                       FNMode,                 &
!!tec                       0,                      &
!!tec                       0,                      &
!!tec                       0,                      &
!!tec                       Null,                   &
!!tec                       valueLoc,               &
!!tec                       Null,                   &
!!tec                       ShrConn)
!!tec
!!tec     un1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%uflu_rms(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%uflu_rms(1:nx_var-1,2:ny_var,2:nz_var)) *0.5
!!tec     vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%vflu_rms(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%vflu_rms(2:nx_var,1:ny_var-1,2:nz_var)) *0.5
!!tec     wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%wflu_rms(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%wflu_rms(2:nx_var,2:ny_var,1:nz_var-1)) *0.5
!!tec
!!tec     nTotalPts = (nx_var)*(ny_var)*(nz_var)
!!tec    !nTotalPts = (nx_var)*(ny_var)*(1)
!!tec      iStat = tecDat142(nTotalPts,      block(g)%xp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!tec      iStat = tecDat142(nTotalPts,      block(g)%yp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!tec      iStat = tecDat142(nTotalPts,      block(g)%zp1(2:nx_var+1,2:ny_var+1,2:nz_var+1) ,   0)
!!tec     !iStat = tecDat142(nTotalPts,     xp(2:nx_var,4)    , 0)
!!tec     !iStat = tecDat142(nTotalPts,     yp(2:ny_var,4)    , 0)
!!tec     !iStat = tecDat142(nTotalPts,     zp(2:nz_var,4)    , 0)
!!tec     iStat = tecDat142(nTotalPts, real(  block(g)%pflu_rms(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!tec     !iStat = tecDat142(nTotalPts, real( dmiu(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),     0)
!!tec     iStat = tecDat142(nTotalPts, real( block(g)%cell(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),0)
!!tec     iStat = tecDat142(nTotalPts, real(un1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!tec     iStat = tecDat142(nTotalPts, real(vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!tec     iStat = tecDat142(nTotalPts, real(wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1),4), 0)
!!tec     iStat = tecEnd142()
!!tec
!!tec       DEALLOCATE(un1,vn1,wn1)
!!tec      end do

           DO g=1,nblocks
           WRITE(filename1,1) ita1, g
1          FORMAT('avgout/fielddata.',i9.9,'.',i3.3,".dat")
           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
           WRITE(786,*)'variables= "x","y","z","u_avg","v_avg","w_avg","p_avg","totime","cellid"'
           WRITE(786,*) 'zone, ', 'i = ', block(g)%nx,' j = ', block(g)%ny, ' k = ', block(g)%nz
           !k = block(g)%nz/2
           DO 30 k = 2, block(g)%nz+1
           DO 30 j = 2, block(g)%ny+1
           DO 30 i = 2, block(g)%nx+1
              u_avg1 = 0.5*(block(g)%u_avg(i,j,k)+block(g)%u_avg(i-1,j,k))
              v_avg1 = 0.5*(block(g)%v_avg(i,j,k)+block(g)%v_avg(i,j-1,k))
              w_avg1 = 0.5*(block(g)%w_avg(i,j,k)+block(g)%w_avg(i,j,k-1))
              WRITE(786,*) block(g)%xp(i), block(g)%yp(j), block(g)%zp(k), u_avg1, v_avg1, w_avg1, block(g)%p_avg(i,j,k), totime, block(g)%cell(i,j,k)
30         CONTINUE
           CLOSE(786)
           END DO
      END SUBROUTINE writeAvgoutput_fl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE write_fl_points

       USE global

       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       CHARACTER(len=150)  :: filename1
       CHARACTER(len=150)  :: filename2
       INTEGER  :: n,k, i, j, g
       REAL (KIND = 8) :: uu, vv, ww, uv, vw, wu

        g=1


           if ( mod(block(g)%nx+1,2) == 0) then
                i=(block(g)%nx+1)/2
           else

                i=((block(g)%nx+1)/2) +1
           endif

           if ( mod(block(g)%nz+1,2) == 0) then
                k=(block(g)%nz+1)/2
           else

                k=((block(g)%nz+1)/2) +1
           endif

 126      FORMAT(' ',7F15.7)
            do n=1,block(g)%ny+2

                if( block(g)%yv(n) >18) then
                        j=n-1
                        exit
                end if
            enddo

                uu= block(g)%u(i,j,k)**2
                vv= block(g)%v(i,j,k)**2
                ww= block(g)%w(i,j,k)**2
                uv= block(g)%u(i,j,k) * block(g)%v(i,j,k)
                vw= block(g)%w(i,j,k) * block(g)%v(i,j,k)
                wu= block(g)%u(i,j,k) * block(g)%w(i,j,k)

            WRITE(filename1,1) j
 1          FORMAT('out/fl_val_',i6.6,'.dat')
         OPEN(111,FILE=filename1,ACCESS='Append',STATUS='unknown')
         WRITE(111,126)  uu,vv,ww,uv,vw,wu, totime
         CLOSE(111)


            do n=1,block(g)%ny+2

                if( block(g)%yv(n) >24) then
                        j=n-1
                        exit
                end if
            enddo

                uu= block(g)%u(i,j,k)**2
                vv= block(g)%v(i,j,k)**2
                ww= block(g)%w(i,j,k)**2
                uv= block(g)%u(i,j,k) * block(g)%v(i,j,k)
                vw= block(g)%w(i,j,k) * block(g)%v(i,j,k)
                wu= block(g)%u(i,j,k) * block(g)%w(i,j,k)

            WRITE(filename1,11) j
 11          FORMAT('out/fl_val_',i6.6,'.dat')
         OPEN(111,FILE=filename1,ACCESS='Append',STATUS='unknown')
         WRITE(111,126)  uu,vv,ww,uv,vw,wu, totime
         CLOSE(111)



            do n=1,block(g)%ny+2

                if( block(g)%yv(n) >31) then
                        j=n-1
                        exit
                end if
            enddo

                uu= block(g)%u(i,j,k)**2
                vv= block(g)%v(i,j,k)**2
                ww= block(g)%w(i,j,k)**2
                uv= block(g)%u(i,j,k) * block(g)%v(i,j,k)
                vw= block(g)%w(i,j,k) * block(g)%v(i,j,k)
                wu= block(g)%u(i,j,k) * block(g)%w(i,j,k)

            WRITE(filename1,111) j
 111          FORMAT('out/fl_val_',i6.6,'.dat')
         OPEN(111,FILE=filename1,ACCESS='Append',STATUS='unknown')
         WRITE(111,126)  uu,vv,ww,uv,vw,wu, totime
         CLOSE(111)

            do n=1,block(g)%ny+2

                if( block(g)%yv(n) >40) then
                        j=n-1
                        exit
                end if
            enddo

                uu= block(g)%u(i,j,k)**2
                vv= block(g)%v(i,j,k)**2
                ww= block(g)%w(i,j,k)**2
                uv= block(g)%u(i,j,k) * block(g)%v(i,j,k)
                vw= block(g)%w(i,j,k) * block(g)%v(i,j,k)
                wu= block(g)%u(i,j,k) * block(g)%w(i,j,k)

            WRITE(filename1,1111) j
 1111    FORMAT('out/fl_val_',i6.6,'.dat')
         OPEN(111,FILE=filename1,ACCESS='Append',STATUS='unknown')
         WRITE(111,126)  uu,vv,ww,uv,vw,wu, totime
         CLOSE(111)



          g=2
           if ( mod(block(g)%nx+1,2) == 0) then
                i=(block(g)%nx+1)/2
           else

                i=((block(g)%nx+1)/2) +1
           endif

           if ( mod(block(g)%nz+1,2) == 0) then
                k=(block(g)%nz+1)/2
           else

                k=((block(g)%nz+1)/2) +1
           endif

            do n=1,block(g)%ny+2

                if( block(g)%yv(n) >46) then
                        j=n-1
                        exit
                end if
            enddo

                uu= block(g)%u(i,j,k)**2
                vv= block(g)%v(i,j,k)**2
                ww= block(g)%w(i,j,k)**2
                uv= block(g)%u(i,j,k) * block(g)%v(i,j,k)
                vw= block(g)%w(i,j,k) * block(g)%v(i,j,k)
                wu= block(g)%u(i,j,k) * block(g)%w(i,j,k)

            WRITE(filename1,11111) j
11111    FORMAT('out/fl_val_',i6.6,'.dat')
         OPEN(111,FILE=filename1,ACCESS='Append',STATUS='unknown')
         WRITE(111,126)  uu,vv,ww,uv,vw,wu, totime
         CLOSE(111)


        end subroutine write_fl_points
!!!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!!      SUBROUTINE writeOutput
!!       USE global
!!       USE mod_tecio
!!       INTEGER, PARAMETER :: rk = selected_real_kind(8)
!!       CHARACTER*120  filename1
!!       CHARACTER*70  filename2
!!       INTEGER  :: k, i, j  ,g
!!       REAL (KIND = 8) :: u1, v1, w1
!!
!!         CHARACTER*1 NULLCHR
!!       INTEGER*4 iStat, FileFormat, FileType, Debug, VIsDouble
!!       integer*4 ZoneType, IMax, JMax, KMax, ICellMax, JCellMax,KCellMax
!!       integer*4 StrandID, unused, IsBlock, NFConns, FNMode, ShrConn
!!       integer,  dimension(:), pointer :: null
!!       integer*4 nTotalPts
!!       integer*4,dimension(8) :: valueLoc
!!       real*8    SolTime
!!       REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:,:) :: un1, vn1, wn1, pn1
!!
!!        IF(mod(ita,2).EQ.0 )THEN
!!        DO g=1,nblocks
!!!           WRITE(filename1,1) ita,g
!!!1          FORMAT('out/w_fielddata.',i9.9,'.',i3.3,".dat")
!!!           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
!!!           WRITE(786,*)'variables= "x","y","z","u","v","w","p","totime","cellid"'
!!!           WRITE(786,*) 'zone, ', 'i = ', block(g)%nx+1,' j = ', block(g)%ny+1, ' k = ', block(g)%nz+1
!!            !k = block(g)%nz/2
!!       nx_var=block(g)%nx
!!       ny_var=block(g)%ny
!!       nz_var=block(g)%nz
!!      ALLOCATE (un1(nx_var+1,ny_var+1,nz_var+1), vn1(nx_var+1,ny_var+1,nz_var+1),wn1(nx_var+1,ny_var+1,nz_var+1),pn1(nx_var+1,ny_var+1,nz_var+1))
!!            DO 30 k = 1, block(g)%nz+1
!!            DO 30 j = 1, block(g)%ny+1
!!            DO 30 i = 1, block(g)%nx+1
!!
!!             bl_intp_y1= block(g)%yp(j)
!!             bl_intp_y2= block(g)%yp(j+1)
!!             bl_intp_x1= block(g)%xp(i)
!!             bl_intp_x2= block(g)%xp(i+1)
!!             bl_intp_f1= block(g)%p(i,j,k)
!!             bl_intp_f2= block(g)%p(i+1,j,k)
!!             bl_intp_f3= block(g)%p(i+1,j+1,k)
!!             bl_intp_f4= block(g)%p(i,j+1,k)
!!             bl_intp_valy=block(g)%y1(j+1)
!!             bl_intp_valx=block(g)%x1(i+1)
!!           call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)
!!
!!             p_new1=bl_interp_ans
!!
!!
!!             bl_intp_y1= block(g)%yp(j)
!!             bl_intp_y2= block(g)%yp(j+1)
!!             bl_intp_x1= block(g)%xp(i)
!!             bl_intp_x2= block(g)%xp(i+1)
!!             bl_intp_f1= block(g)%p(i,j,k+1)
!!             bl_intp_f2= block(g)%p(i+1,j,k+1)
!!             bl_intp_f3= block(g)%p(i+1,j+1,k+1)
!!             bl_intp_f4= block(g)%p(i,j+1,k+1)
!!             bl_intp_valy=block(g)%y1(j+1)
!!             bl_intp_valx=block(g)%x1(i+1)
!!           call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)
!!
!!             p_new2=bl_interp_ans
!!
!!
!!
!!             l_intp_valx=block(g)%z1(k+1)
!!             l_intp_x1= 0.5*(block(g)%z1(k)+block(g)%z1(k+1))
!!             l_intp_x2= 0.5*(block(g)%z1(k+1)+block(g)%z1(k+2))
!!             l_intp_y1=p_new1
!!             l_intp_y2=p_new2
!!             call linearInterp(l_intp_valx,l_intp_x1,l_intp_x2,l_intp_y1,l_intp_y2,l_intp_valy)
!!             p_final=l_intp_valy
!!
!!
!!             bl_intp_y1= 0.5*(block(g)%z1(k)+block(g)%z1(k+1))
!!             bl_intp_y2= 0.5*(block(g)%z1(k+1)+block(g)%z1(k+2))
!!             bl_intp_x1= 0.5*(block(g)%x1(i)+block(g)%x1(i+1))
!!             bl_intp_x2= 0.5*(block(g)%x1(i+1)+block(g)%x1(i+2))
!!             bl_intp_f1= block(g)%v(i,j,k)
!!             bl_intp_f2= block(g)%v(i+1,j,k)
!!             bl_intp_f3= block(g)%v(i+1,j,k+1)
!!             bl_intp_f4= block(g)%v(i,j,k+1)
!!             bl_intp_valy=block(g)%z1(k+1)
!!             bl_intp_valx=block(g)%x1(i+1)
!!           call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)
!!
!!             v_new=bl_interp_ans
!!
!!             bl_intp_y1= 0.5*(block(g)%z1(k)+block(g)%z1(k+1))
!!             bl_intp_y2= 0.5*(block(g)%z1(k+1)+block(g)%z1(k+2))
!!             bl_intp_x1= 0.5*(block(g)%y1(j)+block(g)%y1(j+1))
!!             bl_intp_x2= 0.5*(block(g)%y1(j+1)+block(g)%y1(j+2))
!!             bl_intp_f1= block(g)%u(i,j,k)
!!             bl_intp_f2= block(g)%u(i,j,k+1)
!!             bl_intp_f3= block(g)%u(i,j+1,k+1)
!!             bl_intp_f4= block(g)%u(i,j+1,k)
!!             bl_intp_valy=block(g)%z1(k+1)
!!             bl_intp_valx=block(g)%y1(j+1)
!!           call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)
!!
!!             u_new=bl_interp_ans
!!
!!             bl_intp_y1= 0.5*(block(g)%y1(j)+block(g)%y1(j+1))
!!             bl_intp_y2= 0.5*(block(g)%y1(j+1)+block(g)%y1(j+2))
!!             bl_intp_x1= 0.5*(block(g)%x1(i)+block(g)%x1(i+1))
!!             bl_intp_x2= 0.5*(block(g)%x1(i+1)+block(g)%x1(i+2))
!!             bl_intp_f1= block(g)%w(i,j,k)
!!             bl_intp_f2= block(g)%w(i+1,j,k)
!!             bl_intp_f3= block(g)%w(i+1,j+1,k)
!!             bl_intp_f4= block(g)%w(i,j+1,k)
!!             bl_intp_valy=block(g)%y1(j+1)
!!             bl_intp_valx=block(g)%x1(i+1)
!!           call billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans)
!!
!!             w_new=bl_interp_ans
!!
!!             un1(i,j,k)=u_new
!!             vn1(i,j,k)=v_new
!!             wn1(i,j,k)=w_new
!!             pn1(i,j,k)=p_final
!!
!!
!!              ! WRITE(786,*) block(g)%x1(i+1), block(g)%y1(j+1), block(g)%z1(k+1), u_new, v_new, w_new, p_final, totime, block(g)%cell(i,j,k)
!!
!! 30         CONTINUE
!!             CLOSE(786)
!!      WRITE(filename1,1)ita,g,re
!!1     FORMAT('out/w_fielddata.',i9.9,'.',i3.3,'.',f6.1,".plt")
!!     NULLCHR = CHAR(0)
!!     FileFormat = 0  !0= PLT;  1=SZPLT
!!     FileType = 0    !1= Grid; 2=Field; 0=Grid+Field
!!     Debug = 1
!!     VIsDouble = 0
!!
!!     ZoneType = 0
!!     IMax = block(g)%nx+1
!!     JMax = block(g)%ny+1
!!     KMax = block(g)%nz+1
!!     ICellMax = 0
!!     JCellMax = 0
!!     KCellMax = 0
!!
!!     StrandID = 1
!!     unused = 0
!!     isBlock = 1
!!     NFConns = 0
!!     FNMode = 0
!!     ShrConn = 0
!!     valueLoc(:) = 1
!!     SolTime = totime
!!
!!    !write(filename1,101) ita
!!    !101  format('out/field.',i10.10,'.plt')
!!
!!     iStat = tecIni142('SweptAngleStudy'//NULLCHR,  &
!!                       'x y z p  cell u v w'//NULLCHR,       &
!!                       filename1//NULLCHR,          &
!!                       './out/'//NULLCHR,    &   ! Scratch directory
!!                       FileFormat,                 &
!!                       FileType,                   &
!!                       Debug,                      &
!!                       VIsDouble)
!!
!!     iStat = tecZne142('Instantaneous Plots'//NULLCHR, &
!!                       ZoneType,               &
!!                       IMax,                   &
!!                       JMax,                   &
!!                       KMax,                   &
!!                       ICellMax,               &
!!                       JCellMax,               &
!!                       KCellMax,               &
!!                       SolTime,                &
!!                       StrandID,               &
!!                       unused,                 &
!!                       IsBlock,                &
!!                       NFConns,                &
!!                       FNMode,                 &
!!                       0,                      &
!!                       0,                      &
!!                       0,                      &
!!                       Null,                   &
!!                       valueLoc,               &
!!                       Null,                   &
!!                       ShrConn)
!!
!!   ! un1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%u(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%u(1:nx_var-1,2:ny_var,2:nz_var)) *0.5
!!   ! vn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%v(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%v(2:nx_var,1:ny_var-1,2:nz_var)) *0.5
!!   ! wn1(2:nx_var+1,2:ny_var+1,2:nz_var+1)= (block(g)%w(2:nx_var+1,2:ny_var+1,2:nz_var+1) +block(g)%w(2:nx_var,2:ny_var,1:nz_var-1)) *0.5
!!
!!     nTotalPts = (nx_var+1)*(ny_var+1)*(nz_var+1)
!!    !nTotalPts = (nx_var)*(ny_var)*(1)
!!      iStat = tecDat142(nTotalPts,      block(g)%xpn1(1:nx_var+1,1:ny_var+1,1:nz_var+1) ,   0)
!!      iStat = tecDat142(nTotalPts,      block(g)%ypn1(1:nx_var+1,1:ny_var+1,1:nz_var+1) ,   0)
!!      iStat = tecDat142(nTotalPts,      block(g)%zpn1(1:nx_var+1,1:ny_var+1,1:nz_var+1) ,   0)
!!     !iStat = tecDat142(nTotalPts,     xp(2:nx_var,4)    , 0)
!!     !iStat = tecDat142(nTotalPts,     yp(2:ny_var,4)    , 0)
!!     !iStat = tecDat142(nTotalPts,     zp(2:nz_var,4)    , 0)
!!     iStat = tecDat142(nTotalPts, real(  pn1(1:nx_var+1,1:ny_var+1,1:nz_var+1),4),0)
!!     !iStat = tecDat142(nTotalPts, real( dmiu(2:nx_var+1,2:ny_var+1,2:nz_var+1),4),     0)
!!     iStat = tecDat142(nTotalPts, real( block(g)%cell(1:nx_var+1,1:ny_var+1,1:nz_var+1),4),0)
!!     iStat = tecDat142(nTotalPts, real(un1(1:nx_var+1,1:ny_var+1,1:nz_var+1),4), 0)
!!     iStat = tecDat142(nTotalPts, real(vn1(1:nx_var+1,1:ny_var+1,1:nz_var+1),4), 0)
!!     iStat = tecDat142(nTotalPts, real(wn1(1:nx_var+1,1:ny_var+1,1:nz_var+1),4), 0)
!!     iStat = tecEnd142()
!!
!!       DEALLOCATE(un1,vn1,wn1,pn1)
!!        END DO
!!         ENDIF
!!      END SUBROUTINE writeOutput
!!!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

        SUBROUTINE fineUpdate
        USE global
        IMPLICIT NONE
        REAL (KIND=8) :: bl_intp_valx,bl_intp_valy,bl_intp_valz, bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (KIND=8) :: bl_intp_z1,bl_intp_z2
        REAL (KIND=8) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (KIND=8) :: bl_interp_ans, bl_interp_ans1, bl_interp_ans2
        REAL (KIND=8) :: aa1,aa2,aa3,aa4,aamin

        INTEGER(KIND=8) :: i,j,k,q,s, varx1,varx2, vary1, vary2, l,tar_x, tar_y, loc_x, loc_y, coarse_indx,g, a_blk_no, b_blk_no
        INTEGER(KIND=8) :: varz1,varz2, tar_z, loc_z






        DO g=1,intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk

        DO k=2,intfr(g)%counterzp-1
        DO j=2,intfr(g)%counteryp-1
        DO i=2,intfr(g)%counterxp-1
        varx1=intfr(g)%px_interface_det(2,i)
        varx2=intfr(g)%px_interface_det(3,i)
        vary1=intfr(g)%py_interface_det(2,j)
        vary2=intfr(g)%py_interface_det(3,j)
        varz1=intfr(g)%pz_interface_det(2,k)
        varz2=intfr(g)%pz_interface_det(3,k)

        loc_x=intfr(g)%px_interface_det(1,i)
        loc_y=intfr(g)%py_interface_det(1,j)
        loc_z=intfr(g)%pz_interface_det(1,k)

      ! if ( i .eq. 1 .or. i .eq. intfr(g)%counterxp .or. k .eq. 1 .or. k .eq. intfr(g)%counteryp)then

      ! DO l=vary1,vary2
      ! DO j=varx1,varx2
      !         tar_x=j
      !         !tar_y=1
      !         tar_y=l
      !
      !         block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x,loc_y)
      ! ENDDO
      ! ENDDO
      !else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
               !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) ==0)then
       !if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
       !      block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and.  &
       !      block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
       !      block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then

 !!!!p1
                bl_intp_valx=block(b_blk_no)%xp(tar_x)
                bl_intp_valy=block(b_blk_no)%yp(tar_y)
                bl_intp_valz=block(b_blk_no)%zp(tar_z)
                bl_intp_x1=block(a_blk_no)%xp(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yp(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zp(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zp(loc_z+1)
                bl_intp_f1=block(a_blk_no)%p(loc_x-1,loc_y-1,loc_z-1)
                bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1,loc_z-1)
                bl_intp_f4=block(a_blk_no)%p(loc_x-1,loc_y+1,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%p(loc_x-1,loc_y-1,loc_z+1)
                bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%p(loc_x-1,loc_y+1,loc_z+1)
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans

              ! else

              ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
              ! aa1=(block(a_blk_no)%xp(loc_x-1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y-1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=aa1
              ! elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
              ! aa2=(block(a_blk_no)%xp(loc_x-1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y+1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=dmin1(aamin,aa2)
              ! elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
              ! aa3=(block(a_blk_no)%xp(loc_x+1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y+1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=dmin1(aamin,aa3)
              ! elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
              ! aa4=(block(a_blk_no)%xp(loc_x+1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y-1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=dmin1(aamin,aa4)
              ! endif

              ! if (aamin .eq. aa1)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x-1,loc_y-1)
              ! elseif(aamin .eq. aa2)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x-1,loc_y+1)
              ! elseif(aamin .eq. aa3)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! elseif (aamin .eq. aa4)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x+1,loc_y-1)
              ! endif

              ! endif
                endif
                ENDDO
                ENDDO
                ENDDO
             !   endif
                ENDDO
                ENDDO
                ENDDO


       !DO l=vary1,vary2
       !DO j=varx1,varx2
       !        tar_x=j
       !        !tar_y=1
       !        tar_y=l
       !       !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
       !       if (block(b_blk_no)%cell(tar_x,tar_y) .eq.0)then
!!!!!p1
       !     !! bl_intp_valx=block(b_blk_no)%xp(tar_x)
       !     !! bl_intp_valy=block(b_blk_no)%yp(tar_y)
       !     !! bl_intp_x1=block(a_blk_no)%xp(loc_x)
       !     !! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
       !     !! bl_intp_y1=block(a_blk_no)%yp(loc_y)
       !     !! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
       !     !! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
       !     !! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
       !     !! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
       !     !! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
       !     !!
       !     !! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
       !     !! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
       !     !! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
       !     !! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
       !     !! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
       !     !! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
       !     !! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
       !     !! bl_intp_num= bl_intp_first_term + bl_intp_second_term
       !
       !     !! bl_interp_ans=(bl_intp_num/bl_intp_deno)
       !     !   block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_
       !      endif
              ! tar_x=j
              ! !tar_y=2
              ! tar_y=k
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,1)
              !
!!!!!p2
              ! bl_intp_valx=block(b_blk_no)%xp(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yp(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xp(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yp(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans
              !


              ! tar_x=j
              ! !tar_y=block(b_blk_no)%jtc_en
              ! tar_y=k
              ! coarse_indx=intfr(g)%counteryp
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,coarse_indx)
              !
!!!!!pny
              ! bl_intp_valx=block(b_blk_no)%xp(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yp(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xp(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yp(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans


              ! tar_x=j
              ! tar_y=block(b_blk_no)%jtc_en-1
              ! coarse_indx=intfr(g)%counteryp
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,coarse_indx)
              !
!!!!!pny-1
              ! bl_intp_valx=block(b_blk_no)%xp(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yp(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xp(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yp(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans

!!!!!u
        DO k=2,intfr(g)%counterzu-1
        DO j=2,intfr(g)%counteryu-1
        DO i=2,intfr(g)%counterxu-1
        varx1=intfr(g)%ux_interface_det(2,i)
        varx2=intfr(g)%ux_interface_det(3,i)
        vary1=intfr(g)%uy_interface_det(2,j)
        vary2=intfr(g)%uy_interface_det(3,j)
        varz1=intfr(g)%uz_interface_det(2,k)
        varz2=intfr(g)%uz_interface_det(3,k)


        loc_x=intfr(g)%ux_interface_det(1,i)
        loc_y=intfr(g)%uy_interface_det(1,j)
        loc_z=intfr(g)%uz_interface_det(1,k)
      ! if (k .eq. 1 .or. i .eq.1 .or. k .eq. intfr(g)%counteryu .or. i .eq. intfr(g)%counterxu)then
      !

      ! DO l=vary1,vary2
      ! DO j=varx1,varx2
      !         tar_x=j
      !         !tar_y=1
      !         tar_y=l
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x-1,loc_y)
      !ENDDO
      !ENDDO
      !
      ! else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
              ! if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) ==0)then
      ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then
      !
!!!!!p
                bl_intp_valx=block(b_blk_no)%xu(tar_x)
                bl_intp_valy=block(b_blk_no)%yu(tar_y)
                bl_intp_valz=block(b_blk_no)%zu(tar_z)
                bl_intp_x1=block(a_blk_no)%xu(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yu(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zu(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zu(loc_z+1)
               !bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
                bl_intp_f1=block(a_blk_no)%u(loc_x-2,loc_y-1,loc_z-1)
                bl_intp_f2=block(a_blk_no)%u(loc_x,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%u(loc_x,loc_y+1,loc_z-1)
                bl_intp_f4=block(a_blk_no)%u(loc_x-2,loc_y+1,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%u(loc_x-2,loc_y-1,loc_z+1)
                bl_intp_f2=block(a_blk_no)%u(loc_x,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%u(loc_x,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%u(loc_x-2,loc_y+1,loc_z+1)
                !block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans
!                write(*,1221)'locs',i,k,l,j,block(a_blk_no)%xu(loc_x),block(a_blk_no)%yu(loc_y),block(b_blk_no)%xu(j),block(b_blk_no)%yu(l)
!1221            FORMAT(A4,4I4,4F10.5)
!                write(*,1222) bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans
!1222            FORMAT(' ',5F10.5)
      !        else
      !         if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
      !         aa1=(block(a_blk_no)%xu(loc_x-1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y-1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=aa1
      !         elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
      !         aa2=(block(a_blk_no)%xu(loc_x-1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y+1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=dmin1(aamin,aa2)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
      !         aa3=(block(a_blk_no)%xu(loc_x+1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y+1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=dmin1(aamin,aa3)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
      !         aa4=(block(a_blk_no)%xu(loc_x+1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y-1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=dmin1(aamin,aa4)
      !         endif

      !         if (aamin .eq. aa1)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x,loc_y)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x-2,loc_y-1)
      !         elseif(aamin .eq. aa2)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x,loc_y+1)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x-2,loc_y+1)
      !         elseif(aamin .eq. aa3)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x+1,loc_y+1)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x,loc_y+1)
      !         elseif (aamin .eq. aa4)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x+1,loc_y)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x,loc_y-1)
      !         endif
      ! endif
        endif
        ENDDO
        ENDDO
        ENDDO
        !endif
        ENDDO
        ENDDO
        ENDDO

              ! tar_x=j
              ! tar_y=2
              ! loc_x=intfr(g)%ux_interface_det(1,i)
              ! loc_y=intfr(g)%uy_interface_det(1,1)
              !
!!!!!p
              ! bl_intp_valx=block(b_blk_no)%xu(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yu(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xu(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yu(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
              !


              ! tar_x=j
              ! tar_y=block(b_blk_no)%jtc_en
              ! coarse_indx=intfr(g)%counteryu
              ! loc_x=intfr(g)%ux_interface_det(1,i)
              ! loc_y=intfr(g)%uy_interface_det(1,coarse_indx)
              !
!!!!!p
              ! bl_intp_valx=block(b_blk_no)%xu(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yu(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xu(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yu(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans


              ! tar_x=j
              ! tar_y=block(b_blk_no)%jtc_en-1
              ! coarse_indx=intfr(g)%counteryu
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,coarse_indx)
              !
!!!!!p
              ! bl_intp_valx=block(b_blk_no)%xu(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yu(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xu(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yu(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans

              ! ENDDO
              ! ENDDO
!!!!!v
        DO k=2,intfr(g)%counterzv-1
        DO j=2,intfr(g)%counteryv-1
        DO i=2,intfr(g)%counterxv-1
        varx1=intfr(g)%vx_interface_det(2,i)
        varx2=intfr(g)%vx_interface_det(3,i)
        vary1=intfr(g)%vy_interface_det(2,j)
        vary2=intfr(g)%vy_interface_det(3,j)
        varz1=intfr(g)%vz_interface_det(2,k)
        varz2=intfr(g)%vz_interface_det(3,k)

        loc_x=intfr(g)%vx_interface_det(1,i)
        loc_y=intfr(g)%vy_interface_det(1,j)
        loc_z=intfr(g)%vz_interface_det(1,k)
     !  if (k .eq. 1 .or. i .eq.1 .or. k .eq. intfr(g)%counteryv .or. i .eq. intfr(g)%counterxv)then
     !

     !  DO l=vary1,vary2
     !  DO j=varx1,varx2
     !          tar_x=j
     !          !tar_y=1
     !          tar_y=l
     !          block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x,loc_y-1)
     ! ENDDO
     ! ENDDO
     !
     !  else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
               !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) ==0)then
      ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then


!!!!!p
                bl_intp_valx=block(b_blk_no)%xv(tar_x)
                bl_intp_valy=block(b_blk_no)%yv(tar_y)
                bl_intp_valz=block(b_blk_no)%zv(tar_z)
                bl_intp_x1=block(a_blk_no)%xv(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yv(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zv(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zv(loc_z+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
                bl_intp_f1=block(a_blk_no)%v(loc_x-1,loc_y-2,loc_z-1)
                bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y-2,loc_z-1)
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y,loc_z-1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%v(loc_x-1,loc_y-2,loc_z+1)
                bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y-2,loc_z+1)
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y,loc_z+1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y,loc_z+1)
                !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans

      !         else
      !         if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
      !         aa1=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=aa1
      !         elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
      !         aa2=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa2)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
      !         aa3=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa3)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
      !         aa4=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa4)
      !         endif

      !         if (aamin .eq. aa1)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y-2)
      !         elseif(aamin .eq. aa2)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y)
      !         elseif(aamin .eq. aa3)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         elseif (aamin .eq. aa4)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y-2)
      !         endif

      ! endif
        endif
        ENDDO
        ENDDO
        ENDDO
       ! endif
        ENDDO
        ENDDO
        ENDDO
!!!!!w
        DO k=2,intfr(g)%counterzw-1
        DO j=2,intfr(g)%counteryw-1
        DO i=2,intfr(g)%counterxw-1
        varx1=intfr(g)%wx_interface_det(2,i)
        varx2=intfr(g)%wx_interface_det(3,i)
        vary1=intfr(g)%wy_interface_det(2,j)
        vary2=intfr(g)%wy_interface_det(3,j)
        varz1=intfr(g)%wz_interface_det(2,k)
        varz2=intfr(g)%wz_interface_det(3,k)

        loc_x=intfr(g)%wx_interface_det(1,i)
        loc_y=intfr(g)%wy_interface_det(1,j)
        loc_z=intfr(g)%wz_interface_det(1,k)
     !  if (k .eq. 1 .or. i .eq.1 .or. k .eq. intfr(g)%counteryv .or. i .eq. intfr(g)%counterxv)then
     !

     !  DO l=vary1,vary2
     !  DO j=varx1,varx2
     !          tar_x=j
     !          !tar_y=1
     !          tar_y=l
     !          block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x,loc_y-1)
     ! ENDDO
     ! ENDDO
     !
     !  else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
               !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) ==0)then
      ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then


!!!!!p
                bl_intp_valx=block(b_blk_no)%xw(tar_x)
                bl_intp_valy=block(b_blk_no)%yw(tar_y)
                bl_intp_valz=block(b_blk_no)%zw(tar_z)
                bl_intp_x1=block(a_blk_no)%xw(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xw(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yw(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yw(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zw(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zw(loc_z+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
                bl_intp_f1=block(a_blk_no)%w(loc_x-1,loc_y-1,loc_z-2)
                bl_intp_f2=block(a_blk_no)%w(loc_x+1,loc_y-1,loc_z-2)
                bl_intp_f3=block(a_blk_no)%w(loc_x+1,loc_y+1,loc_z-2)
                bl_intp_f4=block(a_blk_no)%w(loc_x-1,loc_y+1,loc_z-2)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%w(loc_x-1,loc_y-1,loc_z)
                bl_intp_f2=block(a_blk_no)%w(loc_x+1,loc_y-1,loc_z)
                bl_intp_f3=block(a_blk_no)%w(loc_x+1,loc_y+1,loc_z)
                bl_intp_f4=block(a_blk_no)%w(loc_x-1,loc_y+1,loc_z)
                !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans

      !         else
      !         if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
      !         aa1=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=aa1
      !         elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
      !         aa2=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa2)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
      !         aa3=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa3)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
      !         aa4=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa4)
      !         endif

      !         if (aamin .eq. aa1)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y-2)
      !         elseif(aamin .eq. aa2)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y)
      !         elseif(aamin .eq. aa3)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         elseif (aamin .eq. aa4)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y-2)
      !         endif

      ! endif
        endif
        ENDDO
        ENDDO
        ENDDO
       ! endif
        ENDDO
        ENDDO
        ENDDO
               !tar_x=j
               !tar_y=2
               !loc_x=intfr(g)%vx_interface_det(1,i)
               !loc_y=intfr(g)%vy_interface_det(1,1)

!!!!!p
               !bl_intp_valx=block(b_blk_no)%xv(tar_x)
               !bl_intp_valy=block(b_blk_no)%yv(tar_y)
               !bl_intp_x1=block(a_blk_no)%xv(loc_x)
               !bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
               !bl_intp_y1=block(a_blk_no)%yv(loc_y)
               !bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
               !
               !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
               !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
               !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
               !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
               !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
               !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
               !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
               !bl_intp_num= bl_intp_first_term + bl_intp_second_term

               !bl_interp_ans=(bl_intp_num/bl_intp_deno)
               !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
               !


               !tar_x=j
               !tar_y=block(b_blk_no)%jtc_en
               !coarse_indx=intfr(g)%counteryv
               !loc_x=intfr(g)%vx_interface_det(1,i)
               !loc_y=intfr(g)%vy_interface_det(1,coarse_indx)

!!!!!p
               !bl_intp_valx=block(b_blk_no)%xv(tar_x)
               !bl_intp_valy=block(b_blk_no)%yv(tar_y)
               !bl_intp_x1=block(a_blk_no)%xv(loc_x)
               !bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
               !bl_intp_y1=block(a_blk_no)%yv(loc_y)
               !bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
               !
               !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
               !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
               !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
               !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
               !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
               !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
               !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
               !bl_intp_num= bl_intp_first_term + bl_intp_second_term

               !bl_interp_ans=(bl_intp_num/bl_intp_deno)
               !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans


               !tar_x=j
               !tar_y=block(b_blk_no)%jtc_en-1
               !coarse_indx=intfr(g)%counteryv
               !loc_x=intfr(g)%vx_interface_det(1,i)
               !loc_y=intfr(g)%vy_interface_det(1,coarse_indx)

!!!!!p
               !bl_intp_valx=block(b_blk_no)%xv(tar_x)
               !bl_intp_valy=block(b_blk_no)%yv(tar_y)
               !bl_intp_x1=block(a_blk_no)%xv(loc_x)
               !bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
               !bl_intp_y1=block(a_blk_no)%yv(loc_y)
               !bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
               !
               !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
               !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
               !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
               !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
               !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
               !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
               !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
               !bl_intp_num= bl_intp_first_term + bl_intp_second_term

               !bl_interp_ans=(bl_intp_num/bl_intp_deno)
               !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans

               !ENDDO
               !ENDDO



!       DO i=1,intfr(g)%counteryp
!       varx1=intfr(g)%py_interface_det(2,i)
!       varx2=intfr(g)%py_interface_det(3,i)


!       DO j=varx1,varx2
!               tar_y=j
!               tar_x=1
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans
!

!               tar_y=j
!               tar_x=2
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans
!


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en
!               coarse_indx=intfr(g)%counterxp
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en-1
!               coarse_indx=intfr(g)%counterxp
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans

!               ENDDO
!               ENDDO
!!!!!u
!       DO i=1,intfr(g)%counteryu
!       varx1=intfr(g)%uy_interface_det(2,i)
!       varx2=intfr(g)%uy_interface_det(3,i)
!       DO j=varx1,varx2
!               tar_y=j
!               tar_x=1
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
!

!               tar_y=j
!               tar_x=2
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
!


!               tar_y=j
!               tar_x=block(b_blk_no)%jtc_en
!               coarse_indx=intfr(g)%counterxu
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en-1
!               coarse_indx=intfr(g)%counterxu
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans

!               ENDDO
!               ENDDO
!!!!!v
!       DO i=1,intfr(g)%counteryv
!       varx1=intfr(g)%vy_interface_det(2,i)
!       varx2=intfr(g)%vy_interface_det(3,i)
!       DO j=varx1,varx2
!               tar_y=j
!               tar_x=1
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
!

!               tar_y=j
!               tar_x=2
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
!


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en
!               coarse_indx=intfr(g)%counterxv
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en-1
!               coarse_indx=intfr(g)%counterxv
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans

!               ENDDO
!               ENDDO

        ENDDO


        end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE fineUpdate_mv(g)
        USE global
        IMPLICIT NONE
        REAL (KIND=8) :: bl_intp_valx,bl_intp_valy,bl_intp_valz,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (KIND=8) :: bl_intp_z1,bl_intp_z2
        REAL (KIND=8) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (KIND=8) :: bl_interp_ans, bl_interp_ans1, bl_interp_ans2
        REAL (KIND=8) :: aa1,aa2,aa3,aa4,aamin

        INTEGER(KIND=8) :: i,j,k,q,s, varx1,varx2, vary1, vary2, l,tar_x, tar_y, loc_x, loc_y, coarse_indx,a_blk_no, b_blk_no
        INTEGER(KIND=8) :: varz1,varz2, tar_z, loc_z
        INTEGER(KIND=8), INTENT(IN) :: g
        CHARACTER(len=150) :: filename1





        !DO g=1,intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk

        DO k=2,intfr(g)%counterzp-1
        DO j=2,intfr(g)%counteryp-1
        DO i=2,intfr(g)%counterxp-1
        varx1=intfr(g)%px_interface_det(2,i)
        varx2=intfr(g)%px_interface_det(3,i)
        vary1=intfr(g)%py_interface_det(2,j)
        vary2=intfr(g)%py_interface_det(3,j)
        varz1=intfr(g)%pz_interface_det(2,k)
        varz2=intfr(g)%pz_interface_det(3,k)

        loc_x=intfr(g)%px_interface_det(1,i)
        loc_y=intfr(g)%py_interface_det(1,j)
        loc_z=intfr(g)%pz_interface_det(1,k)

      ! if ( i .eq. 1 .or. i .eq. intfr(g)%counterxp .or. k .eq. 1 .or. k .eq. intfr(g)%counteryp)then

      ! DO l=vary1,vary2
      ! DO j=varx1,varx2
      !         tar_x=j
      !         !tar_y=1
      !         tar_y=l
      !
      !         block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x,loc_y)
      ! ENDDO
      ! ENDDO
      !else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
               !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               !if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
       !if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
       !      block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and.  &
       !      block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
       !      block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then

 !!!!p1
                bl_intp_valx=block(b_blk_no)%xp(tar_x)
                bl_intp_valy=block(b_blk_no)%yp(tar_y)
                bl_intp_valz=block(b_blk_no)%zp(tar_z)
                bl_intp_x1=block(a_blk_no)%xp(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yp(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zp(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zp(loc_z+1)
                bl_intp_f1=block(a_blk_no)%p(loc_x-1,loc_y-1,loc_z-1)
                bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1,loc_z-1)
                bl_intp_f4=block(a_blk_no)%p(loc_x-1,loc_y+1,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%p(loc_x-1,loc_y-1,loc_z+1)
                bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%p(loc_x-1,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans

              ! else

              ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
              ! aa1=(block(a_blk_no)%xp(loc_x-1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y-1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=aa1
              ! elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
              ! aa2=(block(a_blk_no)%xp(loc_x-1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y+1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=dmin1(aamin,aa2)
              ! elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
              ! aa3=(block(a_blk_no)%xp(loc_x+1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y+1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=dmin1(aamin,aa3)
              ! elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
              ! aa4=(block(a_blk_no)%xp(loc_x+1)-block(b_blk_no)%xp(tar_x))**2+ (block(a_blk_no)%yp(loc_y-1)-block(b_blk_no)%yp(tar_y))**2
              ! aamin=dmin1(aamin,aa4)
              ! endif

              ! if (aamin .eq. aa1)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x-1,loc_y-1)
              ! elseif(aamin .eq. aa2)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x-1,loc_y+1)
              ! elseif(aamin .eq. aa3)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! elseif (aamin .eq. aa4)then
              ! block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_no)%p(loc_x+1,loc_y-1)
              ! endif

              ! endif
                endif
                ENDDO
                ENDDO
                ENDDO
             !   endif
                ENDDO
                ENDDO
                ENDDO


       !DO l=vary1,vary2
       !DO j=varx1,varx2
       !        tar_x=j
       !        !tar_y=1
       !        tar_y=l
       !       !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
       !       if (block(b_blk_no)%cell(tar_x,tar_y) .eq.0)then
!!!!!p1
       !     !! bl_intp_valx=block(b_blk_no)%xp(tar_x)
       !     !! bl_intp_valy=block(b_blk_no)%yp(tar_y)
       !     !! bl_intp_x1=block(a_blk_no)%xp(loc_x)
       !     !! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
       !     !! bl_intp_y1=block(a_blk_no)%yp(loc_y)
       !     !! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
       !     !! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
       !     !! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
       !     !! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
       !     !! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
       !     !!
       !     !! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
       !     !! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
       !     !! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
       !     !! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
       !     !! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
       !     !! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
       !     !! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
       !     !! bl_intp_num= bl_intp_first_term + bl_intp_second_term
       !
       !     !! bl_interp_ans=(bl_intp_num/bl_intp_deno)
       !     !   block(b_blk_no)%p(tar_x,tar_y)=block(a_blk_
       !      endif
              ! tar_x=j
              ! !tar_y=2
              ! tar_y=k
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,1)
              !
!!!!!p2
              ! bl_intp_valx=block(b_blk_no)%xp(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yp(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xp(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yp(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans
              !


              ! tar_x=j
              ! !tar_y=block(b_blk_no)%jtc_en
              ! tar_y=k
              ! coarse_indx=intfr(g)%counteryp
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,coarse_indx)
              !
!!!!!pny
              ! bl_intp_valx=block(b_blk_no)%xp(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yp(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xp(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yp(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans


              ! tar_x=j
              ! tar_y=block(b_blk_no)%jtc_en-1
              ! coarse_indx=intfr(g)%counteryp
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,coarse_indx)
              !
!!!!!pny-1
              ! bl_intp_valx=block(b_blk_no)%xp(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yp(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xp(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yp(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans

!!!!!u
        DO k=2,intfr(g)%counterzu-1
        DO j=2,intfr(g)%counteryu-1
        DO i=2,intfr(g)%counterxu-1
        varx1=intfr(g)%ux_interface_det(2,i)
        varx2=intfr(g)%ux_interface_det(3,i)
        vary1=intfr(g)%uy_interface_det(2,j)
        vary2=intfr(g)%uy_interface_det(3,j)
        varz1=intfr(g)%uz_interface_det(2,k)
        varz2=intfr(g)%uz_interface_det(3,k)


        loc_x=intfr(g)%ux_interface_det(1,i)
        loc_y=intfr(g)%uy_interface_det(1,j)
        loc_z=intfr(g)%uz_interface_det(1,k)
      ! if (k .eq. 1 .or. i .eq.1 .or. k .eq. intfr(g)%counteryu .or. i .eq. intfr(g)%counterxu)then
      !

      ! DO l=vary1,vary2
      ! DO j=varx1,varx2
      !         tar_x=j
      !         !tar_y=1
      !         tar_y=l
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x-1,loc_y)
      !ENDDO
      !ENDDO
      !
      ! else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
              ! if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
      ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then
      !
!!!!!p

!        write(filename1,1)
!1        format('test.dat')
!        OPEN(919,FILE=filename1,access='append',status='unknown')

!                 write(919,1222)ita,i,loc_x,loc_y,loc_z,tar_x,tar_y,tar_z,block(a_blk_no)%xu(loc_x-1),block(a_blk_no)%xu(loc_x+1),block(b_blk_no)%xu(tar_x)
! 1222            FORMAT(8I6,3F10.5)
!                print*,'!!!',i,loc_x,tar_x
                bl_intp_valx=block(b_blk_no)%xu(tar_x)
                bl_intp_valy=block(b_blk_no)%yu(tar_y)
                bl_intp_valz=block(b_blk_no)%zu(tar_z)
                bl_intp_x1=block(a_blk_no)%xu(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yu(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zu(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zu(loc_z+1)
               !bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
                bl_intp_f1=block(a_blk_no)%u(loc_x-2,loc_y-1,loc_z-1)
                bl_intp_f2=block(a_blk_no)%u(loc_x,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%u(loc_x,loc_y+1,loc_z-1)
                bl_intp_f4=block(a_blk_no)%u(loc_x-2,loc_y+1,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%u(loc_x-2,loc_y-1,loc_z+1)
                bl_intp_f2=block(a_blk_no)%u(loc_x,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%u(loc_x,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%u(loc_x-2,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                !block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans
!                write(*,1221)'locs',i,k,l,j,block(a_blk_no)%xu(loc_x),block(a_blk_no)%yu(loc_y),block(b_blk_no)%xu(j),block(b_blk_no)%yu(l)
!1221            FORMAT(A4,4I4,4F10.5)
!                write(*,1222) bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4,bl_interp_ans
!1222            FORMAT(' ',5F10.5)
      !        else
      !         if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
      !         aa1=(block(a_blk_no)%xu(loc_x-1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y-1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=aa1
      !         elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
      !         aa2=(block(a_blk_no)%xu(loc_x-1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y+1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=dmin1(aamin,aa2)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
      !         aa3=(block(a_blk_no)%xu(loc_x+1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y+1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=dmin1(aamin,aa3)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
      !         aa4=(block(a_blk_no)%xu(loc_x+1)-block(b_blk_no)%xu(tar_x))**2+ (block(a_blk_no)%yu(loc_y-1)-block(b_blk_no)%yu(tar_y))**2
      !         aamin=dmin1(aamin,aa4)
      !         endif

      !         if (aamin .eq. aa1)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x,loc_y)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x-2,loc_y-1)
      !         elseif(aamin .eq. aa2)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x,loc_y+1)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x-2,loc_y+1)
      !         elseif(aamin .eq. aa3)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x+1,loc_y+1)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x,loc_y+1)
      !         elseif (aamin .eq. aa4)then
      !         !block(b_blk_no)%u(tar_x,tar_y)=block(a_blk_no)%u(loc_x+1,loc_y)
      !         block(b_blk_no)%u(tar_x-1,tar_y)=block(a_blk_no)%u(loc_x,loc_y-1)
      !         endif
      ! endif
        endif
        ENDDO
        ENDDO
        ENDDO
        !endif
        ENDDO
        ENDDO
        ENDDO

              ! tar_x=j
              ! tar_y=2
              ! loc_x=intfr(g)%ux_interface_det(1,i)
              ! loc_y=intfr(g)%uy_interface_det(1,1)
              !
!!!!!p
              ! bl_intp_valx=block(b_blk_no)%xu(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yu(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xu(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yu(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
              !


              ! tar_x=j
              ! tar_y=block(b_blk_no)%jtc_en
              ! coarse_indx=intfr(g)%counteryu
              ! loc_x=intfr(g)%ux_interface_det(1,i)
              ! loc_y=intfr(g)%uy_interface_det(1,coarse_indx)
              !
!!!!!p
              ! bl_intp_valx=block(b_blk_no)%xu(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yu(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xu(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yu(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans


              ! tar_x=j
              ! tar_y=block(b_blk_no)%jtc_en-1
              ! coarse_indx=intfr(g)%counteryu
              ! loc_x=intfr(g)%px_interface_det(1,i)
              ! loc_y=intfr(g)%py_interface_det(1,coarse_indx)
              !
!!!!!p
              ! bl_intp_valx=block(b_blk_no)%xu(tar_x)
              ! bl_intp_valy=block(b_blk_no)%yu(tar_y)
              ! bl_intp_x1=block(a_blk_no)%xu(loc_x)
              ! bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
              ! bl_intp_y1=block(a_blk_no)%yu(loc_y)
              ! bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
              ! bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
              ! bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
              ! bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
              ! bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
              !
              ! bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
              ! bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
              ! bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
              ! bl_intp_yty=(bl_intp_y2-bl_intp_valy)
              ! bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
              ! bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
              ! bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
              ! bl_intp_num= bl_intp_first_term + bl_intp_second_term

              ! bl_interp_ans=(bl_intp_num/bl_intp_deno)
              ! block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans

              ! ENDDO
              ! ENDDO
!!!!!v
        DO k=2,intfr(g)%counterzv-1
        DO j=2,intfr(g)%counteryv-1
        DO i=2,intfr(g)%counterxv-1
        varx1=intfr(g)%vx_interface_det(2,i)
        varx2=intfr(g)%vx_interface_det(3,i)
        vary1=intfr(g)%vy_interface_det(2,j)
        vary2=intfr(g)%vy_interface_det(3,j)
        varz1=intfr(g)%vz_interface_det(2,k)
        varz2=intfr(g)%vz_interface_det(3,k)

        loc_x=intfr(g)%vx_interface_det(1,i)
        loc_y=intfr(g)%vy_interface_det(1,j)
        loc_z=intfr(g)%vz_interface_det(1,k)
     !  if (k .eq. 1 .or. i .eq.1 .or. k .eq. intfr(g)%counteryv .or. i .eq. intfr(g)%counterxv)then
     !

     !  DO l=vary1,vary2
     !  DO j=varx1,varx2
     !          tar_x=j
     !          !tar_y=1
     !          tar_y=l
     !          block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x,loc_y-1)
     ! ENDDO
     ! ENDDO
     !
     !  else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
               !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
      ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then


!!!!!p
                bl_intp_valx=block(b_blk_no)%xv(tar_x)
                bl_intp_valy=block(b_blk_no)%yv(tar_y)
                bl_intp_valz=block(b_blk_no)%zv(tar_z)
                bl_intp_x1=block(a_blk_no)%xv(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yv(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zv(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zv(loc_z+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
                bl_intp_f1=block(a_blk_no)%v(loc_x-1,loc_y-2,loc_z-1)
                bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y-2,loc_z-1)
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y,loc_z-1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%v(loc_x-1,loc_y-2,loc_z+1)
                bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y-2,loc_z+1)
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y,loc_z+1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans

      !         else
      !         if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
      !         aa1=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=aa1
      !         elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
      !         aa2=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa2)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
      !         aa3=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa3)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
      !         aa4=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa4)
      !         endif

      !         if (aamin .eq. aa1)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y-2)
      !         elseif(aamin .eq. aa2)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y)
      !         elseif(aamin .eq. aa3)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         elseif (aamin .eq. aa4)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y-2)
      !         endif

      ! endif
        endif
        ENDDO
        ENDDO
        ENDDO
       ! endif
        ENDDO
        ENDDO
        ENDDO
!!!!!w
        DO k=2,intfr(g)%counterzw-1
        DO j=2,intfr(g)%counteryw-1
        DO i=2,intfr(g)%counterxw-1
        varx1=intfr(g)%wx_interface_det(2,i)
        varx2=intfr(g)%wx_interface_det(3,i)
        vary1=intfr(g)%wy_interface_det(2,j)
        vary2=intfr(g)%wy_interface_det(3,j)
        varz1=intfr(g)%wz_interface_det(2,k)
        varz2=intfr(g)%wz_interface_det(3,k)

        loc_x=intfr(g)%wx_interface_det(1,i)
        loc_y=intfr(g)%wy_interface_det(1,j)
        loc_z=intfr(g)%wz_interface_det(1,k)
     !  if (k .eq. 1 .or. i .eq.1 .or. k .eq. intfr(g)%counteryv .or. i .eq. intfr(g)%counterxv)then
     !

     !  DO l=vary1,vary2
     !  DO j=varx1,varx2
     !          tar_x=j
     !          !tar_y=1
     !          tar_y=l
     !          block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x,loc_y-1)
     ! ENDDO
     ! ENDDO
     !
     !  else
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                !tar_y=1
                tar_y=l
                tar_z=q
               !if (block(a_blk_no)%cell(loc_x,loc_y) .eq.0)then
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
      ! if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y-1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x+1,loc_y+1) .ne. 1 .and. &
      !       block(a_blk_no)%cell(loc_x-1,loc_y+1) .ne. 1 ) then


!!!!!p
                bl_intp_valx=block(b_blk_no)%xw(tar_x)
                bl_intp_valy=block(b_blk_no)%yw(tar_y)
                bl_intp_valz=block(b_blk_no)%zw(tar_z)
                bl_intp_x1=block(a_blk_no)%xw(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xw(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yw(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yw(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zw(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zw(loc_z+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
                bl_intp_f1=block(a_blk_no)%w(loc_x-1,loc_y-1,loc_z-2)
                bl_intp_f2=block(a_blk_no)%w(loc_x+1,loc_y-1,loc_z-2)
                bl_intp_f3=block(a_blk_no)%w(loc_x+1,loc_y+1,loc_z-2)
                bl_intp_f4=block(a_blk_no)%w(loc_x-1,loc_y+1,loc_z-2)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%w(loc_x-1,loc_y-1,loc_z)
                bl_intp_f2=block(a_blk_no)%w(loc_x+1,loc_y-1,loc_z)
                bl_intp_f3=block(a_blk_no)%w(loc_x+1,loc_y+1,loc_z)
                bl_intp_f4=block(a_blk_no)%w(loc_x-1,loc_y+1,loc_z)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans

      !         else
      !         if ( block(a_blk_no)%cell(loc_x-1,loc_y-1) .eq. 0 )then
      !         aa1=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=aa1
      !         elseif( block(a_blk_no)%cell(loc_x-1,loc_y+1) .eq. 0)then
      !         aa2=(block(a_blk_no)%xv(loc_x-1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa2)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y+1) .eq. 0 )then
      !         aa3=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y+1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa3)
      !         elseif( block(a_blk_no)%cell(loc_x+1,loc_y-1) .eq. 0)then
      !         aa4=(block(a_blk_no)%xv(loc_x+1)-block(b_blk_no)%xv(tar_x))**2+ (block(a_blk_no)%yv(loc_y-1)-block(b_blk_no)%yv(tar_y))**2
      !         aamin=dmin1(aamin,aa4)
      !         endif

      !         if (aamin .eq. aa1)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y-2)
      !         elseif(aamin .eq. aa2)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x-1,loc_y)
      !         elseif(aamin .eq. aa3)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y+1)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         elseif (aamin .eq. aa4)then
      !         !block(b_blk_no)%v(tar_x,tar_y)=block(a_blk_no)%v(loc_x+1,loc_y)
      !         block(b_blk_no)%v(tar_x,tar_y-1)=block(a_blk_no)%v(loc_x+1,loc_y-2)
      !         endif

      ! endif
        endif
        ENDDO
        ENDDO
        ENDDO
       ! endif
        ENDDO
        ENDDO
        ENDDO
               !tar_x=j
               !tar_y=2
               !loc_x=intfr(g)%vx_interface_det(1,i)
               !loc_y=intfr(g)%vy_interface_det(1,1)

!!!!!p
               !bl_intp_valx=block(b_blk_no)%xv(tar_x)
               !bl_intp_valy=block(b_blk_no)%yv(tar_y)
               !bl_intp_x1=block(a_blk_no)%xv(loc_x)
               !bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
               !bl_intp_y1=block(a_blk_no)%yv(loc_y)
               !bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
               !
               !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
               !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
               !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
               !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
               !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
               !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
               !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
               !bl_intp_num= bl_intp_first_term + bl_intp_second_term

               !bl_interp_ans=(bl_intp_num/bl_intp_deno)
               !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
               !


               !tar_x=j
               !tar_y=block(b_blk_no)%jtc_en
               !coarse_indx=intfr(g)%counteryv
               !loc_x=intfr(g)%vx_interface_det(1,i)
               !loc_y=intfr(g)%vy_interface_det(1,coarse_indx)

!!!!!p
               !bl_intp_valx=block(b_blk_no)%xv(tar_x)
               !bl_intp_valy=block(b_blk_no)%yv(tar_y)
               !bl_intp_x1=block(a_blk_no)%xv(loc_x)
               !bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
               !bl_intp_y1=block(a_blk_no)%yv(loc_y)
               !bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
               !
               !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
               !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
               !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
               !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
               !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
               !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
               !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
               !bl_intp_num= bl_intp_first_term + bl_intp_second_term

               !bl_interp_ans=(bl_intp_num/bl_intp_deno)
               !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans


               !tar_x=j
               !tar_y=block(b_blk_no)%jtc_en-1
               !coarse_indx=intfr(g)%counteryv
               !loc_x=intfr(g)%vx_interface_det(1,i)
               !loc_y=intfr(g)%vy_interface_det(1,coarse_indx)

!!!!!p
               !bl_intp_valx=block(b_blk_no)%xv(tar_x)
               !bl_intp_valy=block(b_blk_no)%yv(tar_y)
               !bl_intp_x1=block(a_blk_no)%xv(loc_x)
               !bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
               !bl_intp_y1=block(a_blk_no)%yv(loc_y)
               !bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
               !bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
               !bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
               !bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
               !bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
               !
               !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
               !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
               !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
               !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
               !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
               !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
               !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
               !bl_intp_num= bl_intp_first_term + bl_intp_second_term

               !bl_interp_ans=(bl_intp_num/bl_intp_deno)
               !block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans

               !ENDDO
               !ENDDO



!       DO i=1,intfr(g)%counteryp
!       varx1=intfr(g)%py_interface_det(2,i)
!       varx2=intfr(g)%py_interface_det(3,i)


!       DO j=varx1,varx2
!               tar_y=j
!               tar_x=1
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans
!

!               tar_y=j
!               tar_x=2
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans
!


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en
!               coarse_indx=intfr(g)%counterxp
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en-1
!               coarse_indx=intfr(g)%counterxp
!               loc_y=intfr(g)%py_interface_det(1,i)
!               loc_x=intfr(g)%px_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xp(tar_x)
!               bl_intp_valy=block(b_blk_no)%yp(tar_y)
!               bl_intp_x1=block(a_blk_no)%xp(loc_x)
!               bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yp(loc_y)
!               bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%p(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%p(tar_x,tar_y)=bl_interp_ans

!               ENDDO
!               ENDDO
!!!!!u
!       DO i=1,intfr(g)%counteryu
!       varx1=intfr(g)%uy_interface_det(2,i)
!       varx2=intfr(g)%uy_interface_det(3,i)
!       DO j=varx1,varx2
!               tar_y=j
!               tar_x=1
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
!

!               tar_y=j
!               tar_x=2
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans
!


!               tar_y=j
!               tar_x=block(b_blk_no)%jtc_en
!               coarse_indx=intfr(g)%counterxu
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en-1
!               coarse_indx=intfr(g)%counterxu
!               loc_y=intfr(g)%uy_interface_det(1,i)
!               loc_x=intfr(g)%ux_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xu(tar_x)
!               bl_intp_valy=block(b_blk_no)%yu(tar_y)
!               bl_intp_x1=block(a_blk_no)%xu(loc_x)
!               bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yu(loc_y)
!               bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%u(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%u(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%u(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%u(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%u(tar_x,tar_y)=bl_interp_ans

!               ENDDO
!               ENDDO
!!!!!v
!       DO i=1,intfr(g)%counteryv
!       varx1=intfr(g)%vy_interface_det(2,i)
!       varx2=intfr(g)%vy_interface_det(3,i)
!       DO j=varx1,varx2
!               tar_y=j
!               tar_x=1
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
!

!               tar_y=j
!               tar_x=2
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,1)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans
!


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en
!               coarse_indx=intfr(g)%counterxv
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans


!               tar_y=j
!               tar_x=block(b_blk_no)%itc_en-1
!               coarse_indx=intfr(g)%counterxv
!               loc_y=intfr(g)%vy_interface_det(1,i)
!               loc_x=intfr(g)%vx_interface_det(1,coarse_indx)
!
!!!!!p
!               bl_intp_valx=block(b_blk_no)%xv(tar_x)
!               bl_intp_valy=block(b_blk_no)%yv(tar_y)
!               bl_intp_x1=block(a_blk_no)%xv(loc_x)
!               bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(a_blk_no)%yv(loc_y)
!               bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
!               bl_intp_f1=block(a_blk_no)%v(loc_x,loc_y)
!               bl_intp_f2=block(a_blk_no)%v(loc_x+1,loc_y)
!               bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y+1)
!               bl_intp_f4=block(a_blk_no)%v(loc_x,loc_y+1)
!
!               bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!               bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!               bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!               bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!               bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term
!
!               bl_interp_ans=(bl_intp_num/bl_intp_deno)
!               block(b_blk_no)%v(tar_x,tar_y)=bl_interp_ans

!               ENDDO
!               ENDDO



        end subroutine

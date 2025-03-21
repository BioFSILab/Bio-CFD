        subroutine coarseUpdate
        use global
        IMPLICIT NONE
        REAL (KIND=8) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (KIND=8) :: bl_intp_valz,bl_intp_z1,bl_intp_z2,bl_intp_f5,bl_intp_f6,bl_intp_f7,bl_intp_f8
        REAL (KIND=8) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (KIND=8) :: bl_interp_ans1, bl_interp_ans2, bl_interp_ansf
        REAL (KIND=8) :: bl_interp_ans

        INTEGER :: i,j,k,g, varx1,varx2, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz
       !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
       !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
       !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
       !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
       !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
       !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
       !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
       !bl_intp_num= bl_intp_first_term + bl_intp_second_term
       !
       !bl_interp_ans=(bl_intp_num/bl_intp_deno)



        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        !counterx= (intfr(g)%xintf_end-intfr(g)%xintf_start +1)/factor
        !countery= (intfr(g)%yintf_end-intfr(g)%yintf_start +1)/factor
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxp-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryp-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzp-block(b_blk_no)%cintp
      !do i=2, intfr(g)%counterxp-1
      ! do j=2,intfr(g)%counteryp-1
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
      !do k=2,intfr(g)%counterzp-1
      !do j=2,intfr(g)%counteryp-1
      !do i=2,intfr(g)%counterxp-1
                tar_x=intfr(g)%px_interface_det(1,i)
                tar_y=intfr(g)%py_interface_det(1,j)
                tar_z=intfr(g)%pz_interface_det(1,k)
                loc_x=(intfr(g)%px_interface_det(2,i)+intfr(g)%px_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%py_interface_det(2,j)+intfr(g)%py_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%pz_interface_det(2,k)+intfr(g)%pz_interface_det(3,k)-1)/2
              !!  if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
              !!if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
 !!!!p
                bl_intp_valx=block(a_blk_no)%xp(tar_x)
                bl_intp_valy=block(a_blk_no)%yp(tar_y)
                bl_intp_valz=block(a_blk_no)%zp(tar_z)
                bl_intp_x1=block(b_blk_no)%xp(loc_x)
                bl_intp_x2=block(b_blk_no)%xp(loc_x+1)
                bl_intp_y1=block(b_blk_no)%yp(loc_y)
                bl_intp_y2=block(b_blk_no)%yp(loc_y+1)
                bl_intp_z1=block(b_blk_no)%zp(loc_z)
                bl_intp_z2=block(b_blk_no)%zp(loc_z+1)
                bl_intp_f1=block(b_blk_no)%p(loc_x,loc_y,loc_z)
                bl_intp_f2=block(b_blk_no)%p(loc_x+1,loc_y,loc_z)
                bl_intp_f3=block(b_blk_no)%p(loc_x+1,loc_y+1,loc_z)
                bl_intp_f4=block(b_blk_no)%p(loc_x,loc_y+1,loc_z)
               !bl_intp_f5=block(b_blk_no)%p(loc_x,loc_y,loc_z+1)
               !bl_intp_f6=block(b_blk_no)%p(loc_x+1,loc_y,loc_z+1)
               !bl_intp_f7=block(b_blk_no)%p(loc_x+1,loc_y+1,loc_z+1)
               !bl_intp_f8=block(b_blk_no)%p(loc_x,loc_y+1,loc_z+1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)

                bl_intp_f1=block(b_blk_no)%p(loc_x,loc_y,loc_z+1)
                bl_intp_f2=block(b_blk_no)%p(loc_x+1,loc_y,loc_z+1)
                bl_intp_f3=block(b_blk_no)%p(loc_x+1,loc_y+1,loc_z+1)
                bl_intp_f4=block(b_blk_no)%p(loc_x,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))

                block(a_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans
 !!!!u
               !! endif


        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxu-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryu-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzu-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
      !do k=2,intfr(g)%counterzu-1
      !do j=2,intfr(g)%counteryu-1
      !do i=2,intfr(g)%counterxu-1
      !do i=2, intfr(g)%counterxu-1
      ! do j=2,intfr(g)%counteryu-1
                tar_x=intfr(g)%ux_interface_det(1,i)
                tar_y=intfr(g)%uy_interface_det(1,j)
                tar_z=intfr(g)%uz_interface_det(1,k)
                !loc_x=(intfr(g)%ux_interface_det(2,i)+intfr(g)%ux_interface_det(3,i)-1)/2
                loc_x=intfr(g)%ux_interface_det(3,i)
                loc_y=(intfr(g)%uy_interface_det(2,j)+intfr(g)%uy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%uz_interface_det(2,k)+intfr(g)%uz_interface_det(3,k)-1)/2
               !! if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
             !! if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then

!!!!!u
                bl_intp_valx=block(a_blk_no)%xu(tar_x)
                bl_intp_valy=block(a_blk_no)%yu(tar_y)
                bl_intp_valz=block(a_blk_no)%zu(tar_z)
                bl_intp_x1=block(b_blk_no)%xu(loc_x-1)
                bl_intp_x2=block(b_blk_no)%xu(loc_x)
                bl_intp_y1=block(b_blk_no)%yu(loc_y)
                bl_intp_y2=block(b_blk_no)%yu(loc_y+1)
                bl_intp_z1=block(b_blk_no)%zu(loc_z)
                bl_intp_z2=block(b_blk_no)%zu(loc_z+1)
                bl_intp_f1=block(b_blk_no)%u(loc_x-2,loc_y,loc_z)
                bl_intp_f2=block(b_blk_no)%u(loc_x-1,loc_y,loc_z)
                bl_intp_f3=block(b_blk_no)%u(loc_x-1,loc_y+1,loc_z)
                bl_intp_f4=block(b_blk_no)%u(loc_x-2,loc_y+1,loc_z)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)

                bl_intp_f1=block(b_blk_no)%u(loc_x-2,loc_y,loc_z+1)
                bl_intp_f2=block(b_blk_no)%u(loc_x-1,loc_y,loc_z+1)
                bl_intp_f3=block(b_blk_no)%u(loc_x-1,loc_y+1,loc_z+1)
                bl_intp_f4=block(b_blk_no)%u(loc_x-2,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans

               !! endif

        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxv-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryv-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzv-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
     ! do k=2,intfr(g)%counterzv-1
     ! do j=2,intfr(g)%counteryv-1
     ! do i=2,intfr(g)%counterxv-1
      !do i=2, intfr(g)%counterxv-1
      ! do j=2,intfr(g)%counteryv-1
                tar_x=intfr(g)%vx_interface_det(1,i)
                tar_y=intfr(g)%vy_interface_det(1,j)
                tar_z=intfr(g)%vz_interface_det(1,k)
                loc_x=(intfr(g)%vx_interface_det(2,i)+intfr(g)%vx_interface_det(3,i)-1)/2
                loc_z=(intfr(g)%vz_interface_det(2,k)+intfr(g)%vz_interface_det(3,k)-1)/2
                !loc_y=(intfr(g)%vy_interface_det(2,j)+intfr(g)%vy_interface_det(3,j)-1)/2
                loc_y=(intfr(g)%vy_interface_det(3,j))
                !!if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
              !!if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then

!!!!!v
                bl_intp_valx=block(a_blk_no)%xv(tar_x)
                bl_intp_valy=block(a_blk_no)%yv(tar_y)
                bl_intp_valz=block(a_blk_no)%zv(tar_z)
                bl_intp_x1=block(b_blk_no)%xv(loc_x)
                bl_intp_x2=block(b_blk_no)%xv(loc_x+1)
                bl_intp_y1=block(b_blk_no)%yv(loc_y-1)
                bl_intp_y2=block(b_blk_no)%yv(loc_y)
                bl_intp_z1=block(b_blk_no)%zv(loc_z)
                bl_intp_z2=block(b_blk_no)%zv(loc_z+1)
                bl_intp_f1=block(b_blk_no)%v(loc_x,loc_y-2,loc_z)
                bl_intp_f2=block(b_blk_no)%v(loc_x+1,loc_y-2,loc_z)
                bl_intp_f3=block(b_blk_no)%v(loc_x+1,loc_y-1,loc_z)
                bl_intp_f4=block(b_blk_no)%v(loc_x,loc_y-1,loc_z)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(b_blk_no)%v(loc_x,loc_y-2,loc_z+1)
                bl_intp_f2=block(b_blk_no)%v(loc_x+1,loc_y-2,loc_z+1)
                bl_intp_f3=block(b_blk_no)%v(loc_x+1,loc_y-1,loc_z+1)
                bl_intp_f4=block(b_blk_no)%v(loc_x,loc_y-1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans
              !! endif


        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxw-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryw-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzw-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
      !do k=2,intfr(g)%counterzw-1
      !do j=2,intfr(g)%counteryw-1
      !do i=2,intfr(g)%counterxw-1
      !do i=2, intfr(g)%counterxv-1
      ! do j=2,intfr(g)%counteryv-1
                tar_x=intfr(g)%wx_interface_det(1,i)
                tar_y=intfr(g)%wy_interface_det(1,j)
                tar_z=intfr(g)%wz_interface_det(1,k)
                loc_x=(intfr(g)%wx_interface_det(2,i)+intfr(g)%wx_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%wy_interface_det(2,j)+intfr(g)%wy_interface_det(3,j)-1)/2
                !loc_y=(intfr(g)%vy_interface_det(2,j)+intfr(g)%vy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%wz_interface_det(3,k))
                !!if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
              !!if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then

!!!!!v
                bl_intp_valx=block(a_blk_no)%xw(tar_x)
                bl_intp_valy=block(a_blk_no)%yw(tar_y)
                bl_intp_valz=block(a_blk_no)%zw(tar_z)
                bl_intp_x1=block(b_blk_no)%xw(loc_x)
                bl_intp_x2=block(b_blk_no)%xw(loc_x+1)
                bl_intp_y1=block(b_blk_no)%yw(loc_y)
                bl_intp_y2=block(b_blk_no)%yw(loc_y+1)
                bl_intp_z1=block(b_blk_no)%zw(loc_z-1)
                bl_intp_z2=block(b_blk_no)%zw(loc_z)
                bl_intp_f1=block(b_blk_no)%w(loc_x,loc_y,loc_z-2)
                bl_intp_f2=block(b_blk_no)%w(loc_x+1,loc_y,loc_z-2)
                bl_intp_f3=block(b_blk_no)%w(loc_x+1,loc_y+1,loc_z-2)
                bl_intp_f4=block(b_blk_no)%w(loc_x,loc_y+1,loc_z-2)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(b_blk_no)%w(loc_x,loc_y,loc_z-1)
                bl_intp_f2=block(b_blk_no)%w(loc_x+1,loc_y,loc_z-1)
                bl_intp_f3=block(b_blk_no)%w(loc_x+1,loc_y+1,loc_z-1)
                bl_intp_f4=block(b_blk_no)%w(loc_x,loc_y+1,loc_z-1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans
              !! endif


        enddo
        enddo
        enddo
        ENDDO
        END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!      SUBROUTINE coarsept_Update
!!      use global

!!      INTEGER :: i,j,k,g,h, varx1,varx2, tar_x, tar_y, loc_x, loc_y, a_blk_no,b_blk_no
!!      CHARACTER*150 :: filename1




!!      DO g=1,intflines
!!      a_blk_no=intfr(g)%a_blk
!!      b_blk_no=intfr(g)%b_blk
!!      factor=b_blk_no/a_blk_no



!       i=intfr(g)%fpx_start
!       j=intfr(g)%fpy_start
!       DO k=intfr(g)%cpy_start,intfr(g)%cpy_end
!       DO h=intfr(g)%cpx_start,intfr(g)%cpx_end
!
!       block(a_blk_no)%p(h,k)=block(b_blk_no)%p(i,j)
!       !print*,block(a_blk_no)%xp(h),block(b_blk_no)%xp(i)
!       i=i+factor
!       ENDDO
!       j=j+factor
!       ENDDO

        !OPEN(UNIT = 786, FILE = 'check.dat', ACCESS= 'append', STATUS = 'unknown')
!!      j=intfr(g)%fcelly_start
!!      DO k=intfr(g)%ccelly_start,intfr(g)%ccelly_end
!!      i=intfr(g)%fcellx_start
!!      DO h=intfr(g)%ccellx_start,intfr(g)%ccellx_end
!!
!!      !write(786,*)h,i,block(a_blk_no)%xu(h)-block(b_blk_no)%xu(i)
!!      !write(786,*)block(a_blk_no)%xu(h),block(b_blk_no)%xu(i)
!!      !!write(786,*)block(a_blk_no)%yv(k)-block(b_blk_no)%yv(j)
!!      !write(786,*)'--------------------------------------'
!!      block(a_blk_no)%cell(h,k)=block(b_blk_no)%cell(i,j)
!!      i=i+factor
!!      ENDDO
!!      j=j+factor
!!      ENDDO
!!      !CLOSE(786)
!!      if ( g .eq. intflines)then
!!      st_rc_x=block(a_blk_no)%irc_st
!!      en_rc_x=block(a_blk_no)%irc_en
!!      st_rc_y=block(a_blk_no)%jrc_st
!!      en_rc_y=block(a_blk_no)%jrc_en

!!      block(a_blk_no)%ibCellCount = 0
!!      block(a_blk_no)%fluidCellCount = 0
!!      block(a_blk_no)%solidCellCount = 0
!!      DO 20 j = st_rc_y, en_rc_y
!!      DO 20 i = st_rc_x, en_rc_x
!!         ! n   = i-1  + (block(g)%nx)*(j-2)
!!          IF (block(a_blk_no)%cell(i,j).eq.1) THEN
!!      	         block(a_blk_no)%solidCellCount = block(a_blk_no)%solidCellCount + 1
!!          ELSEIF (block(a_blk_no)%cell(i,j).eq.0) THEN
!!             block(a_blk_no)%fluidCellCount  = block(a_blk_no)%fluidCellCount + 1
!!          ELSEIF (block(a_blk_no)%cell(i,j).eq.2) THEN
!!             block(a_blk_no)%ibCellCount = block(a_blk_no)%ibCellCount + 1
!!          ENDIF
!! 20     CONTINUE
!!        endif
!       i=intfr(g)%fux_start
!       j=intfr(g)%fuy_start
!       DO k=intfr(g)%cuy_start,intfr(g)%cuy_end
!       DO h=intfr(g)%cux_start,intfr(g)%cux_end
!
!       block(a_blk_no)%u(h,k)=block(b_blk_no)%u(i,j)
!       i=i+factor
!       ENDDO
!       j=j+factor
!       ENDDO


!       i=intfr(g)%fvx_start
!       j=intfr(g)%fvy_start
!       DO k=intfr(g)%cvy_start,intfr(g)%cvy_end
!       DO h=intfr(g)%cvx_start,intfr(g)%cvx_end
!
!       block(a_blk_no)%v(h,k)=block(b_blk_no)%v(i,j)
!       i=i+factor
!       ENDDO
!       j=j+factor
!       ENDDO



!!        ENDDO
!!        END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       subroutine coarseUpdate_newv
!       use global
!       REAL (KIND=8) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
!       REAL (KIND=8) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
!       REAL (KIND=8) :: bl_interp_ans

!       INTEGER :: i,j,k,g, varx1,varx2, tar_x, tar_y, loc_x, loc_y, a_blk_no,b_blk_no
!       INTEGER :: st_idx, en_idx
!       INTEGER :: st_idy, en_idy
!      !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
!      !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
!      !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
!      !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
!      !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
!      !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!      !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!      !bl_intp_num= bl_intp_first_term + bl_intp_second_term
!      !
!      !bl_interp_ans=(bl_intp_num/bl_intp_deno)

!

!       DO g=1,intflines
!       a_blk_no=intfr(g)%a_blk
!       b_blk_no=intfr(g)%b_blk
!       !counterx= (intfr(g)%xintf_end-intfr(g)%xintf_start +1)/factor
!       !countery= (intfr(g)%yintf_end-intfr(g)%yintf_start +1)/factor

!      do i=2, intfr(g)%counterxp-1
!       do j=2,intfr(g)%counteryp-1
!               tar_x=intfr(g)%px_interface_det(1,i)
!               tar_y=intfr(g)%py_interface_det(1,j)
!               loc_x=(intfr(g)%px_interface_det(2,i)+intfr(g)%px_interface_det(3,i)-1)/2
!               loc_y=(intfr(g)%py_interface_det(2,j)+intfr(g)%py_interface_det(3,j)-1)/2
!             !!  if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
!
!!!!!p
!               bl_intp_valx=block(a_blk_no)%xp(tar_x)
!               bl_intp_valy=block(a_blk_no)%yp(tar_y)
!               bl_intp_x1=block(b_blk_no)%xp(loc_x)
!               bl_intp_x2=block(b_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(b_blk_no)%yp(loc_y)
!               bl_intp_y2=block(b_blk_no)%yp(loc_y+1)
!               bl_intp_f1=block(b_blk_no)%p(loc_x,loc_y)
!               bl_intp_f2=block(b_blk_no)%p(loc_x+1,loc_y)
!               bl_intp_f3=block(b_blk_no)%p(loc_x+1,loc_y+1)
!               bl_intp_f4=block(b_blk_no)%p(loc_x,loc_y+1)
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

!               block(a_blk_no)%p(tar_x,tar_y)=bl_interp_ans
!!!!!u
!               !!endif
!

!       enddo
!       enddo
!
!       st_idx=block(b_blk_no)%cintp
!       en_idx=intfr(g)%counterxu-block(b_blk_no)%cintp
!       st_idy=block(b_blk_no)%cintp
!       en_idy=intfr(g)%counteryu-block(b_blk_no)%cintp
!      do i=st_idx, en_idx
!       do j=st_idy,en_idy
!     !do i=2, intfr(g)%counterxu-1
!     ! do j=2,intfr(g)%counteryu-1
!               tar_x=intfr(g)%ux_interface_det(1,i)
!               tar_y=intfr(g)%uy_interface_det(1,j)
!               !loc_x=(intfr(g)%ux_interface_det(2,i)+intfr(g)%ux_interface_det(3,i)-1)/2
!               loc_x=intfr(g)%ux_interface_det(3,i)
!               loc_y=(intfr(g)%uy_interface_det(2,j)+intfr(g)%uy_interface_det(3,j)-1)/2
!             !! if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
!
!!!!!u
!               bl_intp_valx=block(a_blk_no)%xu(tar_x)
!               bl_intp_valy=block(a_blk_no)%yu(tar_y)
!               bl_intp_x1=block(b_blk_no)%xu(loc_x-1)
!               bl_intp_x2=block(b_blk_no)%xu(loc_x)
!               bl_intp_y1=block(b_blk_no)%yu(loc_y)
!               bl_intp_y2=block(b_blk_no)%yu(loc_y+1)
!               bl_intp_f1=block(b_blk_no)%ut(loc_x-2,loc_y)
!               bl_intp_f2=block(b_blk_no)%ut(loc_x-1,loc_y)
!               bl_intp_f3=block(b_blk_no)%ut(loc_x-1,loc_y+1)
!               bl_intp_f4=block(b_blk_no)%ut(loc_x-2,loc_y+1)
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

!               block(a_blk_no)%ut(tar_x-1,tar_y)=bl_interp_ans
!
!              !! endif

!       enddo
!       enddo
!
!       st_idx=block(b_blk_no)%cintp
!       en_idx=intfr(g)%counterxv-block(b_blk_no)%cintp
!       st_idy=block(b_blk_no)%cintp
!       en_idy=intfr(g)%counteryv-block(b_blk_no)%cintp
!      do i=st_idx, en_idx
!       do j=st_idy,en_idy
!     !do i=2, intfr(g)%counterxv-1
!     ! do j=2,intfr(g)%counteryv-1
!               tar_x=intfr(g)%vx_interface_det(1,i)
!               tar_y=intfr(g)%vy_interface_det(1,j)
!               loc_x=(intfr(g)%vx_interface_det(2,i)+intfr(g)%vx_interface_det(3,i)-1)/2
!               !loc_y=(intfr(g)%vy_interface_det(2,j)+intfr(g)%vy_interface_det(3,j)-1)/2
!               loc_y=(intfr(g)%vy_interface_det(3,j)-1)
!              !! if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then

!!!!!v
!               bl_intp_valx=block(a_blk_no)%xv(tar_x)
!               bl_intp_valy=block(a_blk_no)%yv(tar_y)
!               bl_intp_x1=block(b_blk_no)%xv(loc_x)
!               bl_intp_x2=block(b_blk_no)%xv(loc_x+1)
!               bl_intp_y1=block(b_blk_no)%yv(loc_y-1)
!               bl_intp_y2=block(b_blk_no)%yv(loc_y)
!               bl_intp_f1=block(b_blk_no)%vt(loc_x,loc_y-2)
!               bl_intp_f2=block(b_blk_no)%vt(loc_x+1,loc_y-2)
!               bl_intp_f3=block(b_blk_no)%vt(loc_x+1,loc_y-1)
!               bl_intp_f4=block(b_blk_no)%vt(loc_x,loc_y-1)
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

!               block(a_blk_no)%vt(tar_x,tar_y-1)=bl_interp_ans
!              !! endif
!

!       enddo
!       enddo
!       ENDDO
!       END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
        subroutine coarseUpdate_newv
        use global
        IMPLICIT NONE
        REAL (KIND=8) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (KIND=8) :: bl_intp_valz,bl_intp_z1,bl_intp_z2,bl_intp_f5,bl_intp_f6,bl_intp_f7,bl_intp_f8
        REAL (KIND=8) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (KIND=8) :: bl_interp_ans1, bl_interp_ans2, bl_interp_ansf
        REAL (KIND=8) :: bl_interp_ans

        INTEGER :: i,j,k,g, varx1,varx2, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz
        CHARACTER(len=150) filename1

       !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
       !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
       !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
       !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
       !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
       !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
       !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
       !bl_intp_num= bl_intp_first_term + bl_intp_second_term
       !
       !bl_interp_ans=(bl_intp_num/bl_intp_deno)



        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        !counterx= (intfr(g)%xintf_end-intfr(g)%xintf_start +1)/factor
        !countery= (intfr(g)%yintf_end-intfr(g)%yintf_start +1)/factor
!       st_idx=block(b_blk_no)%cintp
!       en_idx=intfr(g)%counterxp-block(b_blk_no)%cintp
!       st_idy=block(b_blk_no)%cintp
!       en_idy=intfr(g)%counteryp-block(b_blk_no)%cintp
!       st_idz=block(b_blk_no)%cintp
!       en_idz=intfr(g)%counterzp-block(b_blk_no)%cintp
!     !do i=2, intfr(g)%counterxp-1
!     ! do j=2,intfr(g)%counteryp-1
!      do k=st_idz,en_idz
!      do j=st_idy,en_idy
!      do i=st_idx, en_idx
!               tar_x=intfr(g)%px_interface_det(1,i)
!               tar_y=intfr(g)%py_interface_det(1,j)
!               tar_z=intfr(g)%pz_interface_det(1,k)
!               loc_x=(intfr(g)%px_interface_det(2,i)+intfr(g)%px_interface_det(3,i)-1)/2
!               loc_y=(intfr(g)%py_interface_det(2,j)+intfr(g)%py_interface_det(3,j)-1)/2
!               loc_z=(intfr(g)%pz_interface_det(2,k)+intfr(g)%pz_interface_det(3,k)-1)/2
!             !!  if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
!             !!if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
!!!!!p
!               bl_intp_valx=block(a_blk_no)%xp(tar_x)
!               bl_intp_valy=block(a_blk_no)%yp(tar_y)
!               bl_intp_valz=block(a_blk_no)%zp(tar_z)
!               bl_intp_x1=block(b_blk_no)%xp(loc_x)
!               bl_intp_x2=block(b_blk_no)%xp(loc_x+1)
!               bl_intp_y1=block(b_blk_no)%yp(loc_y)
!               bl_intp_y2=block(b_blk_no)%yp(loc_y+1)
!               bl_intp_z1=block(b_blk_no)%zp(loc_z)
!               bl_intp_z2=block(b_blk_no)%zp(loc_z+1)
!               bl_intp_f1=block(b_blk_no)%p(loc_x,loc_y,loc_z)
!               bl_intp_f2=block(b_blk_no)%p(loc_x+1,loc_y,loc_z)
!               bl_intp_f3=block(b_blk_no)%p(loc_x+1,loc_y+1,loc_z)
!               bl_intp_f4=block(b_blk_no)%p(loc_x,loc_y+1,loc_z)
!              !bl_intp_f5=block(b_blk_no)%p(loc_x,loc_y,loc_z+1)
!              !bl_intp_f6=block(b_blk_no)%p(loc_x+1,loc_y,loc_z+1)
!              !bl_intp_f7=block(b_blk_no)%p(loc_x+1,loc_y+1,loc_z+1)
!              !bl_intp_f8=block(b_blk_no)%p(loc_x,loc_y+1,loc_z+1)
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
!               bl_interp_ans1=(bl_intp_num/bl_intp_deno)
!
!               bl_intp_f1=block(b_blk_no)%p(loc_x,loc_y,loc_z+1)
!               bl_intp_f2=block(b_blk_no)%p(loc_x+1,loc_y,loc_z+1)
!               bl_intp_f3=block(b_blk_no)%p(loc_x+1,loc_y+1,loc_z+1)
!               bl_intp_f4=block(b_blk_no)%p(loc_x,loc_y+1,loc_z+1)
!               bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
!               bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
!               bl_intp_num= bl_intp_first_term + bl_intp_second_term

!               bl_interp_ans2=(bl_intp_num/bl_intp_deno)
!               bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
!
!               block(a_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans
!!!!!u
!              !! endif
!

!       enddo
!       enddo
!
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxu-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryu-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzu-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
      !do k=2,intfr(g)%counterzu-1
      !do j=2,intfr(g)%counteryu-1
      !do i=2,intfr(g)%counterxu-1
      !do i=2, intfr(g)%counterxu-1
      ! do j=2,intfr(g)%counteryu-1
                tar_x=intfr(g)%ux_interface_det(1,i)
                tar_y=intfr(g)%uy_interface_det(1,j)
                tar_z=intfr(g)%uz_interface_det(1,k)
                !loc_x=(intfr(g)%ux_interface_det(2,i)+intfr(g)%ux_interface_det(3,i)-1)/2
                loc_x=intfr(g)%ux_interface_det(3,i)
                loc_y=(intfr(g)%uy_interface_det(2,j)+intfr(g)%uy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%uz_interface_det(2,k)+intfr(g)%uz_interface_det(3,k)-1)/2
               !! if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
             !! if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then

!!!!!u
                bl_intp_valx=block(a_blk_no)%xu(tar_x)
                bl_intp_valy=block(a_blk_no)%yu(tar_y)
                bl_intp_valz=block(a_blk_no)%zu(tar_z)
                bl_intp_x1=block(b_blk_no)%xu(loc_x-1)
                bl_intp_x2=block(b_blk_no)%xu(loc_x)
                bl_intp_y1=block(b_blk_no)%yu(loc_y)
                bl_intp_y2=block(b_blk_no)%yu(loc_y+1)
                bl_intp_z1=block(b_blk_no)%zu(loc_z)
                bl_intp_z2=block(b_blk_no)%zu(loc_z+1)
                bl_intp_f1=block(b_blk_no)%ut(loc_x-2,loc_y,loc_z)
                bl_intp_f2=block(b_blk_no)%ut(loc_x-1,loc_y,loc_z)
                bl_intp_f3=block(b_blk_no)%ut(loc_x-1,loc_y+1,loc_z)
                bl_intp_f4=block(b_blk_no)%ut(loc_x-2,loc_y+1,loc_z)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)

                bl_intp_f1=block(b_blk_no)%ut(loc_x-2,loc_y,loc_z+1)
                bl_intp_f2=block(b_blk_no)%ut(loc_x-1,loc_y,loc_z+1)
                bl_intp_f3=block(b_blk_no)%ut(loc_x-1,loc_y+1,loc_z+1)
                bl_intp_f4=block(b_blk_no)%ut(loc_x-2,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%ut(tar_x-1,tar_y,tar_z)=bl_interp_ans

!               write(filename1,1917)
! 1917            format('ucinterp_print.dat')
!               OPEN(UNIT=1221,FILE=filename1,ACCESS='append',STATUS='unknown')
!              ! write(1221,391) ita,' ',loc_x,' ',loc_y,' ',loc_z,' ',tar_x,' ',tar_y,' ',tar_z,' ',bl_intp_x1, bl_intp_x2,bl_intp_y1, bl_intp_y2,bl_intp_z1, bl_intp_z2, bl_intp_f1, bl_intp_f2, bl_intp_f3, bl_intp_f4
!               write(1221,391) ita,' ',loc_x,' ',loc_y,' ',loc_z,' ',tar_x,' ',tar_y,' ',tar_z,' ',bl_intp_x1,bl_intp_x2,bl_intp_y1, bl_intp_y2,bl_intp_z1, bl_intp_z2, bl_intp_valx, bl_intp_valy, bl_intp_valz
! !391             format(I4,A1,I4,A1,I4,A1,I4,A1,I4,A1,I4,A1,I4,A1,10F10.7)
! 391             format(I4,A1,I4,A1,I4,A1,I4,A1,I4,A1,I4,A1,I4,A1,9F10.7)
               !! endif

        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxv-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryv-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzv-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
      !do k=2,intfr(g)%counterzv-1
      !do j=2,intfr(g)%counteryv-1
      !do i=2,intfr(g)%counterxv-1
      !do i=2, intfr(g)%counterxv-1
      ! do j=2,intfr(g)%counteryv-1
                tar_x=intfr(g)%vx_interface_det(1,i)
                tar_y=intfr(g)%vy_interface_det(1,j)
                tar_z=intfr(g)%vz_interface_det(1,k)
                loc_x=(intfr(g)%vx_interface_det(2,i)+intfr(g)%vx_interface_det(3,i)-1)/2
                loc_z=(intfr(g)%vz_interface_det(2,k)+intfr(g)%vz_interface_det(3,k)-1)/2
                !loc_y=(intfr(g)%vy_interface_det(2,j)+intfr(g)%vy_interface_det(3,j)-1)/2
                loc_y=(intfr(g)%vy_interface_det(3,j))
                !!if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
              !!if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then

!!!!!v
                bl_intp_valx=block(a_blk_no)%xv(tar_x)
                bl_intp_valy=block(a_blk_no)%yv(tar_y)
                bl_intp_valz=block(a_blk_no)%zv(tar_z)
                bl_intp_x1=block(b_blk_no)%xv(loc_x)
                bl_intp_x2=block(b_blk_no)%xv(loc_x+1)
                bl_intp_y1=block(b_blk_no)%yv(loc_y-1)
                bl_intp_y2=block(b_blk_no)%yv(loc_y)
                bl_intp_z1=block(b_blk_no)%zv(loc_z)
                bl_intp_z2=block(b_blk_no)%zv(loc_z+1)
                bl_intp_f1=block(b_blk_no)%vt(loc_x,loc_y-2,loc_z)
                bl_intp_f2=block(b_blk_no)%vt(loc_x+1,loc_y-2,loc_z)
                bl_intp_f3=block(b_blk_no)%vt(loc_x+1,loc_y-1,loc_z)
                bl_intp_f4=block(b_blk_no)%vt(loc_x,loc_y-1,loc_z)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(b_blk_no)%vt(loc_x,loc_y-2,loc_z+1)
                bl_intp_f2=block(b_blk_no)%vt(loc_x+1,loc_y-2,loc_z+1)
                bl_intp_f3=block(b_blk_no)%vt(loc_x+1,loc_y-1,loc_z+1)
                bl_intp_f4=block(b_blk_no)%vt(loc_x,loc_y-1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%vt(tar_x,tar_y-1,tar_z)=bl_interp_ans
              !! endif


        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxw-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryw-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzw-block(b_blk_no)%cintp
      !do k=2,intfr(g)%counterzw-1
      !do j=2,intfr(g)%counteryw-1
      !do i=2,intfr(g)%counterxw-1
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
      !do i=2, intfr(g)%counterxv-1
      ! do j=2,intfr(g)%counteryv-1
                tar_x=intfr(g)%wx_interface_det(1,i)
                tar_y=intfr(g)%wy_interface_det(1,j)
                tar_z=intfr(g)%wz_interface_det(1,k)
                loc_x=(intfr(g)%wx_interface_det(2,i)+intfr(g)%wx_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%wy_interface_det(2,j)+intfr(g)%wy_interface_det(3,j)-1)/2
                !loc_y=(intfr(g)%vy_interface_det(2,j)+intfr(g)%vy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%wz_interface_det(3,k))
                !!if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
              !!if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then

!!!!!v
                bl_intp_valx=block(a_blk_no)%xw(tar_x)
                bl_intp_valy=block(a_blk_no)%yw(tar_y)
                bl_intp_valz=block(a_blk_no)%zw(tar_z)
                bl_intp_x1=block(b_blk_no)%xw(loc_x)
                bl_intp_x2=block(b_blk_no)%xw(loc_x+1)
                bl_intp_y1=block(b_blk_no)%yw(loc_y)
                bl_intp_y2=block(b_blk_no)%yw(loc_y+1)
                bl_intp_z1=block(b_blk_no)%zw(loc_z-1)
                bl_intp_z2=block(b_blk_no)%zw(loc_z)
                bl_intp_f1=block(b_blk_no)%wt(loc_x,loc_y,loc_z-2)
                bl_intp_f2=block(b_blk_no)%wt(loc_x+1,loc_y,loc_z-2)
                bl_intp_f3=block(b_blk_no)%wt(loc_x+1,loc_y+1,loc_z-2)
                bl_intp_f4=block(b_blk_no)%wt(loc_x,loc_y+1,loc_z-2)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(b_blk_no)%w(loc_x,loc_y,loc_z-1)
                bl_intp_f2=block(b_blk_no)%w(loc_x+1,loc_y,loc_z-1)
                bl_intp_f3=block(b_blk_no)%w(loc_x+1,loc_y+1,loc_z-1)
                bl_intp_f4=block(b_blk_no)%w(loc_x,loc_y+1,loc_z-1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%wt(tar_x,tar_y,tar_z-1)=bl_interp_ans
              !! endif


        enddo
        enddo
        enddo
        ENDDO
        END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
        subroutine coarseUpdate_pc
        use global
        IMPLICIT NONE
        REAL (KIND=8) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (KIND=8) :: bl_intp_valz,bl_intp_z1,bl_intp_z2,bl_intp_f5,bl_intp_f6,bl_intp_f7,bl_intp_f8
        REAL (KIND=8) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (KIND=8) :: bl_interp_ans1, bl_interp_ans2, bl_interp_ansf
        REAL (KIND=8) :: bl_interp_ans

        INTEGER :: i,j,k,g, varx1,varx2, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz
       !bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
       !bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
       !bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
       !bl_intp_yty=(bl_intp_y2-bl_intp_valy)
       !bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
       !bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
       !bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
       !bl_intp_num= bl_intp_first_term + bl_intp_second_term
       !
       !bl_interp_ans=(bl_intp_num/bl_intp_deno)



        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        !counterx= (intfr(g)%xintf_end-intfr(g)%xintf_start +1)/factor
        !countery= (intfr(g)%yintf_end-intfr(g)%yintf_start +1)/factor
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxp-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryp-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzp-block(b_blk_no)%cintp
      !do i=2, intfr(g)%counterxp-1
      ! do j=2,intfr(g)%counteryp-1
     ! do k=2,intfr(g)%counterzp-1
     ! do j=2,intfr(g)%counteryp-1
     ! do i=2,intfr(g)%counterxp-1
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%px_interface_det(1,i)
                tar_y=intfr(g)%py_interface_det(1,j)
                tar_z=intfr(g)%pz_interface_det(1,k)
                loc_x=(intfr(g)%px_interface_det(2,i)+intfr(g)%px_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%py_interface_det(2,j)+intfr(g)%py_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%pz_interface_det(2,k)+intfr(g)%pz_interface_det(3,k)-1)/2
              !!  if( block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
              !!if (block(a_blk_no)%cell(tar_x,tar_y) .eq. 1)then
 !!!!p
                bl_intp_valx=block(a_blk_no)%xp(tar_x)
                bl_intp_valy=block(a_blk_no)%yp(tar_y)
                bl_intp_valz=block(a_blk_no)%zp(tar_z)
                bl_intp_x1=block(b_blk_no)%xp(loc_x)
                bl_intp_x2=block(b_blk_no)%xp(loc_x+1)
                bl_intp_y1=block(b_blk_no)%yp(loc_y)
                bl_intp_y2=block(b_blk_no)%yp(loc_y+1)
                bl_intp_z1=block(b_blk_no)%zp(loc_z)
                bl_intp_z2=block(b_blk_no)%zp(loc_z+1)
                bl_intp_f1=block(b_blk_no)%pc(loc_x,loc_y,loc_z)
                bl_intp_f2=block(b_blk_no)%pc(loc_x+1,loc_y,loc_z)
                bl_intp_f3=block(b_blk_no)%pc(loc_x+1,loc_y+1,loc_z)
                bl_intp_f4=block(b_blk_no)%pc(loc_x,loc_y+1,loc_z)
               !bl_intp_f5=block(b_blk_no)%p(loc_x,loc_y,loc_z+1)
               !bl_intp_f6=block(b_blk_no)%p(loc_x+1,loc_y,loc_z+1)
               !bl_intp_f7=block(b_blk_no)%p(loc_x+1,loc_y+1,loc_z+1)
               !bl_intp_f8=block(b_blk_no)%p(loc_x,loc_y+1,loc_z+1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)

                bl_intp_f1=block(b_blk_no)%pc(loc_x,loc_y,loc_z+1)
                bl_intp_f2=block(b_blk_no)%pc(loc_x+1,loc_y,loc_z+1)
                bl_intp_f3=block(b_blk_no)%pc(loc_x+1,loc_y+1,loc_z+1)
                bl_intp_f4=block(b_blk_no)%pc(loc_x,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
             !  write(*,121)'c',i,j,k,bl_interp_ans1, bl_interp_ans2, bl_interp_ans
!121          format(A1,3I,3F12.6)
                block(a_blk_no)%pc(tar_x,tar_y,tar_z)=bl_interp_ans
                block(a_blk_no)%pco(tar_x,tar_y,tar_z)=bl_interp_ans
 !!!!u
               !! endif


        enddo
        enddo
        enddo
        ENDDO
        END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

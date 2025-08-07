module biocfd_coarse_update
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use global, only : block, intfr, intflines
  implicit none
  private

  public ::  coarseUpdate_newv, coarseUpdate_pc, coarseUpdate

  contains
subroutine coarseUpdate
        REAL (dp) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,&
                         bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (dp) :: bl_intp_valz,bl_intp_z1,bl_intp_z2
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, &
                         bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp) :: bl_interp_ans1, bl_interp_ans2
        REAL (dp) :: bl_interp_ans

        INTEGER :: i,j,k,g, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz

        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxp-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryp-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzp-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%px_interface_det(1,i)
                tar_y=intfr(g)%py_interface_det(1,j)
                tar_z=intfr(g)%pz_interface_det(1,k)
                loc_x=(intfr(g)%px_interface_det(2,i)+intfr(g)%px_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%py_interface_det(2,j)+intfr(g)%py_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%pz_interface_det(2,k)+intfr(g)%pz_interface_det(3,k)-1)/2
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))

                block(a_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans


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
                tar_x=intfr(g)%ux_interface_det(1,i)
                tar_y=intfr(g)%uy_interface_det(1,j)
                tar_z=intfr(g)%uz_interface_det(1,k)
                loc_x=intfr(g)%ux_interface_det(3,i)
                loc_y=(intfr(g)%uy_interface_det(2,j)+intfr(g)%uy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%uz_interface_det(2,k)+intfr(g)%uz_interface_det(3,k)-1)/2
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans


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
                tar_x=intfr(g)%vx_interface_det(1,i)
                tar_y=intfr(g)%vy_interface_det(1,j)
                tar_z=intfr(g)%vz_interface_det(1,k)
                loc_x=(intfr(g)%vx_interface_det(2,i)+intfr(g)%vx_interface_det(3,i)-1)/2
                loc_z=(intfr(g)%vz_interface_det(2,k)+intfr(g)%vz_interface_det(3,k)-1)/2
                loc_y=(intfr(g)%vy_interface_det(3,j))
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans


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
                tar_x=intfr(g)%wx_interface_det(1,i)
                tar_y=intfr(g)%wy_interface_det(1,j)
                tar_z=intfr(g)%wz_interface_det(1,k)
                loc_x=(intfr(g)%wx_interface_det(2,i)+intfr(g)%wx_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%wy_interface_det(2,j)+intfr(g)%wy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%wz_interface_det(3,k))
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans


        enddo
        enddo
        enddo
        ENDDO
        end subroutine coarseUpdate

        subroutine coarseUpdate_newv
        REAL (dp) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,&
               bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (dp) :: bl_intp_valz,bl_intp_z1,bl_intp_z2
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty,&
               bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp) :: bl_interp_ans1, bl_interp_ans2
        REAL (dp) :: bl_interp_ans

        INTEGER :: i,j,k,g, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz

        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxu-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryu-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzu-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%ux_interface_det(1,i)
                tar_y=intfr(g)%uy_interface_det(1,j)
                tar_z=intfr(g)%uz_interface_det(1,k)
                loc_x=intfr(g)%ux_interface_det(3,i)
                loc_y=(intfr(g)%uy_interface_det(2,j)+intfr(g)%uy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%uz_interface_det(2,k)+intfr(g)%uz_interface_det(3,k)-1)/2
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%ut(tar_x-1,tar_y,tar_z)=bl_interp_ans
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
                tar_x=intfr(g)%vx_interface_det(1,i)
                tar_y=intfr(g)%vy_interface_det(1,j)
                tar_z=intfr(g)%vz_interface_det(1,k)
                loc_x=(intfr(g)%vx_interface_det(2,i)+intfr(g)%vx_interface_det(3,i)-1)/2
                loc_z=(intfr(g)%vz_interface_det(2,k)+intfr(g)%vz_interface_det(3,k)-1)/2
                loc_y=(intfr(g)%vy_interface_det(3,j))
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%vt(tar_x,tar_y-1,tar_z)=bl_interp_ans

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
                tar_x=intfr(g)%wx_interface_det(1,i)
                tar_y=intfr(g)%wy_interface_det(1,j)
                tar_z=intfr(g)%wz_interface_det(1,k)
                loc_x=(intfr(g)%wx_interface_det(2,i)+intfr(g)%wx_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%wy_interface_det(2,j)+intfr(g)%wy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%wz_interface_det(3,k))
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%wt(tar_x,tar_y,tar_z-1)=bl_interp_ans


        enddo
        enddo
        enddo
        ENDDO
        end subroutine coarseUpdate_newv
        subroutine coarseUpdate_pc
        REAL (dp) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,&
               bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (dp) :: bl_intp_valz,bl_intp_z1,bl_intp_z2
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, &
             bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp) :: bl_interp_ans1, bl_interp_ans2
        REAL (dp) :: bl_interp_ans

        INTEGER :: i,j,k,g, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz


        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxp-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryp-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzp-block(b_blk_no)%cintp
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%px_interface_det(1,i)
                tar_y=intfr(g)%py_interface_det(1,j)
                tar_z=intfr(g)%pz_interface_det(1,k)
                loc_x=(intfr(g)%px_interface_det(2,i)+intfr(g)%px_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%py_interface_det(2,j)+intfr(g)%py_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%pz_interface_det(2,k)+intfr(g)%pz_interface_det(3,k)-1)/2
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(a_blk_no)%pc(tar_x,tar_y,tar_z)=bl_interp_ans
                block(a_blk_no)%pco(tar_x,tar_y,tar_z)=bl_interp_ans


        enddo
        enddo
        enddo
        ENDDO
        end subroutine coarseUpdate_pc
end module biocfd_coarse_update

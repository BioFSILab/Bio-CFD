module biocfd_fine_interp
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, intfr
  implicit none
  private

  public :: fineUpdate_mv
  contains

        SUBROUTINE fineUpdate_mv(g)
        REAL (dp) :: bl_intp_valx,bl_intp_valy,bl_intp_valz,bl_intp_x1,bl_intp_x2,&
             bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (dp) :: bl_intp_z1,bl_intp_z2
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty,&
             bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp) :: bl_interp_ans, bl_interp_ans1, bl_interp_ans2

        INTEGER(dp) :: i,j,k,q,s,varx1,varx2,vary1,vary2,l,tar_x,tar_y,loc_x,loc_y,&
             a_blk_no,b_blk_no
        INTEGER(int64) :: varz1,varz2, tar_z, loc_z
        INTEGER(int64), INTENT(IN) :: g

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

        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                tar_y=l
                tar_z=q
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans

                end if
                END DO
                END DO
                END DO
                END DO
                END DO
                END DO

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
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                tar_y=l
                tar_z=q
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
                bl_intp_valx=block(b_blk_no)%xu(tar_x)
                bl_intp_valy=block(b_blk_no)%yu(tar_y)
                bl_intp_valz=block(b_blk_no)%zu(tar_z)
                bl_intp_x1=block(a_blk_no)%xu(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yu(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zu(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zu(loc_z+1)
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
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans
        end if
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO

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
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                tar_y=l
                tar_z=q
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
                bl_intp_valx=block(b_blk_no)%xv(tar_x)
                bl_intp_valy=block(b_blk_no)%yv(tar_y)
                bl_intp_valz=block(b_blk_no)%zv(tar_z)
                bl_intp_x1=block(a_blk_no)%xv(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yv(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zv(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zv(loc_z+1)
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
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans
        end if
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO

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
        DO q=varz1,varz2
        DO l=vary1,vary2
        DO s=varx1,varx2
                tar_x=s
                tar_y=l
                tar_z=q
               if (block(b_blk_no)%cell_n(tar_x,tar_y,tar_z) /=1)then
                bl_intp_valx=block(b_blk_no)%xw(tar_x)
                bl_intp_valy=block(b_blk_no)%yw(tar_y)
                bl_intp_valz=block(b_blk_no)%zw(tar_z)
                bl_intp_x1=block(a_blk_no)%xw(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xw(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yw(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yw(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zw(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zw(loc_z+1)
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
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*&
                     ((bl_interp_ans2 - bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans
        end if
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO

        end subroutine fineUpdate_mv
end module biocfd_fine_interp

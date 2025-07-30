module biocfd_fine_interp_bound
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  ! allow(use-all) - TODO: Aim to fix this in the future
  use global
  implicit none

  private

  public :: fineUpdate_bd, fineUpdate_pc_bd, fineUpdate_newv_bd, fineUpdate_bd_mv

  contains
SUBROUTINE fineUpdate_bd
        INTEGER(int64) :: g

        DO g=1,intflines
           call fineUpdate_bd_mv(g)
        ENDDO

      end subroutine fineUpdate_bd

        SUBROUTINE fineUpdate_pc_bd
        REAL (dp) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,&
             bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (dp) :: bl_intp_valz,bl_intp_z1,bl_intp_z2
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, &
             bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp) :: bl_interp_ans, bl_interp_ans1, bl_interp_ans2

        INTEGER(int64) :: i,j,k, varx1,varx2, vary1, vary2, tar_x, tar_y, loc_x, &
             loc_y,g, a_blk_no, b_blk_no
        INTEGER(int64) :: varz1,varz2, tar_z, loc_z
        !> axis and steps control how the looping is performed over the x, y, and z axes
        integer :: axis, steps(3)

        DO g=1,intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk

        ! Loop through the x (1), y (2), and z (3) axes
        DO axis=1, 3

          ! Set all steps to 1
          steps = 1

          ! If intfr%counter[xyz]p was a single variable this would be nicer
          if (axis == 1) then
               steps(1) = intfr(g)%counterxp-1
          else if (axis == 2) then
               steps(2) = intfr(g)%counteryp-1
          else if (axis == 3) then
               steps(3) = intfr(g)%counterzp-1
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ppppppp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k=1, intfr(g)%counterzp, steps(3)
        DO j=1, intfr(g)%counteryp, steps(2)
        DO i=1, intfr(g)%counterxp, steps(1)

        varx1=intfr(g)%px_interface_det(2,i)
        varx2=intfr(g)%px_interface_det(3,i)
        vary1=intfr(g)%py_interface_det(2,j)
        vary2=intfr(g)%py_interface_det(3,j)
        varz1=intfr(g)%pz_interface_det(2,k)
        varz2=intfr(g)%pz_interface_det(3,k)

        loc_x=intfr(g)%px_interface_det(1,i)
        loc_y=intfr(g)%py_interface_det(1,j)
        loc_z=intfr(g)%pz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

                bl_intp_valx=block(b_blk_no)%xp(tar_x)
                bl_intp_valy=block(b_blk_no)%yp(tar_y)
                bl_intp_valz=block(b_blk_no)%zp(tar_z)
                bl_intp_x1=block(a_blk_no)%xp(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xp(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yp(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yp(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zp(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zp(loc_z+1)
                bl_intp_f1=block(a_blk_no)%pc(loc_x-1,loc_y-1,loc_z-1)
                bl_intp_f2=block(a_blk_no)%pc(loc_x+1,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%pc(loc_x+1,loc_y+1,loc_z-1)
                bl_intp_f4=block(a_blk_no)%pc(loc_x-1,loc_y+1,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%pc(loc_x-1,loc_y-1,loc_z+1)
                bl_intp_f2=block(a_blk_no)%pc(loc_x+1,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%pc(loc_x+1,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%pc(loc_x-1,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%pc(tar_x,tar_y,tar_z)=bl_interp_ans
                block(b_blk_no)%pco(tar_x,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

     end do ! axes loop

        ENDDO


        end subroutine fineUpdate_pc_bd

        SUBROUTINE fineUpdate_newv_bd
        REAL (dp) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,&
             bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (dp) :: bl_intp_valz,bl_intp_z1,bl_intp_z2
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, &
             bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp) :: bl_interp_ans, bl_interp_ans1, bl_interp_ans2

        INTEGER(int64) :: i,j,k, varx1,varx2, vary1, vary2, tar_x, tar_y, loc_x, &
             loc_y,g, a_blk_no, b_blk_no
        INTEGER(int64) :: varz1,varz2, tar_z, loc_z
        !> axis and steps control how the looping is performed over the x, y, and z axes
        integer :: axis, steps(3)


        DO g=1,intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!uuuuuuu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Loop through the x (1), y (2), and z (3) axes
          DO axis=1, 3

          ! Set all steps to 1
          steps = 1

          if (axis == 1) then
               steps(1) = intfr(g)%counterxu-1
          else if (axis == 2) then
               steps(2) = intfr(g)%counteryu-1
          else if (axis == 3) then
               steps(3) = intfr(g)%counterzu-1
          end if

        DO k=1, intfr(g)%counterzu, steps(3)
        DO j=1, intfr(g)%counteryu, steps(2)
        DO i=1, intfr(g)%counterxu, steps(1)



        varx1=intfr(g)%ux_interface_det(2,i)
        varx2=intfr(g)%ux_interface_det(3,i)
        vary1=intfr(g)%uy_interface_det(2,j)
        vary2=intfr(g)%uy_interface_det(3,j)
        varz1=intfr(g)%uz_interface_det(2,k)
        varz2=intfr(g)%uz_interface_det(3,k)

        loc_x=intfr(g)%ux_interface_det(1,i)
        loc_y=intfr(g)%uy_interface_det(1,j)
        loc_z=intfr(g)%uz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

                bl_intp_valx=block(b_blk_no)%xu(tar_x)
                bl_intp_valy=block(b_blk_no)%yu(tar_y)
                bl_intp_valz=block(b_blk_no)%zu(tar_z)
                bl_intp_x1=block(a_blk_no)%xu(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xu(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yu(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yu(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zu(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zu(loc_z+1)
                bl_intp_f1=block(a_blk_no)%ut(loc_x-2,loc_y-1,loc_z-1)
                bl_intp_f2=block(a_blk_no)%ut(loc_x,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%ut(loc_x,loc_y+1,loc_z-1)
                bl_intp_f4=block(a_blk_no)%ut(loc_x-2,loc_y+1,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%ut(loc_x-2,loc_y-1,loc_z+1)
                bl_intp_f2=block(a_blk_no)%ut(loc_x,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%ut(loc_x,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%ut(loc_x-2,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%ut(tar_x-1,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

     end do  ! axes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!vvvvvvv!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Loop through the x (1), y (2), and z (3) axes
        DO axis=1, 3

        ! Set all steps to 1
        steps = 1

        if (axis == 1) then
             steps(1) = intfr(g)%counterxv-1
        else if (axis == 2) then
             steps(2) = intfr(g)%counteryv-1
        else if (axis == 3) then
             steps(3) = intfr(g)%counterzv-1
        end if

        DO k=1, intfr(g)%counterzv, steps(3)
        DO j=1, intfr(g)%counteryv, steps(2)
        DO i=1, intfr(g)%counterxv, steps(1)

        varx1=intfr(g)%vx_interface_det(2,i)
        varx2=intfr(g)%vx_interface_det(3,i)
        vary1=intfr(g)%vy_interface_det(2,j)
        vary2=intfr(g)%vy_interface_det(3,j)
        varz1=intfr(g)%vz_interface_det(2,k)
        varz2=intfr(g)%vz_interface_det(3,k)

        loc_x=intfr(g)%vx_interface_det(1,i)
        loc_y=intfr(g)%vy_interface_det(1,j)
        loc_z=intfr(g)%vz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

                bl_intp_valx=block(b_blk_no)%xv(tar_x)
                bl_intp_valy=block(b_blk_no)%yv(tar_y)
                bl_intp_valz=block(b_blk_no)%zv(tar_z)
                bl_intp_x1=block(a_blk_no)%xv(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xv(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yv(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yv(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zv(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zv(loc_z+1)
                bl_intp_f1=block(a_blk_no)%vt(loc_x-1,loc_y-2,loc_z-1)
                bl_intp_f2=block(a_blk_no)%vt(loc_x+1,loc_y-2,loc_z-1)
                bl_intp_f3=block(a_blk_no)%vt(loc_x+1,loc_y  ,loc_z-1)
                bl_intp_f4=block(a_blk_no)%vt(loc_x-1,loc_y  ,loc_z-1)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%vt(loc_x-1,loc_y-2,loc_z+1)
                bl_intp_f2=block(a_blk_no)%vt(loc_x+1,loc_y-2,loc_z+1)
                bl_intp_f3=block(a_blk_no)%vt(loc_x+1,loc_y  ,loc_z+1)
                bl_intp_f4=block(a_blk_no)%vt(loc_x-1,loc_y  ,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%vt(tar_x,tar_y-1,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

     end do  ! axes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!wwwwwww!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Loop through the x (1), y (2), and z (3) axes
        DO axis=1, 3

        ! Set all steps to 1
        steps = 1

        if (axis == 1) then
             steps(1) = intfr(g)%counterxw-1
        else if (axis == 2) then
             steps(2) = intfr(g)%counteryw-1
        else if (axis == 3) then
             steps(3) = intfr(g)%counterzw-1
        end if

        DO k=1, intfr(g)%counterzw, steps(3)
        DO j=1, intfr(g)%counteryw, steps(2)
        DO i=1, intfr(g)%counterxw, steps(1)

        varx1=intfr(g)%wx_interface_det(2,i)
        varx2=intfr(g)%wx_interface_det(3,i)
        vary1=intfr(g)%wy_interface_det(2,j)
        vary2=intfr(g)%wy_interface_det(3,j)
        varz1=intfr(g)%wz_interface_det(2,k)
        varz2=intfr(g)%wz_interface_det(3,k)

        loc_x=intfr(g)%wx_interface_det(1,i)
        loc_y=intfr(g)%wy_interface_det(1,j)
        loc_z=intfr(g)%wz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

                bl_intp_valx=block(b_blk_no)%xw(tar_x)
                bl_intp_valy=block(b_blk_no)%yw(tar_y)
                bl_intp_valz=block(b_blk_no)%zw(tar_z)
                bl_intp_x1=block(a_blk_no)%xw(loc_x-1)
                bl_intp_x2=block(a_blk_no)%xw(loc_x+1)
                bl_intp_y1=block(a_blk_no)%yw(loc_y-1)
                bl_intp_y2=block(a_blk_no)%yw(loc_y+1)
                bl_intp_z1=block(a_blk_no)%zw(loc_z-1)
                bl_intp_z2=block(a_blk_no)%zw(loc_z+1)
                bl_intp_f1=block(a_blk_no)%wt(loc_x-1,loc_y-1,loc_z-2)
                bl_intp_f2=block(a_blk_no)%wt(loc_x+1,loc_y-1,loc_z-2)
                bl_intp_f3=block(a_blk_no)%wt(loc_x+1,loc_y+1,loc_z-2)
                bl_intp_f4=block(a_blk_no)%wt(loc_x-1,loc_y+1,loc_z-2)

                bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
                bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
                bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
                bl_intp_yty=(bl_intp_y2-bl_intp_valy)
                bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term

                bl_interp_ans1=(bl_intp_num/bl_intp_deno)
                bl_intp_f1=block(a_blk_no)%wt(loc_x-1,loc_y-1,loc_z)
                bl_intp_f2=block(a_blk_no)%wt(loc_x+1,loc_y-1,loc_z)
                bl_intp_f3=block(a_blk_no)%wt(loc_x+1,loc_y+1,loc_z)
                bl_intp_f4=block(a_blk_no)%wt(loc_x-1,loc_y+1,loc_z)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%wt(tar_x,tar_y,tar_z-1)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

     end do  ! axes

        ENDDO


        end subroutine fineUpdate_newv_bd

        SUBROUTINE fineUpdate_bd_mv(g)
        REAL (dp) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,&
             bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (dp) :: bl_intp_valz,bl_intp_z1,bl_intp_z2
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, &
             bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp) :: bl_interp_ans, bl_interp_ans1, bl_interp_ans2

        INTEGER(int64) :: i,j,k, varx1,varx2, vary1, vary2, tar_x, tar_y, loc_x, &
             loc_y, a_blk_no, b_blk_no
        INTEGER(int64) :: varz1,varz2, tar_z, loc_z
         INTEGER (int64), INTENT(IN) :: g

           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ppppppp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k=1,intfr(g)%counterzp
        DO j=1,intfr(g)%counteryp
        DO i=1,intfr(g)%counterxp, intfr(g)%counterxp-1

        varx1=intfr(g)%px_interface_det(2,i)
        varx2=intfr(g)%px_interface_det(3,i)
        vary1=intfr(g)%py_interface_det(2,j)
        vary2=intfr(g)%py_interface_det(3,j)
        varz1=intfr(g)%pz_interface_det(2,k)
        varz2=intfr(g)%pz_interface_det(3,k)

        loc_x=intfr(g)%px_interface_det(1,i)
        loc_y=intfr(g)%py_interface_det(1,j)
        loc_z=intfr(g)%pz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1,intfr(g)%counterzp
        DO j=1,intfr(g)%counteryp, intfr(g)%counteryp-1
        DO i=1,intfr(g)%counterxp

        varx1=intfr(g)%px_interface_det(2,i)
        varx2=intfr(g)%px_interface_det(3,i)
        vary1=intfr(g)%py_interface_det(2,j)
        vary2=intfr(g)%py_interface_det(3,j)
        varz1=intfr(g)%pz_interface_det(2,k)
        varz2=intfr(g)%pz_interface_det(3,k)

        loc_x=intfr(g)%px_interface_det(1,i)
        loc_y=intfr(g)%py_interface_det(1,j)
        loc_z=intfr(g)%pz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1,intfr(g)%counterzp, intfr(g)%counterzp-1
        DO j=1,intfr(g)%counteryp
        DO i=1,intfr(g)%counterxp

        varx1=intfr(g)%px_interface_det(2,i)
        varx2=intfr(g)%px_interface_det(3,i)
        vary1=intfr(g)%py_interface_det(2,j)
        vary2=intfr(g)%py_interface_det(3,j)
        varz1=intfr(g)%pz_interface_det(2,k)
        varz2=intfr(g)%pz_interface_det(3,k)

        loc_x=intfr(g)%px_interface_det(1,i)
        loc_y=intfr(g)%py_interface_det(1,j)
        loc_z=intfr(g)%pz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%p(tar_x,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!uuuuuuu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k=1,intfr(g)%counterzu
        DO j=1,intfr(g)%counteryu
        DO i=1,intfr(g)%counterxu, intfr(g)%counterxu-1

        varx1=intfr(g)%ux_interface_det(2,i)
        varx2=intfr(g)%ux_interface_det(3,i)
        vary1=intfr(g)%uy_interface_det(2,j)
        vary2=intfr(g)%uy_interface_det(3,j)
        varz1=intfr(g)%uz_interface_det(2,k)
        varz2=intfr(g)%uz_interface_det(3,k)

        loc_x=intfr(g)%ux_interface_det(1,i)
        loc_y=intfr(g)%uy_interface_det(1,j)
        loc_z=intfr(g)%uz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f2=block(a_blk_no)%u(loc_x  ,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%u(loc_x  ,loc_y+1,loc_z-1)
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
                bl_intp_f2=block(a_blk_no)%u(loc_x  ,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%u(loc_x  ,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%u(loc_x-2,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1,intfr(g)%counterzu
        DO j=1,intfr(g)%counteryu, intfr(g)%counteryu-1
        DO i=1,intfr(g)%counterxu

        varx1=intfr(g)%ux_interface_det(2,i)
        varx2=intfr(g)%ux_interface_det(3,i)
        vary1=intfr(g)%uy_interface_det(2,j)
        vary2=intfr(g)%uy_interface_det(3,j)
        varz1=intfr(g)%uz_interface_det(2,k)
        varz2=intfr(g)%uz_interface_det(3,k)

        loc_x=intfr(g)%ux_interface_det(1,i)
        loc_y=intfr(g)%uy_interface_det(1,j)
        loc_z=intfr(g)%uz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f2=block(a_blk_no)%u(loc_x  ,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%u(loc_x  ,loc_y+1,loc_z-1)
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
                bl_intp_f2=block(a_blk_no)%u(loc_x  ,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%u(loc_x  ,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%u(loc_x-2,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1,intfr(g)%counterzu, intfr(g)%counterzu-1
        DO j=1,intfr(g)%counteryu
        DO i=1,intfr(g)%counterxu

        varx1=intfr(g)%ux_interface_det(2,i)
        varx2=intfr(g)%ux_interface_det(3,i)
        vary1=intfr(g)%uy_interface_det(2,j)
        vary2=intfr(g)%uy_interface_det(3,j)
        varz1=intfr(g)%uz_interface_det(2,k)
        varz2=intfr(g)%uz_interface_det(3,k)

        loc_x=intfr(g)%ux_interface_det(1,i)
        loc_y=intfr(g)%uy_interface_det(1,j)
        loc_z=intfr(g)%uz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f2=block(a_blk_no)%u(loc_x  ,loc_y-1,loc_z-1)
                bl_intp_f3=block(a_blk_no)%u(loc_x  ,loc_y+1,loc_z-1)
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
                bl_intp_f2=block(a_blk_no)%u(loc_x  ,loc_y-1,loc_z+1)
                bl_intp_f3=block(a_blk_no)%u(loc_x  ,loc_y+1,loc_z+1)
                bl_intp_f4=block(a_blk_no)%u(loc_x-2,loc_y+1,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%u(tar_x-1,tar_y,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!vvvvvvv!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k=1,intfr(g)%counterzv
        DO j=1,intfr(g)%counteryv
        DO i=1,intfr(g)%counterxv, intfr(g)%counterxv-1

        varx1=intfr(g)%vx_interface_det(2,i)
        varx2=intfr(g)%vx_interface_det(3,i)
        vary1=intfr(g)%vy_interface_det(2,j)
        vary2=intfr(g)%vy_interface_det(3,j)
        varz1=intfr(g)%vz_interface_det(2,k)
        varz2=intfr(g)%vz_interface_det(3,k)

        loc_x=intfr(g)%vx_interface_det(1,i)
        loc_y=intfr(g)%vy_interface_det(1,j)
        loc_z=intfr(g)%vz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y  ,loc_z-1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y  ,loc_z-1)

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
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y  ,loc_z+1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y  ,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1,intfr(g)%counterzv
        DO j=1,intfr(g)%counteryv, intfr(g)%counteryv-1
        DO i=1,intfr(g)%counterxv

        varx1=intfr(g)%vx_interface_det(2,i)
        varx2=intfr(g)%vx_interface_det(3,i)
        vary1=intfr(g)%vy_interface_det(2,j)
        vary2=intfr(g)%vy_interface_det(3,j)
        varz1=intfr(g)%vz_interface_det(2,k)
        varz2=intfr(g)%vz_interface_det(3,k)

        loc_x=intfr(g)%vx_interface_det(1,i)
        loc_y=intfr(g)%vy_interface_det(1,j)
        loc_z=intfr(g)%vz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y  ,loc_z-1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y  ,loc_z-1)

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
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y  ,loc_z+1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y  ,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1,intfr(g)%counterzv, intfr(g)%counterzv-1
        DO j=1,intfr(g)%counteryv
        DO i=1,intfr(g)%counterxv

        varx1=intfr(g)%vx_interface_det(2,i)
        varx2=intfr(g)%vx_interface_det(3,i)
        vary1=intfr(g)%vy_interface_det(2,j)
        vary2=intfr(g)%vy_interface_det(3,j)
        varz1=intfr(g)%vz_interface_det(2,k)
        varz2=intfr(g)%vz_interface_det(3,k)

        loc_x=intfr(g)%vx_interface_det(1,i)
        loc_y=intfr(g)%vy_interface_det(1,j)
        loc_z=intfr(g)%vz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y  ,loc_z-1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y  ,loc_z-1)

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
                bl_intp_f3=block(a_blk_no)%v(loc_x+1,loc_y  ,loc_z+1)
                bl_intp_f4=block(a_blk_no)%v(loc_x-1,loc_y  ,loc_z+1)
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%v(tar_x,tar_y-1,tar_z)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!wwwwwww!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k=1,intfr(g)%counterzw
        DO j=1,intfr(g)%counteryw
        DO i=1,intfr(g)%counterxw, intfr(g)%counterxw-1

        varx1=intfr(g)%wx_interface_det(2,i)
        varx2=intfr(g)%wx_interface_det(3,i)
        vary1=intfr(g)%wy_interface_det(2,j)
        vary2=intfr(g)%wy_interface_det(3,j)
        varz1=intfr(g)%wz_interface_det(2,k)
        varz2=intfr(g)%wz_interface_det(3,k)

        loc_x=intfr(g)%wx_interface_det(1,i)
        loc_y=intfr(g)%wy_interface_det(1,j)
        loc_z=intfr(g)%wz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f1=block(a_blk_no)%w(loc_x-1,loc_y-1,loc_z  )
                bl_intp_f2=block(a_blk_no)%w(loc_x+1,loc_y-1,loc_z  )
                bl_intp_f3=block(a_blk_no)%w(loc_x+1,loc_y+1,loc_z  )
                bl_intp_f4=block(a_blk_no)%w(loc_x-1,loc_y+1,loc_z  )
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1,intfr(g)%counterzw
        DO j=1,intfr(g)%counteryw, intfr(g)%counteryw-1
        DO i=1,intfr(g)%counterxw

        varx1=intfr(g)%wx_interface_det(2,i)
        varx2=intfr(g)%wx_interface_det(3,i)
        vary1=intfr(g)%wy_interface_det(2,j)
        vary2=intfr(g)%wy_interface_det(3,j)
        varz1=intfr(g)%wz_interface_det(2,k)
        varz2=intfr(g)%wz_interface_det(3,k)

        loc_x=intfr(g)%wx_interface_det(1,i)
        loc_y=intfr(g)%wy_interface_det(1,j)
        loc_z=intfr(g)%wz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f1=block(a_blk_no)%w(loc_x-1,loc_y-1,loc_z  )
                bl_intp_f2=block(a_blk_no)%w(loc_x+1,loc_y-1,loc_z  )
                bl_intp_f3=block(a_blk_no)%w(loc_x+1,loc_y+1,loc_z  )
                bl_intp_f4=block(a_blk_no)%w(loc_x-1,loc_y+1,loc_z  )
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        DO k=1, intfr(g)%counterzw, intfr(g)%counterzw-1
        DO j=1,intfr(g)%counteryw
        DO i=1,intfr(g)%counterxw

        varx1=intfr(g)%wx_interface_det(2,i)
        varx2=intfr(g)%wx_interface_det(3,i)
        vary1=intfr(g)%wy_interface_det(2,j)
        vary2=intfr(g)%wy_interface_det(3,j)
        varz1=intfr(g)%wz_interface_det(2,k)
        varz2=intfr(g)%wz_interface_det(3,k)

        loc_x=intfr(g)%wx_interface_det(1,i)
        loc_y=intfr(g)%wy_interface_det(1,j)
        loc_z=intfr(g)%wz_interface_det(1,k)

        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

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
                bl_intp_f1=block(a_blk_no)%w(loc_x-1,loc_y-1,loc_z  )
                bl_intp_f2=block(a_blk_no)%w(loc_x+1,loc_y-1,loc_z  )
                bl_intp_f3=block(a_blk_no)%w(loc_x+1,loc_y+1,loc_z  )
                bl_intp_f4=block(a_blk_no)%w(loc_x-1,loc_y+1,loc_z  )
                bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
                bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
                bl_intp_num= bl_intp_first_term + bl_intp_second_term
                bl_interp_ans2=(bl_intp_num/bl_intp_deno)
                bl_interp_ans= bl_interp_ans1 + ((bl_intp_valz-bl_intp_z1)*((bl_interp_ans2 - &
                     bl_interp_ans1)/(bl_intp_z2 - bl_intp_z1)))
                block(b_blk_no)%w(tar_x,tar_y,tar_z-1)=bl_interp_ans
                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO

        end subroutine fineUpdate_bd_mv
end module biocfd_fine_interp_bound

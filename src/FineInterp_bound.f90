module biocfd_fine_interp_bound
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, intfr, intflines
  use biocfd_interpolation, only: bilinear_interpolation, linear_interpolation
  use biocfd_blocks, only: Blocks
  use biocfd_interface_type, only: Interface_t

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

      SUBROUTINE fineUpdate_pc_bd(local_intfr,blk_a,blk_b)
        type(Interface_t),intent(in) :: local_intfr
        type(Blocks), intent(inout) :: blk_b
        type(Blocks), intent(in) :: blk_a

        REAL (dp) :: bl_interp_ans
        INTEGER(int64) :: i,j,k, varx1,varx2, vary1, vary2, tar_x, tar_y, loc_x, &
             loc_y
        INTEGER(int64) :: varz1,varz2, tar_z, loc_z
        !> axis and steps control how the looping is performed over the x, y, and z axes
        integer :: axis, steps(3)

        ! Loop through the x (1), y (2), and z (3) axes
        DO axis=1, 3

          ! Set all steps to 1
          steps = 1

          ! If intfr%counter[xyz]p was a single variable this would be nicer
          if (axis == 1) then
               steps(1) = local_intfr%counterxp-1
          else if (axis == 2) then
               steps(2) = local_intfr%counteryp-1
          else if (axis == 3) then
               steps(3) = local_intfr%counterzp-1
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ppppppp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z, bl_interp_ans) &
        !$acc firstprivate(a_blk_no, b_blk_no)
        DO k=1, local_intfr%counterzp, steps(3)
        DO j=1, local_intfr%counteryp, steps(2)
        DO i=1, local_intfr%counterxp, steps(1)

        varx1=local_intfr%px_interface_det(2,i)
        varx2=local_intfr%px_interface_det(3,i)
        vary1=local_intfr%py_interface_det(2,j)
        vary2=local_intfr%py_interface_det(3,j)
        varz1=local_intfr%pz_interface_det(2,k)
        varz2=local_intfr%pz_interface_det(3,k)

        loc_x=local_intfr%px_interface_det(1,i)
        loc_y=local_intfr%py_interface_det(1,j)
        loc_z=local_intfr%pz_interface_det(1,k)

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
          DO tar_y=vary1,vary2
            DO tar_x=varx1,varx2

               bl_interp_ans = trilinear_interpolation(&
                 blk_b%xp(tar_x), blk_b%yp(tar_y), blk_b%zp(tar_z), &
                 loc_x, loc_y, loc_z, blk_a%xp, blk_a%yp, blk_a%zp, &
                 0, blk_a%pc &
               )

                blk_b%pc(tar_x, tar_y, tar_z) = bl_interp_ans
                blk_b%pco(tar_x, tar_y, tar_z) = bl_interp_ans

            ENDDO
          ENDDO
        ENDDO

        ENDDO
        ENDDO
        ENDDO
        !$acc end parallel loop

        end do  ! axes loop

        end subroutine fineUpdate_pc_bd

        SUBROUTINE fineUpdate_newv_bd

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

        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z) &
        !$acc firstprivate(a_blk_no, b_blk_no)
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

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

          block(b_blk_no)%ut(tar_x-1,tar_y,tar_z) = trilinear_interpolation(&
               block(b_blk_no)%xu(tar_x), block(b_blk_no)%yu(tar_y), block(b_blk_no)%zu(tar_z), &
               loc_x, loc_y, loc_z, block(a_blk_no)%xu, block(a_blk_no)%yu, block(a_blk_no)%zu, &
               1, block(a_blk_no)%ut &
          )

          ENDDO
          ENDDO
          ENDDO

        ENDDO
        ENDDO
        ENDDO
        !$acc end parallel loop

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

        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z) &
        !$acc firstprivate(a_blk_no, b_blk_no)
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

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

          block(b_blk_no)%vt(tar_x,tar_y-1,tar_z) = trilinear_interpolation(&
               block(b_blk_no)%xv(tar_x), block(b_blk_no)%yv(tar_y), block(b_blk_no)%zv(tar_z), &
               loc_x, loc_y, loc_z, block(a_blk_no)%xv, block(a_blk_no)%yv, block(a_blk_no)%zv, &
               2, block(a_blk_no)%vt &
          )

          ENDDO
          ENDDO
          ENDDO

        ENDDO
        ENDDO
        ENDDO
        !$acc end parallel loop

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

        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z) &
        !$acc firstprivate(a_blk_no, b_blk_no)
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

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

          block(b_blk_no)%wt(tar_x,tar_y,tar_z-1) = trilinear_interpolation(&
               block(b_blk_no)%xw(tar_x), block(b_blk_no)%yw(tar_y), block(b_blk_no)%zw(tar_z), &
               loc_x, loc_y, loc_z, block(a_blk_no)%xw, block(a_blk_no)%yw, block(a_blk_no)%zw, &
               3, block(a_blk_no)%wt &
          )

          ENDDO
          ENDDO
          ENDDO

        ENDDO
        ENDDO
        ENDDO
        !$acc end parallel loop

     end do  ! axes

        ENDDO


        end subroutine fineUpdate_newv_bd

        SUBROUTINE fineUpdate_bd_mv(g)

        INTEGER(int64) :: i,j,k, varx1,varx2, vary1, vary2, tar_x, tar_y, loc_x, &
             loc_y, a_blk_no, b_blk_no
        INTEGER(int64) :: varz1,varz2, tar_z, loc_z
         INTEGER (int64), INTENT(IN) :: g
         !> axis and steps control how the looping is performed over the x, y, and z axes
        integer :: axis, steps(3)

           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ppppppp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Loop through the x (1), y (2), and z (3) axes
        DO axis=1, 3

        ! Set all steps to 1
        steps = 1

        if (axis == 1) then
             steps(1) = intfr(g)%counterxp-1
        else if (axis == 2) then
             steps(2) = intfr(g)%counteryp-1
        else if (axis == 3) then
             steps(3) = intfr(g)%counterzp-1
        end if


        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z) &
        !$acc firstprivate(a_blk_no, b_blk_no)
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

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

          block(b_blk_no)%p(tar_x,tar_y,tar_z) = trilinear_interpolation(&
               block(b_blk_no)%xp(tar_x), block(b_blk_no)%yp(tar_y), block(b_blk_no)%zp(tar_z), &
               loc_x, loc_y, loc_z, block(a_blk_no)%xp, block(a_blk_no)%yp, block(a_blk_no)%zp, &
               0, block(a_blk_no)%p &
          )

                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO
         !$acc end parallel loop

        end do  ! axes

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

        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z) &
        !$acc firstprivate(a_blk_no, b_blk_no)
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

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

          block(b_blk_no)%u(tar_x-1,tar_y,tar_z) = trilinear_interpolation(&
               block(b_blk_no)%xu(tar_x), block(b_blk_no)%yu(tar_y), block(b_blk_no)%zu(tar_z), &
               loc_x, loc_y, loc_z, block(a_blk_no)%xu, block(a_blk_no)%yu, block(a_blk_no)%zu, &
               1, block(a_blk_no)%u &
          )

                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO
         !$acc end parallel loop

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

        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z) &
        !$acc firstprivate(a_blk_no, b_blk_no)
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

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

          block(b_blk_no)%v(tar_x,tar_y-1,tar_z) = trilinear_interpolation(&
               block(b_blk_no)%xv(tar_x), block(b_blk_no)%yv(tar_y), block(b_blk_no)%zv(tar_z), &
               loc_x, loc_y, loc_z, block(a_blk_no)%xv, block(a_blk_no)%yv, block(a_blk_no)%zv, &
               2, block(a_blk_no)%v &
          )

                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO
         !$acc end parallel loop

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

        !$acc parallel loop collapse(3) private(varx1, varx2, vary1, vary2, varz1, varz2) &
        !$acc private(loc_x, loc_y, loc_z) &
        !$acc firstprivate(a_blk_no, b_blk_no)
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

        !$acc loop collapse(3) seq
        DO tar_z=varz1,varz2
        DO tar_y=vary1,vary2
        DO tar_x=varx1,varx2

          block(b_blk_no)%w(tar_x,tar_y,tar_z-1) = trilinear_interpolation(&
               block(b_blk_no)%xw(tar_x), block(b_blk_no)%yw(tar_y), block(b_blk_no)%zw(tar_z), &
               loc_x, loc_y, loc_z, block(a_blk_no)%xw, block(a_blk_no)%yw, block(a_blk_no)%zw, &
               3, block(a_blk_no)%w &
          )

                ENDDO
                ENDDO
                ENDDO

        ENDDO
        ENDDO
        ENDDO
        !$acc end parallel loop

        end do  ! axes

        end subroutine fineUpdate_bd_mv

      !> Perform trilinear interpolation. Implemented as two bilinear
      !> intepolations followed by a linear interpolation of the
      !> results. See
      !> https://en.wikipedia.org/wiki/Trilinear_interpolation.
      !>
      !> WARNING: that this lives in this module and not in the
      !> Interpolation module as it contains some specific logic for
      !> this module, namely the grid indices
      pure function trilinear_interpolation(x, y, z, i, j, k, xgrid, ygrid, zgrid, offset, var) &
           result(out)

        !> The target position of the trilinear interpolation
        real(dp), intent(in) :: x, y, z
        !> Indices used to determine the grid locations
        integer(int64), intent(in) :: i, j, k
        !> The grids
        real(dp), intent(in), dimension(:) :: xgrid, ygrid, zgrid
        !> The offset value changes depending on whether we are computing p, u, v, or w
        !> p=0, u=1, v=2, w=3
        integer, intent(in) :: offset
        !> The variable that is to be interpolated
        real(dp), intent(in), dimension(:, :, :) :: var

        !> The interpolated result
        real(dp) :: out
        ! Intermediate results used in the calculation
        real(dp) :: z1, z2

        ! Indices used to determine the variable indices after offsets
        integer :: vi, vj, vk

        vi = i
        vj = j
        vk = k

        ! Set the indices appropiately
        if (offset == 0) then
           ! Leave everything as is
        else if (offset == 1) then
           vi = i - 1
        else if (offset == 2) then
           vj = j - 1
        else if (offset == 3) then
           vk = k - 1
        else
           ! Anything else is an error, but how to handle it?
        end if

        z1 = bilinear_interpolation(x, y, xgrid(i-1), xgrid(i+1), ygrid(j-1), ygrid(j+1), &
                                   [var(vi-1, vj-1, vk-1), var(vi+1, vj-1, vk-1), &
                                   var(vi-1, vj+1, vk-1), var(vi+1, vj+1, vk-1)])

        z2 = bilinear_interpolation(x, y, xgrid(i-1), xgrid(i+1), ygrid(j-1), ygrid(j+1), &
                                   [var(vi-1, vj-1, vk+1), var(vi+1, vj-1, vk+1), &
                                   var(vi-1, vj+1, vk+1), var(vi+1, vj+1, vk+1)])

        out = linear_interpolation(z, zgrid(k-1), zgrid(k+1), z1, z2)
      end function trilinear_interpolation
end module biocfd_fine_interp_bound

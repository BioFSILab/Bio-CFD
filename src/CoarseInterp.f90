module biocfd_coarse_update
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use global, only : block, intfr
  use biocfd_interpolation, only: bilinear_interpolation, linear_interpolation

  implicit none
  private

  public ::  coarseUpdate_newv, coarseUpdate_pc, coarseUpdate

  contains
subroutine coarseUpdate
        INTEGER :: i,j,k,g, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz

        DO g=1,size(intfr)
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxp-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryp-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzp-block(b_blk_no)%cintp

        !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
        do k=st_idz,en_idz
          do j=st_idy,en_idy
            do i=st_idx, en_idx
              tar_x = intfr(g)%px_interface_det(1,i)
              tar_y = intfr(g)%py_interface_det(1,j)
              tar_z = intfr(g)%pz_interface_det(1,k)
              loc_x = (intfr(g)%px_interface_det(2,i) + intfr(g)%px_interface_det(3,i)-1) / 2
              loc_y = (intfr(g)%py_interface_det(2,j) + intfr(g)%py_interface_det(3,j)-1) / 2
              loc_z = (intfr(g)%pz_interface_det(2,k) + intfr(g)%pz_interface_det(3,k)-1) / 2

              block(a_blk_no)%p(tar_x, tar_y, tar_z) = trilinear_interpolation( &
                  block(a_blk_no)%xp(tar_x), block(a_blk_no)%yp(tar_y), block(a_blk_no)%zp(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xp, block(b_blk_no)%yp, block(b_blk_no)%zp, &
                  0, block(b_blk_no)%p)

            enddo
          enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxu-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryu-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzu-block(b_blk_no)%cintp
       !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%ux_interface_det(1,i)
                tar_y=intfr(g)%uy_interface_det(1,j)
                tar_z=intfr(g)%uz_interface_det(1,k)
                loc_x=intfr(g)%ux_interface_det(3,i)
                loc_y=(intfr(g)%uy_interface_det(2,j)+intfr(g)%uy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%uz_interface_det(2,k)+intfr(g)%uz_interface_det(3,k)-1)/2

                block(a_blk_no)%u(tar_x-1,tar_y,tar_z) = trilinear_interpolation( &
                  block(a_blk_no)%xu(tar_x), block(a_blk_no)%yu(tar_y), block(a_blk_no)%zu(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xu, block(b_blk_no)%yu, block(b_blk_no)%zu, &
                  1, block(b_blk_no)%u)

        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxv-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryv-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzv-block(b_blk_no)%cintp
       !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%vx_interface_det(1,i)
                tar_y=intfr(g)%vy_interface_det(1,j)
                tar_z=intfr(g)%vz_interface_det(1,k)
                loc_x=(intfr(g)%vx_interface_det(2,i)+intfr(g)%vx_interface_det(3,i)-1)/2
                loc_z=(intfr(g)%vz_interface_det(2,k)+intfr(g)%vz_interface_det(3,k)-1)/2
                loc_y=(intfr(g)%vy_interface_det(3,j))

                block(a_blk_no)%v(tar_x,tar_y-1,tar_z) = trilinear_interpolation( &
                  block(a_blk_no)%xv(tar_x), block(a_blk_no)%yv(tar_y), block(a_blk_no)%zv(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xv, block(b_blk_no)%yv, block(b_blk_no)%zv, &
                  2, block(b_blk_no)%v)

        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxw-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryw-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzw-block(b_blk_no)%cintp
       !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%wx_interface_det(1,i)
                tar_y=intfr(g)%wy_interface_det(1,j)
                tar_z=intfr(g)%wz_interface_det(1,k)
                loc_x=(intfr(g)%wx_interface_det(2,i)+intfr(g)%wx_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%wy_interface_det(2,j)+intfr(g)%wy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%wz_interface_det(3,k))

                block(a_blk_no)%w(tar_x,tar_y,tar_z-1) = trilinear_interpolation( &
                  block(a_blk_no)%xw(tar_x), block(a_blk_no)%yw(tar_y), block(a_blk_no)%zw(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xw, block(b_blk_no)%yw, block(b_blk_no)%zw, &
                  3, block(b_blk_no)%w)

        enddo
        enddo
        enddo
        ENDDO
        end subroutine coarseUpdate

        subroutine coarseUpdate_newv

        INTEGER :: i,j,k,g, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz

        DO g=1,size(intfr)
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxu-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryu-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzu-block(b_blk_no)%cintp
       !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%ux_interface_det(1,i)
                tar_y=intfr(g)%uy_interface_det(1,j)
                tar_z=intfr(g)%uz_interface_det(1,k)
                loc_x=intfr(g)%ux_interface_det(3,i)
                loc_y=(intfr(g)%uy_interface_det(2,j)+intfr(g)%uy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%uz_interface_det(2,k)+intfr(g)%uz_interface_det(3,k)-1)/2

                block(a_blk_no)%ut(tar_x-1,tar_y,tar_z) = trilinear_interpolation( &
                  block(a_blk_no)%xu(tar_x), block(a_blk_no)%yu(tar_y), block(a_blk_no)%zu(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xu, block(b_blk_no)%yu, block(b_blk_no)%zu, &
                  1, block(b_blk_no)%ut)

        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxv-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryv-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzv-block(b_blk_no)%cintp
       !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%vx_interface_det(1,i)
                tar_y=intfr(g)%vy_interface_det(1,j)
                tar_z=intfr(g)%vz_interface_det(1,k)
                loc_x=(intfr(g)%vx_interface_det(2,i)+intfr(g)%vx_interface_det(3,i)-1)/2
                loc_z=(intfr(g)%vz_interface_det(2,k)+intfr(g)%vz_interface_det(3,k)-1)/2
                loc_y=(intfr(g)%vy_interface_det(3,j))

                block(a_blk_no)%vt(tar_x, tar_y-1, tar_z) = trilinear_interpolation( &
                  block(a_blk_no)%xv(tar_x), block(a_blk_no)%yv(tar_y), block(a_blk_no)%zv(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xv, block(b_blk_no)%yv, block(b_blk_no)%zv, &
                  2, block(b_blk_no)%vt)

        enddo
        enddo
        enddo

        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxw-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryw-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzw-block(b_blk_no)%cintp
       !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%wx_interface_det(1,i)
                tar_y=intfr(g)%wy_interface_det(1,j)
                tar_z=intfr(g)%wz_interface_det(1,k)
                loc_x=(intfr(g)%wx_interface_det(2,i)+intfr(g)%wx_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%wy_interface_det(2,j)+intfr(g)%wy_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%wz_interface_det(3,k))

                block(a_blk_no)%wt(tar_x, tar_y, tar_z-1) = trilinear_interpolation( &
                  block(a_blk_no)%xw(tar_x), block(a_blk_no)%yw(tar_y), block(a_blk_no)%zw(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xw, block(b_blk_no)%yw, block(b_blk_no)%zw, &
                  3, block(b_blk_no)%wt)

        enddo
        enddo
        enddo
        ENDDO
        end subroutine coarseUpdate_newv
        subroutine coarseUpdate_pc

        REAL (dp) :: bl_interp_ans

        INTEGER :: i,j,k,g, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz


        DO g=1,size(intfr)
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        st_idx=block(b_blk_no)%cintp
        en_idx=intfr(g)%counterxp-block(b_blk_no)%cintp
        st_idy=block(b_blk_no)%cintp
        en_idy=intfr(g)%counteryp-block(b_blk_no)%cintp
        st_idz=block(b_blk_no)%cintp
        en_idz=intfr(g)%counterzp-block(b_blk_no)%cintp
       !$acc parallel loop collapse(3) private(tar_x, tar_y, tar_z, loc_x, loc_y, loc_z)
       do k=st_idz,en_idz
       do j=st_idy,en_idy
       do i=st_idx, en_idx
                tar_x=intfr(g)%px_interface_det(1,i)
                tar_y=intfr(g)%py_interface_det(1,j)
                tar_z=intfr(g)%pz_interface_det(1,k)
                loc_x=(intfr(g)%px_interface_det(2,i)+intfr(g)%px_interface_det(3,i)-1)/2
                loc_y=(intfr(g)%py_interface_det(2,j)+intfr(g)%py_interface_det(3,j)-1)/2
                loc_z=(intfr(g)%pz_interface_det(2,k)+intfr(g)%pz_interface_det(3,k)-1)/2

                bl_interp_ans = trilinear_interpolation( &
                  block(a_blk_no)%xp(tar_x), block(a_blk_no)%yp(tar_y), block(a_blk_no)%zp(tar_z), &
                  loc_x, loc_y, loc_z, &
                  block(b_blk_no)%xp, block(b_blk_no)%yp, block(b_blk_no)%zp, &
                  0, block(b_blk_no)%pc)

                block(a_blk_no)%pc(tar_x, tar_y, tar_z) = bl_interp_ans
                block(a_blk_no)%pco(tar_x, tar_y, tar_z) = bl_interp_ans



        enddo
        enddo
        enddo
        ENDDO
      end subroutine coarseUpdate_pc

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
        integer, intent(in) :: i, j, k
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

        ! Indices used to determine the grid locations after offsets
        integer :: ii, jj, kk
        ! Indices used to determine the variable indices after offsets
        integer :: vi, vj, vk

        ii = i
        jj = j
        kk = k
        vi = i
        vj = j
        vk = k

        ! Set the indices appropiately
        if (offset == 0) then
           ! Leave everything as is
        else if (offset == 1) then
           ii = i - 1
           vi = i - 2
        else if (offset == 2) then
           jj = j - 1
           vj = j - 2
        else if (offset == 3) then
           kk = k - 1
           vk = k - 2
        else
           ! Anything else is an error, but how to handle it?
        end if

        z1 = bilinear_interpolation(x, y, xgrid(ii), xgrid(ii+1), ygrid(jj), ygrid(jj+1), &
             [var(vi, vj, vk), var(vi+1, vj, vk), var(vi, vj+1, vk), var(vi+1, vj+1, vk)])

        z2 = bilinear_interpolation(x, y, xgrid(ii), xgrid(ii+1), ygrid(jj), ygrid(jj+1), &
             [var(vi, vj, vk+1), var(vi+1, vj, vk+1), var(vi, vj+1, vk+1), var(vi+1, vj+1, vk+1)])

        out = linear_interpolation(z, zgrid(kk), zgrid(kk+1), z1, z2)
      end function trilinear_interpolation
end module biocfd_coarse_update

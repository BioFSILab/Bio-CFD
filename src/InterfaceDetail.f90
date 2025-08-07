module biocfd_interface_detail
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, intfr, intflines
  implicit none
  private

  public :: interfaceDetail

  contains
SUBROUTINE interfaceDetail
        INTEGER (int64) :: j, g, factor, a_maxx, a_maxy,a_maxz, a_mm, &
        b_maxx, b_maxy, b_maxz, b_mm, a_max_intf_length, b_max_intf_length, &
        a_blk_no, b_blk_no, xx1_p, xx2_p, &
        a_mm_x, a_mm_y,a_mm_z

        print*, 'before allocation in interface'
        a_max_intf_length=0
        b_max_intf_length=0
        DO g=1, intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk
           factor=intfr(g)%b_msh/intfr(g)%a_msh
           a_maxx=block(a_blk_no)%nx+3
           a_maxy=block(a_blk_no)%ny+3
           a_maxz=block(a_blk_no)%nz+3
           a_mm=max(a_maxx,a_maxy,a_maxz)
           b_maxx=block(b_blk_no)%nx+3
           b_maxy=block(b_blk_no)%ny+3
           b_maxz=block(b_blk_no)%nz+3
           b_mm=max(b_maxx,b_maxy,b_maxz)
        a_mm_x=int(block(b_blk_no)%nx/factor) +3
        a_mm_y=int(block(b_blk_no)%ny/factor) +3
        a_mm_z=int(block(b_blk_no)%nz/factor) +3
        print*,a_mm_x,a_mm_y,a_mm_z
        print*,block(b_blk_no)%nx, factor
        a_mm=max(a_mm,b_mm)
        ALLOCATE(intfr(g)%px_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%ux_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%vx_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%wx_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%pz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%uz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%vz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%wz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%py_interface_det(3,a_mm_y))
        ALLOCATE(intfr(g)%uy_interface_det(3,a_mm_y))
        ALLOCATE(intfr(g)%vy_interface_det(3,a_mm_y))
        ALLOCATE(intfr(g)%wy_interface_det(3,a_mm_y))
        END DO

!coarse mesh start and end index

        print*,'After allocation in interface detail'
        DO g=1, intflines
        factor=intfr(g)%b_msh/intfr(g)%a_msh
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

        ! p
        call set_interface_detail(block(a_blk_no)%xp, intfr(g)%xintf_start, intfr(g)%xintf_end, &
                         factor, intfr(g)%px_interface_det, intfr(g)%counterxp)

        call set_interface_detail(block(a_blk_no)%yp, intfr(g)%yintf_start, intfr(g)%yintf_end, &
                         factor, intfr(g)%py_interface_det, intfr(g)%counteryp)

        call set_interface_detail(block(a_blk_no)%zp, intfr(g)%zintf_start, intfr(g)%zintf_end, &
                         factor, intfr(g)%pz_interface_det, intfr(g)%counterzp)

        ! u
        call set_interface_detail(block(a_blk_no)%xu, intfr(g)%xintf_start, intfr(g)%xintf_end, &
                         factor, intfr(g)%ux_interface_det, intfr(g)%counterxu, 3)

        call set_interface_detail(block(a_blk_no)%yu, intfr(g)%yintf_start, intfr(g)%yintf_end, &
                         factor, intfr(g)%uy_interface_det, intfr(g)%counteryu)

        call set_interface_detail(block(a_blk_no)%zu, intfr(g)%zintf_start, intfr(g)%zintf_end, &
                         factor, intfr(g)%uz_interface_det, intfr(g)%counterzu)

        ! v
        call set_interface_detail(block(a_blk_no)%xv, intfr(g)%xintf_start, intfr(g)%xintf_end, &
                         factor, intfr(g)%vx_interface_det, intfr(g)%counterxv)

        call set_interface_detail(block(a_blk_no)%yv, intfr(g)%yintf_start, intfr(g)%yintf_end, &
                         factor, intfr(g)%vy_interface_det, intfr(g)%counteryv, 3)

        call set_interface_detail(block(a_blk_no)%zv, intfr(g)%zintf_start, intfr(g)%zintf_end, &
                         factor, intfr(g)%vz_interface_det, intfr(g)%counterzv)

        ! w
        call set_interface_detail(block(a_blk_no)%xw, intfr(g)%xintf_start, intfr(g)%xintf_end, &
                         factor, intfr(g)%wx_interface_det, intfr(g)%counterxw)

        call set_interface_detail(block(a_blk_no)%yw, intfr(g)%yintf_start, intfr(g)%yintf_end, &
                         factor, intfr(g)%wy_interface_det, intfr(g)%counteryw)

        call set_interface_detail(block(a_blk_no)%zw, intfr(g)%zintf_start, intfr(g)%zintf_end, &
                         factor, intfr(g)%wz_interface_det, intfr(g)%counterzw, 3)

        ! TN: I expect that the print here has been left over by
        ! mistake, but to get a one-to-one match with the reference
        ! I'm putting this in. TODO: Remove in future if we confirm
        ! not needed.
        xx1_p = intfr(g)%ux_interface_det(1, 1)
        xx2_p = intfr(g)%ux_interface_det(1, intfr(g)%counterxu - 1)
        print*,block(a_blk_no)%xu(xx1_p),block(a_blk_no)%xu(xx2_p),block(b_blk_no)%xu(2)

        END DO

        ! p
        DO j=1,intflines
          call print_interface_detail("px", intfr(j)%counterxp, intfr(j)%px_interface_det, &
                                      block(a_blk_no)%xp, block(b_blk_no)%xp)
        end do
        DO j=1,intfLines
          call print_interface_detail("py", intfr(j)%counteryp, intfr(j)%py_interface_det, &
                                      block(a_blk_no)%yp, block(b_blk_no)%yp)
        end do
        DO j=1,intflines
          call print_interface_detail("pz", intfr(j)%counterzp, intfr(j)%pz_interface_det, &
                                      block(a_blk_no)%zp, block(b_blk_no)%zp)
        end do

        ! u
        DO j=1,intflines
          call print_interface_detail("ux", intfr(j)%counterxu, intfr(j)%ux_interface_det, &
                                      block(a_blk_no)%xu, block(b_blk_no)%xu)
        end do
        DO j=1,intfLines
          call print_interface_detail("uy", intfr(j)%counteryu, intfr(j)%uy_interface_det, &
                                      block(a_blk_no)%yu, block(b_blk_no)%yu)
        end do
        DO j=1,intflines
          call print_interface_detail("uz", intfr(j)%counterzu, intfr(j)%uz_interface_det, &
                                      block(a_blk_no)%zu, block(b_blk_no)%zu)
        end do

        ! v
        DO j=1,intflines
          call print_interface_detail("vx", intfr(j)%counterxv, intfr(j)%vx_interface_det, &
                                      block(a_blk_no)%xv, block(b_blk_no)%xv)
        end do
        DO j=1,intfLines
          call print_interface_detail("vy", intfr(j)%counteryv, intfr(j)%vy_interface_det, &
                                      block(a_blk_no)%yv, block(b_blk_no)%yv)
        end do
        DO j=1,intflines
          call print_interface_detail("vz", intfr(j)%counterzv, intfr(j)%vz_interface_det, &
                                      block(a_blk_no)%zv, block(b_blk_no)%zv)
        end do

        ! w
        DO j=1,intflines
          call print_interface_detail("wx", intfr(j)%counterxw, intfr(j)%wx_interface_det, &
                                      block(a_blk_no)%xw, block(b_blk_no)%xw)
        end do
        DO j=1,intfLines
          call print_interface_detail("wy", intfr(j)%counteryw, intfr(j)%wy_interface_det, &
                                      block(a_blk_no)%yw, block(b_blk_no)%yw)
        end do
        DO j=1,intflines
          call print_interface_detail("wz", intfr(j)%counterzw, intfr(j)%wz_interface_det, &
                                      block(a_blk_no)%zw, block(b_blk_no)%zw)
        end do

        end subroutine interfaceDetail

        !> Set each interface detail array e.g [puvw][xyz]
        subroutine set_interface_detail(block_var, intf_start, intf_end, factor, &
                                        interface_details, interface_counter, in_starter)
          !> The block variable that is being set e.g. xp
          real(dp), dimension(:), intent(in) :: block_var
          !> The interface start and end values
          real(dp), intent(in) :: intf_start, intf_end
          !> The factor i.e. ratio of two meshes
          integer(int64), intent(in) :: factor
          !> The interface details that are actually being set e.g. px_interface_det
          integer(int64), dimension(:, :), intent(inout) :: interface_details
          !> TODO: Check whether the counter really needs to be
          !> returned. Can probably be determined by the size of
          !> interface details. It is weirdly one less than the size
          !> of the interface_details array (dim=2), which seems like
          !> a potential source of bugs.
          integer(int64), intent(inout) :: interface_counter
          ! TN: I don't understand why this works, need to read up on
          ! it. It is not clear to me whether this will work for all
          ! of our compilers, online discussion is not clear about
          ! whether this is really part of the standard
          integer, intent(in), optional, value :: in_starter

          ! TODO: Probably no need for these to be int64 (starter is
          ! either 2 or 3)
          integer(int64) :: starter, ender
          integer(int64) :: i
          integer(int64) :: xx1_p, xx2_p

          ! In most of the examples this subroutine was derived from,
          ! starter starts at 2, however, in some cases it is 3... not
          ! sure why.
          !
          ! TODO: We aren't really sure why we can't directly use
          ! starter but this satisfies the compiler and the linter
          if (.not. present(in_starter)) then
            starter = 2
          else
            starter = in_starter
          end if

          ! The following two do loops can be replaced with findloc
          ! calls
          DO i=1, size(block_var)
            if(block_var(i) > intf_start)then
                    xx1_p=i
                    exit
            endif
          enddo

          DO i=xx1_p, size(block_var)
            if(block_var(i) > intf_end)then
                    xx2_p=i-1
                    exit
            endif
          enddo

          interface_details(1,1) = xx1_p-1
          interface_details(2,1) = starter - 1
          interface_details(3,1) = starter - 1
          ! starter=1+1
          ender=starter+factor-1
          interface_counter = 2

          DO
            interface_details(1, interface_counter) = xx1_p
            interface_details(2, interface_counter) = starter
            interface_details(3, interface_counter) = ender

            starter=ender+1
            ender=ender+factor
            interface_counter = interface_counter + 1
            xx1_p=xx1_p+1
            if(xx1_p > xx2_p) EXIT
          END DO


          interface_details(1, interface_counter) = xx2_p+1
          interface_details(2, interface_counter) = starter
          interface_details(3, interface_counter) = starter
        end subroutine set_interface_detail

        !> Print the details of the interface
        subroutine print_interface_detail(label, counter, interface_details, &
                                          block_var_a, block_var_b)
          !> A two letter label for the details being printed e.g. px
          character(len=2), intent(in) :: label
          !> The counter to iterate through the interface (see TODO in set_interface_detail)
          integer(int64), intent(in) :: counter
          !> The actual interface details being printed
          integer(int64), dimension(:, :), intent(in) :: interface_details
          !> The "block" variables e.g. block(a_blk_no)%zv and block(b_blk_no)%zv
          real(dp), dimension(:), intent(in) :: block_var_a, block_var_b

          integer(int64) :: i
          ! Initial space kept for a perfect match with old code, but
          ! could probably be safely removed
          print "(' **********************',A,'***************************')", label
          DO i=1, counter
          print "(A3,I5,3I5,3F10.5)", label, i, &
          interface_details(1,i), interface_details(2,i), &
          interface_details(3,i), &
          block_var_a(interface_details(1,i)), &
          block_var_b(interface_details(2,i)), &
          block_var_b(interface_details(3,i))
          end do
        end subroutine print_interface_detail

end module biocfd_interface_detail

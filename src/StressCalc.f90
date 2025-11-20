module biocfd_stress_calculation
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use biocfd_forcing, only: compute_value_and_derivatives
  use biocfd_block_type, only: Block_t
  IMPLICIT NONE

  private

  public :: stressCal1

contains

  SUBROUTINE stressCal1(blk, mu_f, rho_f)
    !> The block to perform the stress calculation on (not this is
    !> inout only because of setting thetaDot and thetaDDot, which I
    !> think might not be needed)
    type(Block_t), intent(inout) :: blk
    !> Parameters used in the stress calculation
    real(dp), intent(in) :: mu_f, rho_f

    INTEGER :: ielem, i_x1, i_y1, i_z1

    INTEGER :: i_cell, j_cell, k_cell

    REAL(dp) :: diagdis, normdis, aval, bval, cval, del_X, del_Y, del_Z

    REAL(dp) :: pos1_x, pos1_y, pos1_z, p_pos1, dpdn, dpdn_e, usurf, u_pos1, vsurf, v_pos1, &
               wsurf, w_pos1,  dudn_e, dvdn_e, dwdn_e, ddn_s(3)

    REAL(dp) :: alen, area, area_xz, area_yz, area_xy

    REAL(dp) :: shear_force(3), f_surf(3)

    REAL(dp) :: pressureDrag, viscousDrag, viscousLift, pressureLift, area_Sx, area_Sy, surf_area

    real(dp) :: ac_y, ac_z, at_y, at_z
    real(dp) :: derivatives(3)
    ! In this case I think it makes sense to store the normal - there
    ! is an argument to make cosAlpha, cosBeta, and cosGamma a length
    ! 3 array in Block_t
    real(dp) :: normal(3)

    viscousDrag = 0.
    pressureDrag = 0.
    viscousLift = 0.
    PressureLift = 0.
    surf_area = 0.
    area_Sx = 0.
    area_Sy = 0
    !$acc parallel loop gang vector &
    !$acc private (i_x1, i_y1, i_z1, i_cell, j_cell, k_cell, diagdis, normdis, &
    !$acc          aval, bval, cval, del_X, del_Y, del_Z,                       &
    !$acc          pos1_x, pos1_y, pos1_z,                                   &
    !$acc          p_pos1, dpdn, dpdn_e, &
    !$acc          usurf, u_pos1, vsurf, v_pos1, &
    !$acc          wsurf, w_pos1, &
    !$acc          dudn_e, dvdn_e, dwdn_e, &
    !$acc          ddn_s, &
    !$acc          alen, area, area_xz, area_yz, area_xy, &
    !$acc          shear_force, f_surf, ac_z, ac_y, at_y, at_z)  &
    !$acc default(present)    &
    !$acc firstprivate (rho_f, mu_f) &
    !$acc private(derivatives, normal) &
    !$acc reduction(+: pressureDrag, viscousDrag, viscousLift, pressureLift, area_Sx, area_Sy, surf_area)
    !DO ielem = 1, blk%ibElemCnt
    DO ielem = 1, blk%ibElems
       !  IF((blk%zcent(ielem).ge.0.0).and.(blk%zcent(ielem).le.2.0)) THEN
       !***********************interpolation points****************************

       normal = [blk%cosAlpha(ielem), blk%cosBeta(ielem), blk%cosGamma(ielem)]

       i_cell = find_index_in_array(blk%xcent(ielem), blk%x1, 2_int64, blk%nx+1)
       j_cell = find_index_in_array(blk%ycent(ielem), blk%y1, 2_int64, blk%ny+1)
       k_cell = find_index_in_array(blk%zcent(ielem), blk%z1, 2_int64, blk%nz+1)

       del_X = blk%x1(i_cell+1)-blk%x1(i_cell)
       del_Y = blk%y1(j_cell+1)-blk%y1(j_cell)
       del_Z = blk%z1(k_cell+1)-blk%z1(k_cell)

       diagdis = sqrt(del_X**2 + del_Y**2 + del_Z**2)

       normdis = diagdis
       pos1_x = blk%xcent(ielem) + normdis * normal(1)
       pos1_y = blk%ycent(ielem) + normdis * normal(2)
       pos1_z = blk%zcent(ielem) + normdis * normal(3)

       !**************************velocity and pressure at the surface**********************
       IF (blk%ibSurfID(ielem)==50) THEN
          blk% thetaDot  =  0.
          blk% thetaDDot =  0.
       ELSE IF (blk%ibSurfID(ielem)==51) THEN
          blk% thetaDot  = blk% thetaDot1
          blk% thetaDDot = blk% thetaDDot1
       ELSE IF (blk%ibSurfId(ielem)==52) THEN
          blk% thetaDot  = blk% thetaDot2
          blk% thetaDDot = blk% thetaDDot2
       END IF

       usurf = blk%xdot
       vsurf = -blk%thetaDot * (blk%zcent(ielem) - blk%piv_z) + blk%ydot
       wsurf = blk%thetaDot * (blk%ycent(ielem) - blk%piv_y)

       ac_z = -blk%thetaDot**2 * (blk%zcent(ielem) - blk%piv_z)
       ac_y = -blk%thetaDot**2 * (blk%ycent(ielem) - blk%piv_y)
       at_z = blk%thetaDDot * (blk%ycent(ielem) - blk%piv_y)
       at_y = -blk%thetaDDot * (blk%zcent(ielem) - blk%piv_z)

       dpdn = -((ac_z + at_z)*blk%cosGamma(ielem) &
            + (ac_y + at_y)*blk%cosBeta(ielem))-blk%yddot*blk%cosBeta(ielem)

       !*******************velocity interpolation at point 2******************

       !******************u velocity interpolation at point 2******************

       i_x1 = find_index_in_array(pos1_x, blk%xu, 2_int64, blk%nx+1)
       i_y1 = find_index_in_array(pos1_y, blk%yu, 2_int64, blk%ny+1)
       i_z1 = find_index_in_array(pos1_z, blk%zu, 2_int64, blk%nz+1)

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
            blk%xu, blk%yu, blk%zu, 1, &
            blk%u, u_pos1, derivatives)

       dudn_e = dot_product(derivatives, normal)
       ddn_s(1) = (2._dp/normdis)*(u_pos1 - usurf) - dudn_e

       !******************v velocity interpolation in point 2******************
       i_x1 = find_index_in_array(pos1_x, blk%xv, 2_int64, blk%nx+1)
       i_y1 = find_index_in_array(pos1_y, blk%yv, 2_int64, blk%ny+1)
       i_z1 = find_index_in_array(pos1_z, blk%zv, 2_int64, blk%nz+1)

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
            blk%xv, blk%yv, blk%zv, 2, &
            blk%v, v_pos1, derivatives)

       dvdn_e = dot_product(derivatives, normal)
       ddn_s(2) = (2._dp/normdis)*(v_pos1 - vsurf) - dvdn_e

       !******************w velocity interpolation in point 2******************
       i_x1 = find_index_in_array(pos1_x, blk%xw, 2_int64, blk%nx+1)
       i_y1 = find_index_in_array(pos1_y, blk%yw, 2_int64, blk%ny+1)
       i_z1 = find_index_in_array(pos1_z, blk%zw, 2_int64, blk%nz+1)

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
            blk%xw, blk%yw, blk%zw, 3, &
            blk%w, w_pos1, derivatives)

       dwdn_e = dot_product(derivatives, normal)
       ddn_s(3) = (2._dp/normdis)*(w_pos1 - wsurf) - dwdn_e

       !***********************calculate area of the elements******************
       alen = blk%element_length(ielem)
       area = alen/2._dp
       area_yz = 0.5_dp*abs(alen * normal(1))
       area_xz = 0.5_dp*abs(alen * normal(2))
       area_xy = 0.5_dp*abs(alen * normal(3))

       !*********non-dimensional viscous stress & force calculation************
       shear_force = (ddn_s - dot_product(ddn_s, normal) * normal) * mu_f * area
       !************************presssure interpolation************************

       !*******************pressure interpolation at point 2*******************
       i_x1 = find_index_in_array(pos1_x, blk%xp, 2_int64, blk%nx+1)
       i_y1 = find_index_in_array(pos1_y, blk%yp, 2_int64, blk%ny+1)
       i_z1 = find_index_in_array(pos1_z, blk%zp, 2_int64, blk%nz+1)

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
            blk%xp, blk%yp, blk%zp, 0, &
            blk%p, p_pos1, derivatives)

       dpdn_e = dot_product(derivatives, normal)

       bval = dpdn  !dpdn=-dudt
       aval = (dpdn_e - dpdn)/(2*diagdis)
       cval = p_pos1 - (dpdn_e + dpdn)*diagdis*0.5_dp

       f_surf = -cval * area * normal * rho_f

       !***********************drag calculation********************************
       viscousDrag = viscousDrag + shear_force(1)
       pressureDrag = pressureDrag + f_surf(1)
       viscousLift = viscousLift + shear_force(2)
       PressureLift = PressureLift + f_surf(2)
       surf_area = surf_area + area
       area_Sx = area_Sx + area_xz
       area_Sy = area_Sy + area_yz

    END DO
    !$acc end parallel loop

    area_Sx= 0.5_dp * area_Sx
    area_Sy= 0.5_dp * area_Sy

    blk%viscous_drag_coefficient = 2 * (viscousDrag/area_Sx)
    blk%pressure_drag_coefficient = 2 * (pressureDrag/area_Sy)
    blk%viscous_lift_coefficient = 2 * (viscousLift/area_Sx)
    blk%pressure_lift_coefficient = 2 * (pressureLift/area_Sy)

  end subroutine stressCal1

!> Find the position in the array where the value is greater than
!> element i but less than element i+1
pure function find_index_in_array(value, array, start, end) result(index)
   real(dp), intent(in) :: value
   real(dp), intent(in) :: array(:)
   ! TODO: No need for these to be int64
   integer(int64), intent(in) :: start, end

   !> The resulting index
   integer :: index

   ! Internal counter
   integer :: i
   !$acc routine seq
   do i=start, end
      if (value >= array(i) .and. value < array(i+1)) then
         index = i
         return  ! As soon as we find a value we can return
      end if
   end do
end function find_index_in_array

end module biocfd_stress_calculation

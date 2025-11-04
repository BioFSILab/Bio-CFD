module biocfd_stress_calculation
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  USE global, only: block, blk_start, mu_f, rho_f, totime
  use biocfd_forcing, only: compute_value_and_derivatives
  use biocfd_block_type, only: Block_t
  IMPLICIT NONE

  private

  public :: stressCal1

contains

  SUBROUTINE stressCal1(blk, id)

    type(Block_t), intent(inout) :: blk
    integer(int64), intent(in) :: id

    INTEGER:: i, j, k, ielem, i_x1, i_y1, i_z1

    INTEGER:: i_cell, j_cell, k_cell

    REAL(dp):: diagdis, normdis, aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z

    REAL(dp):: pos1_x, pos1_y, pos1_z, psurf, p_pos1, dpdn, dpdn_e, usurf, u_pos1, vsurf, v_pos1, &
               wsurf, w_pos1,  dudn_e, dvdn_e, dwdn_e, dudn_s, dvdn_s, dwdn_s

    REAL(dp):: alen, area, area_xz, area_yz, area_xy

    REAL(dp):: shear_x_force, shear_y_force, shear_z_force, f_surf, f_surf_x, f_surf_y, f_surf_z

    REAL(dp):: pressureDrag, viscousDrag, viscousLift, pressureLift, viscousDragcoefficient, &
               pressureDragcoefficient, viscousLiftcoefficient, PressureLiftcoefficient, &
               area_Sx, area_Sy, surf_area

    real(dp) :: ac_y, ac_z, at_y, at_z
    real(dp) :: derivatives(3)

    CHARACTER(len=150) :: filename1

    viscousDrag = 0.
    pressureDrag = 0.
    viscousLift = 0.
    PressureLift = 0.
    surf_area = 0.
    area_Sx = 0.
    area_Sy = 0
    !$acc parallel loop gang vector reduction(+: pressureDrag, viscousDrag, viscousLift, pressureLift, area_Sx, area_Sy, surf_area)  &
    !$acc private (i_x1, i_y1, i_z1, i_cell, j_cell, k_cell, diagdis, normdis,                    &
    !$acc          aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z,                       &
    !$acc          pos1_x, pos1_y, pos1_z,                                   &
    !$acc          psurf, p_pos1, dpdn, dpdn_e, &
    !$acc          usurf, u_pos1, vsurf, v_pos1, &
    !$acc          wsurf, w_pos1, &
    !$acc          dudn_e, dvdn_e, dwdn_e, &
    !$acc          dudn_s, dvdn_s, dwdn_s, &
    !$acc          alen, area, area_xz, area_yz, area_xy, shear_x_force, shear_y_force,           &
    !$acc          shear_z_force, f_surf, f_surf_x, f_surf_y, f_surf_z,ac_z,ac_y,ac_x,at_y,at_z)  &
    !$acc default(present)    &
    !$acc firstprivate (blk%nx, blk%ny, blk%nz, deltat, rho_f, re, mu_f)
    !$acc private(derivatives)
    !DO ielem = 1, blk%ibElemCnt
    DO ielem = 1, blk%ibElems
       !  IF((blk%zcent(ielem).ge.0.0).and.(blk%zcent(ielem).le.2.0)) THEN
       !***********************interpolation points****************************

       i_cell = find_index_in_array(blk%xcent(ielem), blk%x1, 2_int64, blk%nx+1)
       j_cell = find_index_in_array(blk%ycent(ielem), blk%y1, 2_int64, blk%ny+1)
       k_cell = find_index_in_array(blk%zcent(ielem), blk%z1, 2_int64, blk%nz+1)

       del_X = blk%x1(i_cell+1)-blk%x1(i_cell)
       del_Y = blk%y1(j_cell+1)-blk%y1(j_cell)
       del_Z = blk%z1(k_cell+1)-blk%z1(k_cell)

       diagdis = dsqrt(del_X**2 + del_Y**2 + del_Z**2)

       normdis = diagdis
       pos1_x = blk%xcent(ielem) + normdis*blk%cosAlpha(ielem)
       pos1_y = blk%ycent(ielem) + normdis*blk%cosBeta(ielem)
       pos1_z = blk%zcent(ielem) + normdis*blk%cosGamma(ielem)

       !**************************velocity and pressure at the surface**********************
       IF (blk%ibSurfID(ielem)==50) THEN
          blk% thetaDot  =  0.
          blk% thetaDDot =  0.
          usurf = 0._dp + blk%xdot
          vsurf = 0._dp + blk%ydot
          wsurf = 0._dp
          ac_z      =  0._dp  !-thetaDot**2*(zcent(nelp(index_ts(n))) - piv_z)
          ac_y      =  0._dp  !-thetaDot**2*(ycent(nelp(index_ts(n))) - piv_y)
          at_z      =  0._dp  ! thetaDDot*(ycent(nelp(index_ts(n))) - piv_y)
          at_y      =  0._dp  !
       ELSEIF (blk%ibSurfID(ielem)==51) THEN
          blk% thetaDot  = blk% thetaDot1
          blk% thetaDDot = blk% thetaDDot1
          usurf    = 0._dp + blk%xdot
          vsurf    = -blk%thetaDot*(blk%zcent(ielem) - blk%piv_z) + blk%ydot
          wsurf    = blk%thetaDot*(blk%ycent(ielem) - blk%piv_y)
          ac_z      = -blk%thetaDot**2*(blk%zcent(ielem) -blk% piv_z)
          ac_y      = -blk%thetaDot**2*(blk%ycent(ielem) -blk% piv_y)
          at_z      = blk% thetaDDot*(blk%ycent(ielem)-blk% piv_y)
          at_y      = -blk%thetaDDot*(blk%zcent(ielem)-blk% piv_z)
       ELSEIF (blk%ibSurfId(ielem)==52) THEN
          blk% thetaDot  = blk% thetaDot2
          blk% thetaDDot = blk% thetaDDot2
          usurf = 0._dp +blk%xdot
          vsurf    = -blk%thetaDot*(blk%zcent(ielem) - blk%piv_z)+ blk%ydot  ! + ydot
          wsurf    = blk%thetaDot*(blk%ycent(ielem) - blk%piv_y)  ! + ydot
          ac_z      = -blk%thetaDot**2*(blk%zcent(ielem) -blk% piv_z)
          ac_y      = -blk%thetaDot**2*(blk%ycent(ielem) - blk%piv_y)
          at_z      =  blk%thetaDDot*(blk%ycent(ielem) -blk% piv_y)
          at_y      = -blk%thetaDDot*(blk%zcent(ielem) - blk%piv_z)
       ENDIF
       !       usurf     =  0.
       !       wsurf     = blk% thetaDot*(blk%ycent(ielem) -blk% piv_y)
       !       vsurf     = -blk%thetaDot*(blk%zcent(ielem) -blk% piv_z)
       !       dpdn      = -((ac_z + at_z)*blk%cosAlpha(ielem)  + (ac_y + at_y)*blk%cosBeta(ielem))
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

       dudn_e = derivatives(1) * blk%cosAlpha(ielem) &
            + derivatives(2) * blk%cosBeta(ielem) &
            + derivatives(3) * blk%cosGamma(ielem)

       dudn_s = (2._dp/normdis)*(u_pos1 - usurf) - dudn_e

       !******************v velocity interpolation in point 2******************
       i_x1 = find_index_in_array(pos1_x, blk%xv, 2_int64, blk%nx+1)
       i_y1 = find_index_in_array(pos1_y, blk%yv, 2_int64, blk%ny+1)
       i_z1 = find_index_in_array(pos1_z, blk%zv, 2_int64, blk%nz+1)

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
            blk%xv, blk%yv, blk%zv, 2, &
            blk%v, v_pos1, derivatives)

       dvdn_e = derivatives(1) * blk%cosAlpha(ielem) &
            + derivatives(2) * blk%cosBeta(ielem) &
            + derivatives(3) * blk%cosGamma(ielem)

       dvdn_s = (2._dp/normdis)*(v_pos1 - vsurf) - dvdn_e

       !******************w velocity interpolation in point 2******************
       i_x1 = find_index_in_array(pos1_x, blk%xw, 2_int64, blk%nx+1)
       i_y1 = find_index_in_array(pos1_y, blk%yw, 2_int64, blk%ny+1)
       i_z1 = find_index_in_array(pos1_z, blk%zw, 2_int64, blk%nz+1)

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
            blk%xw, blk%yw, blk%zw, 3, &
            blk%w, w_pos1, derivatives)

       dwdn_e = derivatives(1) * blk%cosAlpha(ielem) &
            + derivatives(2) * blk%cosBeta(ielem) &
            + derivatives(3) * blk%cosGamma(ielem)

       dwdn_s = (2._dp/normdis)*(w_pos1 - wsurf) - dwdn_e

       !***********************calculate area of the elements******************
       ! alen = sqrt(blk%alpha3(ielem)**2.+blk%beta3(ielem)**2.+ blk%gamma3(ielem)**2.)
       alen = blk%element_length(ielem)
       area = alen/2._dp
       area_yz = 0.5_dp*abs(alen * blk%cosAlpha(ielem))
       area_xz = 0.5_dp*abs(alen * blk%cosBeta(ielem))
       area_xy = 0.5_dp*abs(alen * blk%cosGamma(ielem))

       !*********non-dimensional viscous stress & force calculation************
       stx1 = dudn_s - (dudn_s*blk%cosAlpha(ielem) + dvdn_s*blk%cosBeta(ielem) &
            + dwdn_s*blk%cosGamma(ielem))*blk%cosAlpha(ielem)

       sty1 = dvdn_s - (dudn_s*blk%cosAlpha(ielem) + dvdn_s*blk%cosBeta(ielem) &
            + dwdn_s*blk%cosGamma(ielem))*blk%cosBeta(ielem)

       stz1 = dwdn_s - (dudn_s*blk%cosAlpha(ielem) + dvdn_s*blk%cosBeta(ielem) &
            + dwdn_s*blk%cosGamma(ielem))*blk%cosGamma(ielem)

       shear_x_force = mu_f*stx1*area
       shear_y_force = mu_f*sty1*area
       shear_z_force = mu_f*stz1*area

       !************************presssure interpolation************************

       !*******************pressure interpolation at point 2*******************
       i_x1 = find_index_in_array(pos1_x, blk%xp, 2_int64, blk%nx+1)
       i_y1 = find_index_in_array(pos1_y, blk%yp, 2_int64, blk%ny+1)
       i_z1 = find_index_in_array(pos1_z, blk%zp, 2_int64, blk%nz+1)

       call compute_value_and_derivatives(pos1_x, pos1_y, pos1_z, i_x1, i_y1, i_z1, &
            blk%xp, blk%yp, blk%zp, 0, &
            blk%p, p_pos1, derivatives)

       dpdn_e = derivatives(1)*blk%cosAlpha(ielem) &
            + derivatives(2)*blk%cosBeta(ielem) &
            + derivatives(3)*blk%cosGamma(ielem)

       bval = dpdn  !dpdn=-dudt
       aval = (dpdn_e - dpdn)/(2*diagdis)
       cval = p_pos1 - (dpdn_e + dpdn)*diagdis*0.5_dp

       psurf = cval

       f_surf = cval*area

       f_surf_x = -f_surf*blk%cosAlpha(ielem)*rho_f

       f_surf_y = -f_surf*blk%cosBeta(ielem)*rho_f

       f_surf_z = -f_surf*blk%cosGamma(ielem)*rho_f

       !***********************drag calculation********************************
       viscousDrag = viscousDrag + shear_x_force
       pressureDrag = pressureDrag + f_surf_x
       viscousLift = viscousLift + shear_y_force
       PressureLift = PressureLift + f_surf_y
       surf_area = surf_area + area
       area_Sx = area_Sx + area_xz
       area_Sy = area_Sy + area_yz

    END DO
    !$acc end parallel

    area_Sx= 0.5_dp * area_Sx
    area_Sy= 0.5_dp * area_Sy
    ! area_Sx= 1.!0.5 * area_Sx
    ! area_Sy= 1.!0.5 * area_Sy
    viscousDragcoefficient= 2* (viscousDrag/area_Sx)
    pressureDragcoefficient= 2* (pressureDrag/area_Sy)
    viscousLiftcoefficient= 2* (viscousLift/area_Sx)
    pressureLiftcoefficient= 2* (pressureLift/area_Sy)

    !*********************drag file writing*********************************

    WRITE(filename1,19) id
19  FORMAT('dragcoff_',I4.4,'.dat')
    OPEN(899,file=filename1,Access='Append',status='unknown')
    WRITE(899,*) viscousDragcoefficient, pressureDragcoefficient, totime
    CLOSE(899)
    WRITE(filename1,29) id
29  FORMAT('liftcoff_',I4.4,'.dat')
    OPEN(999,file=filename1,Access='Append',status='unknown')
    WRITE(999,*) viscousLiftcoefficient,  pressureLiftcoefficient, totime
    CLOSE(999)

  end subroutine stressCal1

!> Find the position in the array where the value is greater than
!> element i but less than element i-1
pure function find_index_in_array(value, array, start, end) result(index)
   real(dp), intent(in) :: value
   real(dp), intent(in) :: array(:)
   ! TODO: No need for these to be int64
   integer(int64), intent(in) :: start, end

   !> The resulting index
   integer :: index

   ! Internal counter
   integer :: i
   !$acc loop seq
   do i=start, end
      if (value >= array(i) .and. value < array(i+1)) then
         index = i
         return  ! As soon as we find a value we can return
      end if
   end do
end function find_index_in_array

end module biocfd_stress_calculation

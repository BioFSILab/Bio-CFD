module biocfd_navier_stokes
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use global, only : block, al, alpha, deltat, nblocks, re, rev
  implicit none
  private

  public :: non_uni_coeff, nsmomentum2order

contains

       SUBROUTINE non_uni_coeff
       INTEGER  (dp) :: i, j, k, g, nx_var, ny_var, nz_var
       REAL (dp)   :: tmp_dx1, tmp_dx2, tmp_dx3, tmp_dx4, tmp_dy1, tmp_dy2, tmp_dy3,         &
                      tmp_dy4, tmp_dz1, tmp_dz2, tmp_dz3, tmp_dz4

       real(dp) :: theta(3)
       real(dp) :: f(8)
       !> For s we will use special indexing, note that original code
       !> doesn't have s77 or 87, but we'll just live with that for
       !> now
       real(dp) :: s(5:10, 1:7)
       real(dp) :: ak(7)


        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz

        allocate(block(g)%ca_uu(6, nx_var+2), block(g)%ck_uu(6, nx_var+2))
        allocate(block(g)%ca_vv(6, ny_var+2), block(g)%ck_vv(6, ny_var+2))
        allocate(block(g)%ca_ww(6, nz_var+2), block(g)%ck_ww(6, nz_var+2))
        allocate(block(g)%ca_uv(6, nx_var+2), block(g)%ck_uv(6, nx_var+2))
        allocate(block(g)%ca_uw(6, nx_var+2), block(g)%ck_uw(6, nx_var+2))
        allocate(block(g)%ca_vu(6, ny_var+2), block(g)%ck_vu(6, ny_var+2))
        allocate(block(g)%ca_vw(6, ny_var+2), block(g)%ck_vw(6, ny_var+2))
        allocate(block(g)%ca_wu(6, nz_var+2), block(g)%ck_wu(6, nz_var+2))
        allocate(block(g)%ca_wv(6, nz_var+2), block(g)%ck_wv(6, nz_var+2))

       do i=2, block(g)%nx
       theta = block(g)%deltax(i:i+2) / block(g)%deltax(i-1:i+1)
       f = compute_f(theta)
       s = compute_s(theta(2), f)

       block(g)%ca_uu(1, i) = s(9, 7) * s(10, 1) - s(10, 7) * s(9, 1)
       block(g)%ca_uu(2, i) = s(9, 7) * s(10, 2) + s(10, 7) * s(9, 3)
       block(g)%ca_uu(3, i) = s(9, 5) * s(10, 7) - s(10, 5) * s(9, 7)
       block(g)%ca_uu(4, i) = s(9, 4) * s(10, 7) + s(10, 4) * s(9, 7)
       block(g)%ca_uu(5, i) = s(9, 7) * s(10, 3) - s(10, 7) * s(9, 2)
       block(g)%ca_uu(6, i) = s(9, 7) * s(10, 6) - s(10, 7) * s(9, 6)

       ak = compute_ak(theta)

       block(g)%ck_uu(1, i) = ak(3)
       block(g)%ck_uu(2, i) = -ak(4)
       block(g)%ck_uu(3, i) = ak(1) + ak(2)
       block(g)%ck_uu(4, i) = -ak(5)
       block(g)%ck_uu(5, i) = ak(6)
       block(g)%ck_uu(6, i) = ak(7)
       enddo

       do j=2, block(g)%ny

       theta = block(g)%deltay(j:j+2) / block(g)%deltay(j-1:j+1)

       f = compute_f(theta)
       s = compute_s(theta(2), f)

       block(g)%ca_vv(1, j) = s(9, 7) * s(10, 1) - s(10, 7) * s(9, 1)
       block(g)%ca_vv(2, j) = s(9, 7) * s(10, 2) + s(10, 7) * s(9, 3)
       block(g)%ca_vv(3, j) = s(9, 5) * s(10, 7) - s(10, 5) * s(9, 7)
       block(g)%ca_vv(4, j) = s(9, 4) * s(10, 7) + s(10, 4) * s(9, 7)
       block(g)%ca_vv(5, j) = s(9, 7) * s(10, 3) - s(10, 7) * s(9, 2)
       block(g)%ca_vv(6, j) = s(9, 7) * s(10, 6) - s(10, 7) * s(9, 6)

       ak = compute_ak(theta)

       block(g)%ck_vv(1, j) = ak(3)
       block(g)%ck_vv(2, j) = -ak(4)
       block(g)%ck_vv(3, j) = ak(1) + ak(2)
       block(g)%ck_vv(4, j) = -ak(5)
       block(g)%ck_vv(5, j) = ak(6)
       block(g)%ck_vv(6, j) = ak(7)
       enddo

       do k=2, block(g)%nz

       theta = block(g)%deltaz(k:k+2) / block(g)%deltaz(k-1:k+1)

       f = compute_f(theta)
       s = compute_s(theta(2), f)

       block(g)%ca_ww(1, k) = s(9, 7) * s(10, 1) - s(10, 7) * s(9, 1)
       block(g)%ca_ww(2, k) = s(9, 7) * s(10, 2) + s(10, 7) * s(9, 3)
       block(g)%ca_ww(3, k) = s(9, 5) * s(10, 7) - s(10, 5) * s(9, 7)
       block(g)%ca_ww(4, k) = s(9, 4) * s(10, 7) + s(10, 4) * s(9, 7)
       block(g)%ca_ww(5, k) = s(9, 7) * s(10, 3) - s(10, 7) * s(9, 2)
       block(g)%ca_ww(6, k) = s(9, 7) * s(10, 6) - s(10, 7) * s(9, 6)

       ak = compute_ak(theta)

       block(g)%ck_ww(1, k) = ak(3)
       block(g)%ck_ww(2, k) = -ak(4)
       block(g)%ck_ww(3, k) = ak(1) + ak(2)
       block(g)%ck_ww(4, k) = -ak(5)
       block(g)%ck_ww(5, k) = ak(6)
       block(g)%ck_ww(6, k) = ak(7)
       enddo

       do i=2, block(g)%nx
       if(i==2)then
       tmp_dx1=block(g)%deltax(i-1)
       else
       tmp_dx1=0.5_dp*(block(g)%deltax(i-1)+block(g)%deltax(i-2))
       endif
       tmp_dx2=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i-1))
       tmp_dx3=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i+1))
       tmp_dx4=0.5_dp*(block(g)%deltax(i+1)+block(g)%deltax(i+2))

       theta = [tmp_dx2/tmp_dx1, tmp_dx3/tmp_dx2, tmp_dx4/tmp_dx3]

       f = compute_f(theta)
       s = compute_s(theta(2), f)

       block(g)%ca_uv(1, i) = s(9, 7) * s(10, 1) - s(10, 7) * s(9, 1)
       block(g)%ca_uv(2, i) = s(9, 7) * s(10, 2) + s(10, 7) * s(9, 3)
       block(g)%ca_uv(3, i) = s(9, 5) * s(10, 7) - s(10, 5) * s(9, 7)
       block(g)%ca_uv(4, i) = s(9, 4) * s(10, 7) + s(10, 4) * s(9, 7)
       block(g)%ca_uv(5, i) = s(9, 7) * s(10, 3) - s(10, 7) * s(9, 2)
       block(g)%ca_uv(6, i) = s(9, 7) * s(10, 6) - s(10, 7) * s(9, 6)

       block(g)%ca_uw(1, i) = block(g)%ca_uv(1, i)
       block(g)%ca_uw(2, i) = block(g)%ca_uv(2, i)
       block(g)%ca_uw(3, i) = block(g)%ca_uv(3, i)
       block(g)%ca_uw(4, i) = block(g)%ca_uv(4, i)
       block(g)%ca_uw(5, i) = block(g)%ca_uv(5, i)
       block(g)%ca_uw(6, i) = block(g)%ca_uv(6, i)

       ak = compute_ak(theta)

       block(g)%ck_uv(1, i) = ak(3)
       block(g)%ck_uv(2, i) = -ak(4)
       block(g)%ck_uv(3, i) = ak(1)+ak(2)
       block(g)%ck_uv(4, i) = -ak(5)
       block(g)%ck_uv(5, i) = ak(6)
       block(g)%ck_uv(6, i) = ak(7)

       block(g)%ck_uw(1, i)=block(g)%ck_uv(1, i)
       block(g)%ck_uw(2, i)=block(g)%ck_uv(2, i)
       block(g)%ck_uw(3, i)=block(g)%ck_uv(3, i)
       block(g)%ck_uw(4, i)=block(g)%ck_uv(4, i)
       block(g)%ck_uw(5, i)=block(g)%ck_uv(5, i)
       block(g)%ck_uw(6, i)=block(g)%ck_uv(6, i)
       enddo

       do j=2, block(g)%ny
       if(j==2)then
       tmp_dy1=block(g)%deltay(j-1)
       else
       tmp_dy1=0.5_dp*(block(g)%deltay(j-1)+block(g)%deltay(j-2))
       endif
       tmp_dy2=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j-1))
       tmp_dy3=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j+1))
       tmp_dy4=0.5_dp*(block(g)%deltay(j+1)+block(g)%deltay(j+2))

       theta= [tmp_dy2/tmp_dy1, tmp_dy3/tmp_dy2, tmp_dy4/tmp_dy3]

       f = compute_f(theta)
       s = compute_s(theta(2), f)

       block(g)%ca_vu(1, j) = s(9, 7) * s(10, 1) - s(10, 7) * s(9, 1)
       block(g)%ca_vu(2, j) = s(9, 7) * s(10, 2) + s(10, 7) * s(9, 3)
       block(g)%ca_vu(3, j) = s(9, 5) * s(10, 7) - s(10, 5) * s(9, 7)
       block(g)%ca_vu(4, j) = s(9, 4) * s(10, 7) + s(10, 4) * s(9, 7)
       block(g)%ca_vu(5, j) = s(9, 7) * s(10, 3) - s(10, 7) * s(9, 2)
       block(g)%ca_vu(6, j) = s(9, 7) * s(10, 6) - s(10, 7) * s(9, 6)

       block(g)%ca_vw(1, j)=block(g)%ca_vu(1, j)
       block(g)%ca_vw(2, j)=block(g)%ca_vu(2, j)
       block(g)%ca_vw(3, j)=block(g)%ca_vu(3, j)
       block(g)%ca_vw(4, j)=block(g)%ca_vu(4, j)
       block(g)%ca_vw(5, j)=block(g)%ca_vu(5, j)
       block(g)%ca_vw(6, j)=block(g)%ca_vu(6, j)

       ak = compute_ak(theta)

       block(g)%ck_vu(1, j) = ak(3)
       block(g)%ck_vu(2, j) = -ak(4)
       block(g)%ck_vu(3, j) = ak(1) + ak(2)
       block(g)%ck_vu(4, j) = -ak(5)
       block(g)%ck_vu(5, j) = ak(6)
       block(g)%ck_vu(6, j) = ak(7)

       block(g)%ck_vw(1, j)=block(g)%ck_vu(1, j)
       block(g)%ck_vw(2, j)=block(g)%ck_vu(2, j)
       block(g)%ck_vw(3, j)=block(g)%ck_vu(3, j)
       block(g)%ck_vw(4, j)=block(g)%ck_vu(4, j)
       block(g)%ck_vw(5, j)=block(g)%ck_vu(5, j)
       block(g)%ck_vw(6, j)=block(g)%ck_vu(6, j)
       enddo

       do k=2, block(g)%nz
       if(k==2)then
       tmp_dz1=block(g)%deltaz(k-1)
       else
       tmp_dz1=0.5_dp*(block(g)%deltaz(k-1)+block(g)%deltaz(k-2))
       endif
       tmp_dz2=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
       tmp_dz3=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k+1))
       tmp_dz4=0.5_dp*(block(g)%deltaz(k+1)+block(g)%deltaz(k+2))

       theta = [tmp_dz2/tmp_dz1, tmp_dz3/tmp_dz2, tmp_dz4/tmp_dz3]

       f = compute_f(theta)
       s = compute_s(theta(2), f)

       block(g)%ca_wu(1, k) = s(9, 7) * s(10, 1) - s(10, 7) * s(9, 1)
       block(g)%ca_wu(2, k) = s(9, 7) * s(10, 2) + s(10, 7) * s(9, 3)
       block(g)%ca_wu(3, k) = s(9, 5) * s(10, 7) - s(10, 5) * s(9, 7)
       block(g)%ca_wu(4, k) = s(9, 4) * s(10, 7) + s(10, 4) * s(9, 7)
       block(g)%ca_wu(5, k) = s(9, 7) * s(10, 3) - s(10, 7) * s(9, 2)
       block(g)%ca_wu(6, k) = s(9, 7) * s(10, 6) - s(10, 7) * s(9, 6)

       block(g)%ca_wv(1, k)=block(g)%ca_wu(1, k)
       block(g)%ca_wv(2, k)=block(g)%ca_wu(2, k)
       block(g)%ca_wv(3, k)=block(g)%ca_wu(3, k)
       block(g)%ca_wv(4, k)=block(g)%ca_wu(4, k)
       block(g)%ca_wv(5, k)=block(g)%ca_wu(5, k)
       block(g)%ca_wv(6, k)=block(g)%ca_wu(6, k)

       ak = compute_ak(theta)

       block(g)%ck_wu(1, k) = ak(3)
       block(g)%ck_wu(2, k) = -ak(4)
       block(g)%ck_wu(3, k) = ak(1) + ak(2)
       block(g)%ck_wu(4, k) = -ak(5)
       block(g)%ck_wu(5, k) = ak(6)
       block(g)%ck_wu(6, k) = ak(7)

       block(g)%ck_wv(1, k)=block(g)%ck_wu(1, k)
       block(g)%ck_wv(2, k)=block(g)%ck_wu(2, k)
       block(g)%ck_wv(3, k)=block(g)%ck_wu(3, k)
       block(g)%ck_wv(4, k)=block(g)%ck_wu(4, k)
       block(g)%ck_wv(5, k)=block(g)%ck_wu(5, k)
       block(g)%ck_wv(6, k)=block(g)%ck_wu(6, k)
       enddo
        END DO

       write(*,*)'leaving non_uni_coeff'

       END SUBROUTINE non_uni_coeff

      subroutine nsMomentum2order
!c***********************************************************************
!c     navier-stokes equations for constant properties
!c***********************************************************************
      INTEGER (dp) :: i, j, k, g, n ,nx_var,ny_var,nz_var, n1, nn, i11, j11, k11, &
           index_ip1, index_im1, index_jp1, index_jm1, index_kp1, index_km1
         REAL (dp) :: dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14, &
                 u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16, &
                       w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr, &
                       dx2xl,dy2ye,dy2yw, dxr, dx, dxl, dye, dy, dyw, dzt, dz, dzb, &
                       dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,dwudz, &
                       wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy, &
                       dwutdz,duudx,dvudy,d2udx2,d2udy2,d2udz2, &
                       uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz, &
                       duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w, &
                       u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz, &
                       duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,&
                          ztt2,residw, temp_u1dotn, temp_u2dotn, temp_v1dotn, temp_v2dotn, &
                          temp_w1dotn, temp_w2dotn, temp_pdotn
     al = 1.

        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz

!$acc parallel loop gang vector  &
!$acc private (i, j, k, n1, dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,   &
!$acc	        u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
!$acc		 w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
!$acc	        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
!$acc		 dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,          &
!$acc		 wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
!$acc	        dwutdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
!$acc	        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,            &
!$acc	        duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
!$acc	        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,        &
!$acc		 duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw, &
!$acc           index_ip1,index_im1,index_jp1,index_jm1,index_kp1,index_km1, &
!$acc           i11,j11,k11,temp_u2dotn,temp_u1dotn, &
!$acc           temp_v2dotn,temp_v1dotn, &
!$acc           temp_w2dotn,temp_w1dotn, temp_pdotn) &
!$acc default(present)   &
!$acc firstprivate(nx_var, ny_var, nz_var, rev, deltat, al)
     DO n = 1, block(g)%fluidCellCount

      i = block(g)%fluidIndexPtr(n, 1)
      j = block(g)%fluidIndexPtr(n, 2)
      k = block(g)%fluidIndexPtr(n, 3)
      index_ip1 = 0
      index_im1 = 0
      index_jp1 = 0
      index_jm1 = 0
      index_kp1 = 0
      index_km1 = 0
      n1 = i-1 + nx_var*(j-2) + nx_var*ny_var*(k-2)
      IF (block(g)%cell2(i+1, j, k)==2 .OR. block(g)%cell2(i-1, j, k)==2 .OR. &
           block(g)%cell2(i, j+1, k)==2 .OR. block(g)%cell2(i, j-1, k)==2 .OR. &
           block(g)%cell2(i, j, k+1)==2.OR. block(g)%cell2(i, j, k-1)==2) THEN
        !$acc loop seq
             DO nn = 1, block(g)%TSCellCount
                i11 = block(g)%TSIndexPtr(nn, 1)
                j11 = block(g)%TSIndexPtr(nn, 2)
                k11 = block(g)%TSIndexPtr(nn, 3)

                IF (block(g)%cell2(i+1, j, k)==2) THEN
                   index_ip1 = nn
                   temp_u2dotn = 0.
                   temp_u1dotn = 0.
                   temp_v2dotn = 0.
                   temp_v1dotn = 0.
                   temp_w2dotn = 0.
                   temp_w1dotn = 0.
                   temp_pdotn  = 0.
                   temp_u2dotn = (block(g)%x1(i+1) - &
                        block(g)%xcent(block(g)%nelu2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(index_ip1))) +  &
                     (block(g)%yp(j)-block(g)%ycent(block(g)%nelu2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(index_ip1)))  +  &
                     (block(g)%zp(k)-block(g)%zcent(block(g)%nelu2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(index_ip1)))
                   temp_v2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(index_ip1))) +  &
                 (block(g)%y1(j+1) - block(g)%ycent(block(g)%nelv2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(index_ip1)))  +  &
                    (block(g)%zp(k)- block(g)%zcent(block(g)%nelv2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(index_ip1)))
                   temp_v1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv1(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(index_ip1))) +  &
                     (block(g)%y1(j)-block(g)%ycent(block(g)%nelv1(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(index_ip1)))  +  &
                     (block(g)%zp(k)-block(g)%zcent(block(g)%nelv1(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(index_ip1)))
                   temp_w2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(index_ip1))) +  &
                     (block(g)%yp(j)-block(g)%ycent(block(g)%nelw2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(index_ip1)))  +  &
                   (block(g)%z1(k+1)-block(g)%zcent(block(g)%nelw2(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(index_ip1)))
                   temp_w1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw1(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(index_ip1))) +  &
                  (block(g)%yp(j)-block(g)%ycent(block(g)%nelw1(block(g)%index_ts(index_ip1))))*&
                        block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(index_ip1)))  +  &
                    (block(g)%z1(k)-block(g)%zcent(block(g)%nelw1(block(g)%index_ts(index_ip1))))*&
                       block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(index_ip1)))
                   temp_pdotn  = (block(g)%xp(i)- &
                        block(g)%xcent(block(g)%nelp(block(g)%index_ts(index_ip1))))&
                        *block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(index_ip1)))  +  &
                   (block(g)%yp(j)   - block(g)%ycent(block(g)%nelp(block(g)%index_ts(index_ip1))))&
                        *block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(index_ip1)))   +  &
                   (block(g)%zp(k)- block(g)%zcent(block(g)%nelp(block(g)%index_ts(index_ip1)))) &
                        *block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(index_ip1)))
                   IF (temp_u2dotn<0) THEN
                      block(g)%u(i+1,j,k) = block(g)%u2_ghost(index_ip1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i+1,j,k) = block(g)%v2_ghost(index_ip1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i+1,j-1,k) = block(g)%v1_ghost(index_ip1)
                   ENDIF
                   IF (temp_w2dotn<0) THEN
                      block(g)%w(i+1,j,k) = block(g)%w2_ghost(index_ip1)
                   ENDIF
                   IF (temp_w1dotn<0) THEN
                      block(g)%w(i+1,j,k-1) = block(g)%w1_ghost(index_ip1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i+1,j,k) = block(g)%p_ghost(index_ip1)
                   ENDIF
                ENDIF

                IF (block(g)%cell2(i-1, j, k)==2) THEN
                   index_im1 = nn
                   temp_u2dotn = 0.
                   temp_u1dotn = 0.
                   temp_v2dotn = 0.
                   temp_v1dotn = 0.
                   temp_w2dotn = 0.
                   temp_w1dotn = 0.
                   temp_pdotn  = 0.
                   temp_u1dotn = (block(g)%x1(i)   - &
                        block(g)%xcent(block(g)%nelu1(block(g)%index_ts(index_im1))))*&
                        block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(index_im1))) +  &
                     (block(g)%yp(j)-block(g)%ycent(block(g)%nelu1(block(g)%index_ts(index_im1))))&
                        *block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(index_im1)))  +  &
                    (block(g)%zp(k)-block(g)%zcent(block(g)%nelu1(block(g)%index_ts(index_im1))))*&
                        block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(index_im1)))
                   temp_v2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv2(block(g)%index_ts(index_im1))))*&
                        block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(index_im1))) +  &
                   (block(g)%y1(j+1)-block(g)%ycent(block(g)%nelv2(block(g)%index_ts(index_im1))))*&
                        block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(index_im1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelv2(block(g)%index_ts(index_im1))))*&
                        block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(index_im1)))
                   temp_v1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv1(block(g)%index_ts(index_im1))))*&
                        block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(index_im1))) +  &
                 (block(g)%y1(j)   - block(g)%ycent(block(g)%nelv1(block(g)%index_ts(index_im1))))*&
                   block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(index_im1)))  +  &
                   (block(g)%zp(k)   - &
                   block(g)%zcent(block(g)%nelv1(block(g)%index_ts(index_im1))))*&
                   block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(index_im1)))
                   temp_w2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw2(block(g)%index_ts(index_im1))))*&
                        block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(index_im1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw2(block(g)%index_ts(index_im1))))*&
                        block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(index_im1)))  +  &
                 (block(g)%z1(k+1) - block(g)%zcent(block(g)%nelw2(block(g)%index_ts(index_im1))))*&
                        block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(index_im1)))
                   temp_w1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw1(block(g)%index_ts(index_im1))))*&
                        block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(index_im1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw1(block(g)%index_ts(index_im1))))*&
                        block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(index_im1)))  +  &
                 (block(g)%z1(k)   - block(g)%zcent(block(g)%nelw1(block(g)%index_ts(index_im1))))*&
                        block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(index_im1)))
                   temp_pdotn  = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelp(block(g)%index_ts(index_im1))))*&
                        block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(index_im1)))  +  &
                  (block(g)%yp(j)   - block(g)%ycent(block(g)%nelp(block(g)%index_ts(index_im1))))*&
                       block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(index_im1)))   +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelp(block(g)%index_ts(index_im1))))*&
                       block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(index_im1)))
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-2,j,k) = block(g)%u1_ghost(index_im1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i-1,j,k) = block(g)%v2_ghost(index_im1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i-1,j-1,k) = block(g)%v1_ghost(index_im1)
                   ENDIF
                           IF (temp_w2dotn<0) THEN
                      block(g)%w(i-1,j,k) = block(g)%w2_ghost(index_im1)
                   ENDIF
                           IF (temp_w1dotn<0) THEN
                      block(g)%w(i-1,j,k-1) = block(g)%w1_ghost(index_im1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i-1,j,k) = block(g)%p_ghost(index_im1)
                   ENDIF
                ENDIF

                IF (block(g)%cell2(i, j+1, k)==2) THEN
                   index_jp1 = nn
                   temp_u2dotn = 0.
                   temp_u1dotn = 0.
                   temp_v2dotn = 0.
                   temp_v1dotn = 0.
                           temp_w2dotn = 0.
                   temp_w1dotn = 0.
                   temp_pdotn  = 0.
                   temp_u2dotn = (block(g)%x1(i+1) - &
                        block(g)%xcent(block(g)%nelu2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(index_jp1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelu2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(index_jp1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(index_jp1)))
                   temp_u1dotn = (block(g)%x1(i)   - &
                        block(g)%xcent(block(g)%nelu1(block(g)%index_ts(index_jp1))))*&
                   block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(index_jp1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelu1(block(g)%index_ts(index_jp1))))*&
                   block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(index_jp1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu1(block(g)%index_ts(index_jp1))))*&
                   block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(index_jp1)))
                   temp_v2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(index_jp1))) +  &
                 (block(g)%y1(j+1) - block(g)%ycent(block(g)%nelv2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(index_jp1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelv2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(index_jp1)))
                   temp_w2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(index_jp1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(index_jp1)))  +  &
                 (block(g)%z1(k+1) - block(g)%zcent(block(g)%nelw2(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(index_jp1)))
                   temp_w1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw1(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(index_jp1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw1(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(index_jp1)))  +  &
                 (block(g)%z1(k)   - block(g)%zcent(block(g)%nelw1(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(index_jp1)))
                   temp_pdotn  = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelp(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(index_jp1)))  +  &
                  (block(g)%yp(j)   - block(g)%ycent(block(g)%nelp(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(index_jp1)))   +  &
                  (block(g)%zp(k)   - block(g)%zcent(block(g)%nelp(block(g)%index_ts(index_jp1))))*&
                        block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(index_jp1)))
                   IF (temp_u2dotn<0) THEN
                      block(g)%u(i,j+1,k) = block(g)%u2_ghost(index_jp1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j+1,k) = block(g)%u1_ghost(index_jp1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i,j+1,k) = block(g)%v2_ghost(index_jp1)
                   ENDIF
                           IF (temp_w2dotn<0) THEN
                      block(g)%w(i,j+1,k) = block(g)%w2_ghost(index_jp1)
                   ENDIF
                           IF (temp_w1dotn<0) THEN
                      block(g)%w(i,j+1,k-1) = block(g)%w1_ghost(index_jp1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j+1,k) = block(g)%p_ghost(index_jp1)
                   ENDIF
                ENDIF

                        IF (block(g)%cell2(i, j-1, k)==2) THEN
                   index_jm1 = nn
                   temp_u2dotn = 0.
                   temp_u1dotn = 0.
                   temp_v2dotn = 0.
                   temp_v1dotn = 0.
                           temp_w2dotn = 0.
                   temp_w1dotn = 0.
                   temp_pdotn  = 0.
                   temp_u2dotn = (block(g)%x1(i+1) - &
                        block(g)%xcent(block(g)%nelu2(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(index_jm1))) +  &
                 (block(g)%yp(j) - block(g)%ycent(block(g)%nelu2(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(index_jm1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu2(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(index_jm1)))
                   temp_u1dotn = (block(g)%x1(i)   - &
                        block(g)%xcent(block(g)%nelu1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(index_jm1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelu1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(index_jm1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(index_jm1)))
                   temp_v1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(index_jm1))) +  &
                 (block(g)%y1(j)   - block(g)%ycent(block(g)%nelv1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(index_jm1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelv1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(index_jm1)))
                   temp_w2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw2(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(index_jm1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw2(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(index_jm1)))  +  &
                 (block(g)%z1(k+1) - block(g)%zcent(block(g)%nelw2(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(index_jm1)))
                   temp_w1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(index_jm1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(index_jm1)))  +  &
                 (block(g)%z1(k)   - block(g)%zcent(block(g)%nelw1(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(index_jm1)))
                   temp_pdotn  = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelp(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(index_jm1)))  +  &
                  (block(g)%yp(j)   - block(g)%ycent(block(g)%nelp(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(index_jm1)))   +  &
                  (block(g)%zp(k)   - block(g)%zcent(block(g)%nelp(block(g)%index_ts(index_jm1))))*&
                        block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(index_jm1)))
                   IF (temp_u2dotn<0) THEN
                      block(g)%u(i,j-1,k) = block(g)%u2_ghost(index_jm1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j-1,k) = block(g)%u1_ghost(index_jm1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i,j-2,k) = block(g)%v1_ghost(index_jm1)
                   ENDIF
                   IF (temp_w2dotn<0) THEN
                      block(g)%w(i,j-1,k) = block(g)%w2_ghost(index_jm1)
                   ENDIF
                   IF (temp_w1dotn<0) THEN
                      block(g)%w(i,j-1,k-1) = block(g)%w1_ghost(index_jm1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j-1,k) = block(g)%p_ghost(index_jm1)
                   ENDIF
                ENDIF

                IF (block(g)%cell2(i, j, k+1)==2) THEN
                   index_kp1 = nn
                   temp_u2dotn = 0.
                   temp_u1dotn = 0.
                   temp_v2dotn = 0.
                   temp_v1dotn = 0.
                           temp_w2dotn = 0.
                   temp_w1dotn = 0.
                   temp_pdotn  = 0.
                   temp_u2dotn = (block(g)%x1(i+1) - &
                        block(g)%xcent(block(g)%nelu2(block(g)%index_ts(index_kp1))))*&
                        block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(index_kp1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelu2(block(g)%index_ts(index_kp1))))*&
                        block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(index_kp1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu2(block(g)%index_ts(index_kp1))))*&
                        block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(index_kp1)))
                   temp_u1dotn = (block(g)%x1(i)   - &
                        block(g)%xcent(block(g)%nelu1(block(g)%index_ts(index_kp1))))*&
                   block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(index_kp1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelu1(block(g)%index_ts(index_kp1))))*&
                        block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(index_kp1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu1(block(g)%index_ts(index_kp1))))*&
                        block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(index_kp1)))
                     temp_v2dotn = (block(g)%xp(i)   - &
                          block(g)%xcent(block(g)%nelv2(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(index_kp1))) +  &
                 (block(g)%y1(j+1) - block(g)%ycent(block(g)%nelv2(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(index_kp1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelv2(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(index_kp1)))
                     temp_v1dotn = (block(g)%xp(i)   - &
                          block(g)%xcent(block(g)%nelv1(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(index_kp1))) +  &
                 (block(g)%y1(j)   - block(g)%ycent(block(g)%nelv1(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(index_kp1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelv1(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(index_kp1)))
                     temp_w2dotn = (block(g)%xp(i)   - &
                          block(g)%xcent(block(g)%nelw2(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(index_kp1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw2(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(index_kp1)))  +  &
                 (block(g)%z1(k+1) - block(g)%zcent(block(g)%nelw2(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(index_kp1)))
                     temp_pdotn  = (block(g)%xp(i)   - &
                          block(g)%xcent(block(g)%nelp(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(index_kp1)))  +  &
                  (block(g)%yp(j)   - block(g)%ycent(block(g)%nelp(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(index_kp1)))   +  &
                  (block(g)%zp(k)   - block(g)%zcent(block(g)%nelp(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(index_kp1)))
                   IF (temp_u2dotn<0) THEN
                      block(g)%u(i,j,k+1) = block(g)%u2_ghost(index_kp1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j,k+1) = block(g)%u1_ghost(index_kp1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i,j,k+1) = block(g)%v2_ghost(index_kp1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i,j-1,k+1) = block(g)%v1_ghost(index_kp1)
                   ENDIF
                   IF (temp_w2dotn<0) THEN
                      block(g)%w(i,j,k+1) = block(g)%w2_ghost(index_kp1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j,k+1) = block(g)%p_ghost(index_kp1)
                   ENDIF
                ENDIF

                IF (block(g)%cell2(i, j, k-1)==2) THEN
                   index_km1 = nn
                   temp_u2dotn = 0.
                   temp_u1dotn = 0.
                   temp_v2dotn = 0.
                   temp_v1dotn = 0.
                   temp_w2dotn = 0.
                   temp_w1dotn = 0.
                   temp_pdotn  = 0.
                   temp_u2dotn = (block(g)%x1(i+1) - &
                        block(g)%xcent(block(g)%nelu2(block(g)%index_ts(index_km1))))*&
                        block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(index_km1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelu2(block(g)%index_ts(index_km1))))*&
                        block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(index_km1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu2(block(g)%index_ts(index_km1))))*&
                        block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(index_km1)))
                   temp_u1dotn = (block(g)%x1(i)   - &
                        block(g)%xcent(block(g)%nelu1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(index_km1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelu1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(index_km1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelu1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(index_km1)))
                   temp_v2dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv2(block(g)%index_ts(index_km1))))*&
                        block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(index_km1))) +  &
                 (block(g)%y1(j+1) - block(g)%ycent(block(g)%nelv2(block(g)%index_ts(index_km1))))*&
                        block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(index_km1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelv2(block(g)%index_ts(index_km1))))*&
                        block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(index_km1)))
                   temp_v1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelv1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(index_km1))) +  &
                 (block(g)%y1(j)   - block(g)%ycent(block(g)%nelv1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(index_km1)))  +  &
                 (block(g)%zp(k)   - block(g)%zcent(block(g)%nelv1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(index_km1)))
                   temp_w1dotn = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelw1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(index_km1))) +  &
                 (block(g)%yp(j)   - block(g)%ycent(block(g)%nelw1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(index_km1)))  +  &
                 (block(g)%z1(k)   - block(g)%zcent(block(g)%nelw1(block(g)%index_ts(index_km1))))*&
                        block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(index_km1)))
                   temp_pdotn  = (block(g)%xp(i)   - &
                        block(g)%xcent(block(g)%nelp(block(g)%index_ts(index_km1))))*&
                        block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(index_km1)))  +  &
                  (block(g)%yp(j)   - block(g)%ycent(block(g)%nelp(block(g)%index_ts(index_km1))))*&
                        block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(index_km1)))   +  &
                  (block(g)%zp(k)   - block(g)%zcent(block(g)%nelp(block(g)%index_ts(index_km1))))*&
                        block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(index_km1)))
                   IF (temp_u2dotn<0) THEN
                      block(g)%u(i,j,k-1) = block(g)%u2_ghost(index_km1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j,k-1) = block(g)%u1_ghost(index_km1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i,j,k-1) = block(g)%v2_ghost(index_km1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i,j-1,k-1) = block(g)%v1_ghost(index_km1)
                   ENDIF
                   IF (temp_w1dotn<0) THEN
                      block(g)%w(i,j,k-2) = block(g)%w1_ghost(index_km1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j,k-1) = block(g)%p_ghost(index_km1)
                   ENDIF
                ENDIF

             ENDDO
      ENDIF
          dxr=block(g)%deltax(i+1)
          dx=block(g)%deltax(i)
          dxl=block(g)%deltax(i-1)
          dye=block(g)%deltay(j+1)
          dy=block(g)%deltay(j)
          dyw=block(g)%deltay(j-1)
          dzt=block(g)%deltaz(k+1)
          dz=block(g)%deltaz(k)
          dzb=block(g)%deltaz(k-1)
!cccccccccccccccccccccccccc  diff-u     ccccccccccccccccccccccccccccccccc
!c     duu / dx
       u1a = block(g)%u(i-1,j,k) + block(g)%u(i,j,k)
       u22 = block(g)%u(i-1,j,k) - block(g)%u(i,j,k)
       u3  = block(g)%u(i,j,k)   + block(g)%u(i+1,j,k)
       u4  = block(g)%u(i,j,k)   - block(g)%u(i+1,j,k)

!c     duv  / dy
       u5 = block(g)%u(i,j-1,k)  + block(g)%u(i,j,k)
       u6 = block(g)%u(i,j-1,k)  - block(g)%u(i,j,k)
       u7 = block(g)%u(i,j,k)    + block(g)%u(i,j+1,k)
       u8 = block(g)%u(i,j,k)    - block(g)%u(i,j+1,k)

!c     dwu  /  dz
       u9 = block(g)%u(i,j,k-1)  + block(g)%u(i,j,k)
       u10= block(g)%u(i,j,k-1)  - block(g)%u(i,j,k)
       u11= block(g)%u(i,j,k)    + block(g)%u(i,j,k+1)
       u12= block(g)%u(i,j,k)    - block(g)%u(i,j,k+1)

!c     duv/dx
       u13 = block(g)%u(i-1,j,k) + block(g)%u(i-1,j+1,k)
       u14 = u7

!c     duw / dx
       u15 = block(g)%u(i-1,j,k) + block(g)%u(i-1,j,k+1)
       u16 = u11
!cccccccccccccccccccccccccc    diff -v    ccccccccccccccccccccccccccccccc
!c     dvu / dx
       v1a = block(g)%v(i,j-1,k)  + block(g)%v(i+1,j-1,k)
       v22 = block(g)%v(i,j,k)    + block(g)%v(i+1,j,k)

!c     duv / dx
       v3 = block(g)%v(i-1,j,k)   + block(g)%v(i,j,k)
       v4 = block(g)%v(i-1,j,k)   - block(g)%v(i,j,k)
       v5 = block(g)%v(i,j,k)     + block(g)%v(i+1,j,k)
       v6 = block(g)%v(i,j,k)     - block(g)%v(i+1,j,k)

!c     dvv / dy
       v7 = block(g)%v(i,j-1,k)   + block(g)%v(i,j,k)
       v8 = block(g)%v(i,j-1,k)   - block(g)%v(i,j,k)
       v9 = block(g)%v(i,j,k)     + block(g)%v(i,j+1,k)
       v10= block(g)%v(i,j,k)     - block(g)%v(i,j+1,k)

!c     dwu / dz
       v11 = block(g)%v(i,j,k-1)  + block(g)%v(i,j,k)
       v12 = block(g)%v(i,j,k-1)  - block(g)%v(i,j,k)
       v13 = block(g)%v(i,j,k)    + block(g)%v(i,j,k+1)
       v14 = block(g)%v(i,j,k)    - block(g)%v(i,j,k+1)

!c     dvw / dy
       v15 = block(g)%v(i,j-1,k)  + block(g)%v(i,j-1,k+1)
       v16 = v13
!cccccccccccccccccccccccccccc  diff - w  cccccccccccccccccccccccccccccccc
!c     dwu / dz
       w1a = block(g)%w(i,j,k-1)  + block(g)%w(i+1,j,k-1)
       w22 = block(g)%w(i,j,k)    + block(g)%w(i+1,j,k)
!c
!c     dwv / dz
       w3 = block(g)%w(i,j,k-1)   + block(g)%w(i,j+1,k-1)
       w4 = block(g)%w(i,j,k)     + block(g)%w(i,j+1,k)

!c     duw / dx
       w5 = block(g)%w(i-1,j,k)   + block(g)%w(i,j,k)
       w6 = block(g)%w(i-1,j,k)   - block(g)%w(i,j,k)
       w7 = block(g)%w(i,j,k)     + block(g)%w(i+1,j,k)
       w8 = block(g)%w(i,j,k)     - block(g)%w(i+1,j,k)

!c     dvw / dy
       w9 = block(g)%w(i,j-1,k)   + block(g)%w(i,j,k)
       w10 = block(g)%w(i,j-1,k)  - block(g)%w(i,j,k)
       w11 = block(g)%w(i,j,k)    + block(g)%w(i,j+1,k)
       w12 = block(g)%w(i,j,k)    - block(g)%w(i,j+1,k)

!c     dww / dz
       w13 = block(g)%w(i,j,k-1)  + block(g)%w(i,j,k)
       w14 = block(g)%w(i,j,k-1)  - block(g)%w(i,j,k)
       w15 = block(g)%w(i,j,k)    + block(g)%w(i,j,k+1)
       w16 = block(g)%w(i,j,k)    - block(g)%w(i,j,k+1)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       dpdx = (block(g)%p(i,j,k) - block(g)%p(i+1,j,k))/(0.5_dp*(dxr+dx))
       dpdy = (block(g)%p(i,j,k) - block(g)%p(i,j+1,k))/(0.5_dp*(dye+dy))
       dpdz = (block(g)%p(i,j,k) - block(g)%p(i,j,k+1))/(0.5_dp*(dzt+dz))

            dx2xr = dx + dxr
       dx2xl = dx + dxl
           dy2ye = dy + dye
           dy2yw = dy + dyw
           dz2zt = dz + dzt
           dz2zb = dz + dzb
!c*********************** U - Momentum **********************************
           vu_e=block(g)%v(i+1,j,k)+(dxr/(dx2xr))*(block(g)%v(i,j,k)-block(g)%v(i+1,j,k))
           vu_w=block(g)%v(i+1,j-1,k)+(dxr/(dx2xr))*(block(g)%v(i,j-1,k)-block(g)%v(i+1,j-1,k))
           v_in_um=0.5_dp*(vu_e+vu_w)

           wu_n=block(g)%w(i+1,j,k)+(dxr/(dx2xr))*(block(g)%w(i,j,k)-block(g)%w(i+1,j,k))
           wu_s=block(g)%w(i+1,j,k-1)+(dxr/(dx2xr))*(block(g)%w(i,j,k-1)-block(g)%w(i+1,j,k-1))
           w_in_um=0.5_dp*(wu_n+wu_s)

!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
       if(i>2.and.i<nx_var.and.j>2.and.j<ny_var+1.and. &
      k>2.and.k<nz_var+1.and.block(g)%cell(i+1,j,k)/=2.and.     &
      block(g)%cell(i-1,j,k)/=2.and.block(g)%cell(i,j+1,k)/=2.and.       &
      block(g)%cell(i,j-1,k)/=2.and.block(g)%cell(i,j,k+1)/=2.and.       &
      block(g)%cell(i,j,k-1)/=2.and.block(g)%cell2(i+2,j,k)/=2.and. &
      block(g)%cell2(i,j+2,k)/=2.and.block(g)%cell2(i,j,k+2)/=2.and. &
      block(g)%cell2(i-2,j,k)/=2.and.block(g)%cell2(i,j-2,k)/=2.and. &
      block(g)%cell2(i,j,k-2)/=2) THEN


           ddy=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j-1))
       ddye=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j+1))
       ddz=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
       ddzr=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k+1))

       duutdx=block(g)%u(i,j,k)*(block(g)%ca_uu(1, i)*block(g)%u(i+2,j,k)+&
              block(g)%ca_uu(2, i)*block(g)%u(i+1,j,k) &
             +block(g)%ca_uu(3, i)*block(g)%u(i,j,k)+block(g)%ca_uu(4, i)&
             *block(g)%u(i-1,j,k)+block(g)%ca_uu(5, i)*&
      block(g)%u(i-2,j,k))/(block(g)%ca_uu(6, i)*block(g)%deltax(i+1))+dabs(block(g)%u(i,j,k))* &
      (block(g)%ck_uu(1, i)*block(g)%u(i+2,j,k)+block(g)%ck_uu(2, i)*block(g)%u(i+1,j,k) &
             +block(g)%ck_uu(3, i)*block(g)%u(i,j,k)+block(g)%ck_uu(4, i)&
             *block(g)%u(i-1,j,k)+block(g)%ck_uu(5, i)*&
      block(g)%u(i-2,j,k))/(2.0_dp*block(g)%ck_uu(6, i)*block(g)%deltax(i))

       dvutdy=v_in_um*(block(g)%ca_vu(1, j)*block(g)%u(i,j+2,k)+&
            block(g)%ca_vu(2, j)*block(g)%u(i,j+1,k) &
            +block(g)%ca_vu(3, j)*block(g)%u(i,j,k)+block(g)%ca_vu(4, j)*&
            block(g)%u(i,j-1,k)+block(g)%ca_vu(5, j)* &
      block(g)%u(i,j-2,k))/(block(g)%ca_vu(6, j)*ddye)+dabs(v_in_um)* &
      (block(g)%ck_vu(1, j)*block(g)%u(i,j+2,k)+block(g)%ck_vu(2, j)*block(g)%u(i,j+1,k) &
      +block(g)%ck_vu(3, j)*block(g)%u(i,j,k)+block(g)%ck_vu(4, j)*&
      block(g)%u(i,j-1,k)+block(g)%ck_vu(5, j)* &
      block(g)%u(i,j-2,k))/(2.0_dp*block(g)%ck_vu(6, j)*ddy)

       dwutdz=w_in_um*(block(g)%ca_wu(1, k)*block(g)%u(i,j,k+2)+&
            block(g)%ca_wu(2, k)*block(g)%u(i,j,k+1) &
            +block(g)%ca_wu(3, k)*block(g)%u(i,j,k)+block(g)%ca_wu(4, k)*&
            block(g)%u(i,j,k-1)+block(g)%ca_wu(5, k)* &
      block(g)%u(i,j,k-2))/(block(g)%ca_wu(6, k)*ddzr)+dabs(w_in_um)* &
      (block(g)%ck_wu(1, k)*block(g)%u(i,j,k+2)+block(g)%ck_wu(2, k)*block(g)%u(i,j,k+1) &
      +block(g)%ck_wu(3, k)*block(g)%u(i,j,k)+block(g)%ck_wu(4, k)*block(g)%u(i,j,k-1)+&
      block(g)%ck_wu(5, k)* &
      block(g)%u(i,j,k-2))/(2.0_dp*block(g)%ck_wu(6, k)*ddz)

            duudx=duutdx
       dvudy=dvutdy
       dwudz=dwutdz

       else
 !cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
   duudx = -0.25_dp*( u1a*u1a + alpha*dabs(u1a)*u22-u3*u3-alpha*dabs(u3)*u4)/block(g)%deltax(i)

   dvudy = -0.25_dp*( v1a*u5 + alpha*dabs(v1a)*u6 - v22*u7 - alpha*dabs(v22)*u8)/block(g)%deltay(j)

   dwudz = -0.25_dp*(w1a*u9 + alpha*dabs(w1a)*u10-w22*u11 - alpha*dabs(w22)*u12) /block(g)%deltaz(k)

       endif

!ccccccccccccccccccccccccc  grad of u part cccccccccccccccccccccccccccccc
       d2udx2=(2.0_dp/dx2xr)*((-u4/dxr)+(u22/dx))

       d2udy2=(2.0_dp/(0.5_dp*(dy2ye+dy2yw)))*((-u8/(0.5_dp*dy2ye))+ &
      (u6/(0.5_dp*dy2yw)))

       d2udz2=(2.0_dp/(0.5_dp*(dz2zt+dz2zb)))*((-u12/(0.5_dp*dz2zt))+ &
      (u10/(0.5_dp*dz2zb)))

       xtt2=(d2udx2+d2udy2+d2udz2)/re

       residu=(-duudx-dvudy-dwudz+xtt2)
        block(g)%ut(i,j,k)=block(g)%u(i,j,k)+deltat*(residu+dpdx)

!c*********************** V - Momentum *********************************
       uv_e=block(g)%u(i,j+1,k)+(dye/(dy2ye))*(block(g)%u(i,j,k)-block(g)%u(i,j+1,k))
       uv_w=block(g)%u(i-1,j+1,k)+(dye/(dy2ye))*(block(g)%u(i-1,j,k)-block(g)%u(i-1,j+1,k))
       u_in_vm=0.5_dp*(uv_e+uv_w)

       wv_n=block(g)%w(i,j+1,k)+(dye/(dy2ye))*(block(g)%w(i,j,k)-block(g)%w(i,j+1,k))
       wv_s=block(g)%w(i,j+1,k-1)+(dye/(dy2ye))*(block(g)%w(i,j,k-1)-block(g)%w(i,j+1,k-1))
       w_in_vm=0.5_dp*(wv_n+wv_s)

!cccccccccccccc---Third Order Upwinding ----cccccccccccccccccccccccccccc
     if(i>2.and.i<nx_var+1.and.j>2.and.j<ny_var.and. &
      k>2.and.k<nz_var+1.and.block(g)%cell(i+1,j,k)/=2.and.  &
      block(g)%cell(i-1,j,k)/=2.and.block(g)%cell(i,j+1,k)/=2.and. &
      block(g)%cell(i,j-1,k)/=2.and.block(g)%cell(i,j,k+1)/=2.and. &
      block(g)%cell(i,j,k-1)/=2.and.block(g)%cell2(i+2,j,k)/=2.and. &
      block(g)%cell2(i,j+2,k)/=2.and.block(g)%cell2(i,j,k+2)/=2.and. &
      block(g)%cell2(i-2,j,k)/=2.and.block(g)%cell2(i,j-2,k)/=2.and. &
      block(g)%cell2(i,j,k-2)/=2) THEN
       ddx=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i-1))
       ddxr=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i+1))
       ddz=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
       ddzr=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k+1))

       ! TODO: Note that this calculation, which is duvt/dx, is using
       ! ca_uw not ca_uv... ca_uw == ca_uv, however, I think for
       ! completeness it would make more sense to use e.g. ca_uv
       ! (which isn't used anywhere except assignment)
     duvtdx=u_in_vm*(block(g)%ca_uw(1, i)*block(g)%v(i+2,j,k)&
           +block(g)%ca_uw(2, i)*block(g)%v(i+1,j,k) &
           +block(g)%ca_uw(3, i)*block(g)%v(i,j,k)+block(g)%ca_uw(4, i)&
           *block(g)%v(i-1,j,k)+block(g)%ca_uw(5, i)* &
      block(g)%v(i-2,j,k))/(block(g)%ca_uw(6, i)*ddxr)+dabs(u_in_vm)* &
      (block(g)%ck_uw(1, i)*block(g)%v(i+2,j,k)+block(g)%ck_uw(2, i)*block(g)%v(i+1,j,k) &
      +block(g)%ck_uw(3, i)*block(g)%v(i,j,k)+block(g)%ck_uw(4, i)*block(g)%v(i-1,j,k)&
      +block(g)%ck_uw(5, i)* &
      block(g)%v(i-2,j,k))/(2.0_dp*block(g)%ck_uw(6, i)*ddx)

     dvvtdy=block(g)%v(i,j,k)*(block(g)%ca_vv(1, j)*block(g)%v(i,j+2,k)+&
          block(g)%ca_vv(2, j)*block(g)%v(i,j+1,k) &
          +block(g)%ca_vv(3, j)*block(g)%v(i,j,k)+block(g)%ca_vv(4, j)*block(g)%v(i,j-1,k)+&
          block(g)%ca_vv(5, j)* &
      block(g)%v(i,j-2,k))/(block(g)%ca_vv(6, j)*block(g)%deltay(j+1))+dabs(block(g)%v(i,j,k))* &
      (block(g)%ck_vv(1, j)*block(g)%v(i,j+2,k)+block(g)%ck_vv(2, j)*block(g)%v(i,j+1,k) &
      +block(g)%ck_vv(3, j)*block(g)%v(i,j,k)+block(g)%ck_vv(4, j)*block(g)%v(i,j-1,k)+&
      block(g)%ck_vv(5, j)* &
      block(g)%v(i,j-2,k))/(2.0_dp*block(g)%ck_vv(6, j)*block(g)%deltay(j))

       ! See TODO about mismatch of variables... this is dwv but is
       ! using e.g. ca_wu rather than wv
     dwvtdz=w_in_vm*(block(g)%ca_wu(1, k)*block(g)%v(i,j,k+2)+block(g)%ca_wu(2, k)&
           *block(g)%v(i,j,k+1) &
           +block(g)%ca_wu(3, k)*block(g)%v(i,j,k)+block(g)%ca_wu(4, k)&
           *block(g)%v(i,j,k-1)+block(g)%ca_wu(5, k)* &
      block(g)%v(i,j,k-2))/(block(g)%ca_wu(6, k)*ddzr)+dabs(w_in_vm)* &
      (block(g)%ck_wu(1, k)*block(g)%v(i,j,k+2)+block(g)%ck_wu(2, k)*block(g)%v(i,j,k+1) &
      +block(g)%ck_wu(3, k)*block(g)%v(i,j,k)+block(g)%ck_wu(4, k)*block(g)%v(i,j,k-1)+&
      block(g)%ck_wu(5, k)* &
      block(g)%v(i,j,k-2))/(2.0_dp*block(g)%ck_wu(6, k)*ddz)

       duvdx=duvtdx
       dvvdy=dvvtdy
       dwvdz=dwvtdz

       else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
      duvdx = -0.25_dp*(u13*v3 + alpha*dabs(u13)*v4- u14*v5 - alpha*dabs(u14)*v6)/block(g)%deltax(i)

      dvvdy = -0.25_dp*(v7*v7 + alpha*dabs(v7)*v8 -v9*v9 - alpha*dabs(v9)*v10)/block(g)%deltay(j)

      dwvdz = -0.25_dp*(w3*v11 + alpha*dabs(w3)*v12 -w4*v13 - alpha*dabs(w4)*v14)/block(g)%deltaz(k)

       endif

!ccccccccccccccccccccccccc  grad of v part cccccccccccccccccccccccccccccc
       d2vdx2=(2.0_dp/(0.5_dp*(dx2xr+dx2xl)))*((-v6/(0.5_dp*dx2xr))+ &
       (v4/(0.5_dp*dx2xl)))

       d2vdy2=(2.0_dp/dy2ye)*((-v10/dye)+(v8/dy))

       d2vdz2=(2.0_dp/(0.5_dp*(dz2zt+dz2zb)))*((-v14/(0.5_dp*dz2zt))+ &
       (v12/(0.5_dp*dz2zb)))

       ytt2=(d2vdx2+d2vdy2+d2vdz2)/re

       residv=(-duvdx-dvvdy-dwvdz+ytt2)

       block(g)%vt(i,j,k)=block(g)%v(i,j,k)+deltat*(residv+dpdy)

!c*********************** W - Momentum **********************************
       uw_e=block(g)%u(i,j,k+1)+(dzt/(dz2zt))*(block(g)%u(i,j,k)-block(g)%u(i,j,k+1))
       uw_w=block(g)%u(i-1,j,k+1)+(dzt/(dz2zt))*(block(g)%u(i-1,j,k)-block(g)%u(i-1,j,k+1))
       u_in_wm=0.5_dp*(uw_e+uw_w)

       vw_n=block(g)%v(i,j,k+1)+(dzt/(dz2zt))*(block(g)%v(i,j,k)-block(g)%v(i,j,k+1))
       vw_s=block(g)%v(i,j-1,k+1)+(dzt/(dz2zt))*(block(g)%v(i,j-1,k)-block(g)%v(i,j-1,k+1))
       v_in_wm=0.5_dp*(vw_n+vw_s)

!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
     if(i>2.and.i<nx_var+1.and.j>2.and.j<ny_var+1.and. &
      k>2.and.k<nz_var.and.block(g)%cell(i+1,j,k)/=2.and. &
      block(g)%cell(i-1,j,k)/=2.and.block(g)%cell(i,j+1,k)/=2.and. &
      block(g)%cell(i,j-1,k)/=2.and.block(g)%cell(i,j,k+1)/=2.and. &
      block(g)%cell(i,j,k-1)/=2.and.block(g)%cell2(i+2,j,k)/=2.and. &
      block(g)%cell2(i,j+2,k)/=2.and.block(g)%cell2(i,j,k+2)/=2.and. &
      block(g)%cell2(i-2,j,k)/=2.and.block(g)%cell2(i,j-2,k)/=2.and. &
      block(g)%cell2(i,j,k-2)/=2) THEN

       ddx=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i-1))
       ddxr=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i+1))
       ddy=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j-1))
       ddye=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j+1))

       duwtdx=u_in_wm*(block(g)%ca_uw(1, i)*block(g)%w(i+2,j,k)+block(g)%ca_uw(2, i)*&
            block(g)%w(i+1,j,k) &
            +block(g)%ca_uw(3, i)*block(g)%w(i,j,k)+block(g)%ca_uw(4, i)*block(g)%w(i-1,j,k)+&
            block(g)%ca_uw(5, i)* &
      block(g)%w(i-2,j,k))/(block(g)%ca_uw(6, i)*ddxr)+dabs(u_in_wm)* &
      (block(g)%ck_uw(1, i)*block(g)%w(i+2,j,k)+block(g)%ck_uw(2, i)*block(g)%w(i+1,j,k) &
      +block(g)%ck_uw(3, i)*block(g)%w(i,j,k)+block(g)%ck_uw(4, i)*block(g)%w(i-1,j,k)+&
      block(g)%ck_uw(5, i)* &
      block(g)%w(i-2,j,k))/(2.0_dp*block(g)%ck_uw(6, i)*ddx)

       dvwtdy=v_in_wm*(block(g)%ca_vw(1, j)*block(g)%w(i,j+2,k)+block(g)%ca_vw(2, j)*&
            block(g)%w(i,j+1,k) &
            +block(g)%ca_vw(3, j)*block(g)%w(i,j,k)+block(g)%ca_vw(4, j)*block(g)%w(i,j-1,k)+&
            block(g)%ca_vw(5, j)* &
      block(g)%w(i,j-2,k))/(block(g)%ca_vw(6, j)*ddye)+dabs(v_in_wm)* &
      (block(g)%ck_vw(1, j)*block(g)%w(i,j+2,k)+block(g)%ck_vw(2, j)*block(g)%w(i,j+1,k) &
      +block(g)%ck_vw(3, j)*block(g)%w(i,j,k)+block(g)%ck_vw(4, j)*block(g)%w(i,j-1,k)+&
       block(g)%ck_vw(5, j)* &
      block(g)%w(i,j-2,k))/(2.0_dp*block(g)%ck_vw(6, j)*ddy)

       dwwtdz=block(g)%w(i,j,k)*(block(g)%ca_ww(1, k)*block(g)%w(i,j,k+2)+block(g)%ca_ww(2, k)*&
            block(g)%w(i,j,k+1) &
            +block(g)%ca_ww(3, k)*block(g)%w(i,j,k)+block(g)%ca_ww(4, k)*block(g)%w(i,j,k-1)+&
            block(g)%ca_ww(5, k)* &
      block(g)%w(i,j,k-2))/(block(g)%ca_ww(6, k)*block(g)%deltaz(k+1))+dabs(block(g)%w(i,j,k))* &
      (block(g)%ck_ww(1, k)*block(g)%w(i,j,k+2)+block(g)%ck_ww(2, k)*block(g)%w(i,j,k+1) &
      +block(g)%ck_ww(3, k)*block(g)%w(i,j,k)+block(g)%ck_ww(4, k)*block(g)%w(i,j,k-1)+&
      block(g)%ck_ww(5, k)* &
      block(g)%w(i,j,k-2))/(2.0_dp*block(g)%ck_ww(6, k)*block(g)%deltaz(k))

            duwdx=duwtdx
       dvwdy=dvwtdy
       dwwdz=dwwtdz

       else
!cccccccccccccc---First Order Upwinding ----cccccccccccccccccccccccccccc

   duwdx=-0.5_dp*(u15*w5 + alpha*dabs(u15)*w6 -  u16*w7 - alpha*dabs(u16)*w8)/block(g)%deltax(i)

   dvwdy=-0.5_dp*(v15*w9 + alpha*dabs(v15)*w10 - v16*w11 -alpha*dabs(v16)*w12)/block(g)%deltay(j)

   dwwdz=-0.5_dp*(w13*w13 + alpha*dabs(w13)*w14 - w15*w15 - alpha*dabs(w15)*w16)/block(g)%deltaz(k)

       endif

!ccccccccccccccccccccccccc  grad of w part cccccccccccccccccccccccccccccc

       d2wdx2=(2.0_dp/(0.5_dp*(dx2xr+dx2xl)))*((-W8/(0.5_dp*dx2xr))+ &
      (W6/(0.5_dp*dx2xl)))

       d2wdy2=(2.0_dp/(0.5_dp*(dy2ye+dy2yw)))*((-W12/(0.5_dp*dy2ye))+ &
      (W10/(0.5_dp*dy2yw)))

       d2wdz2=(2.0_dp/dz2zt)*((-w16/dzt)+(w14/dz))

       ztt2=(d2wdx2+d2wdy2+d2wdz2)/re

       residw =(-duwdx-dvwdy-dwwdz+ztt2)

       block(g)%wt(i,j,k)=block(g)%w(i,j,k)+deltat*(residw+dpdz)

!c***********************************************************************

IF (block(g)%cell2(i+1, j, k)==2) THEN
   block(g)%u(i+1,j,k)   = block(g)%u2t_ghost(index_ip1)
   block(g)%v(i+1,j,k)   = block(g)%v2t_ghost(index_ip1)
   block(g)%v(i+1,j-1,k) = block(g)%v1t_ghost(index_ip1)
   block(g)%w(i+1,j,k)   = block(g)%w2t_ghost(index_ip1)
   block(g)%w(i+1,j,k-1) = block(g)%w1t_ghost(index_ip1)
   block(g)%p(i+1,j,k)   = block(g)%pt_ghost(index_ip1)
ENDIF
IF (block(g)%cell2(i-1, j, k)==2) THEN
   block(g)%u(i-2,j,k)   = block(g)%u1t_ghost(index_im1)
   block(g)%v(i-1,j,k)   = block(g)%v2t_ghost(index_im1)
   block(g)%v(i-1,j-1,k) = block(g)%v1t_ghost(index_im1)
   block(g)%w(i-1,j,k)   = block(g)%w2t_ghost(index_im1)
   block(g)%w(i-1,j,k-1) = block(g)%w1t_ghost(index_im1)
   block(g)%p(i-1,j,k)   = block(g)%pt_ghost(index_im1)
ENDIF
IF (block(g)%cell2(i, j+1, k)==2) THEN
   block(g)%u(i,j+1,k)   = block(g)%u2t_ghost(index_jp1)
   block(g)%u(i-1,j+1,k) = block(g)%u1t_ghost(index_jp1)
   block(g)%v(i,j+1,k)   = block(g)%v2t_ghost(index_jp1)
   block(g)%w(i,j+1,k)   = block(g)%w2t_ghost(index_jp1)
   block(g)%w(i,j+1,k-1) = block(g)%w1t_ghost(index_jp1)
   block(g)%p(i,j+1,k)   = block(g)%pt_ghost(index_jp1)
ENDIF
IF (block(g)%cell2(i, j-1, k)==2) THEN
   block(g)%u(i,j-1,k)   = block(g)%u2t_ghost(index_jm1)
   block(g)%u(i-1,j-1,k) = block(g)%u1t_ghost(index_jm1)
   block(g)%v(i,j-2,k)   = block(g)%v1t_ghost(index_jm1)
   block(g)%w(i,j-1,k)   = block(g)%w2t_ghost(index_jm1)
   block(g)%w(i,j-1,k-1) = block(g)%w1t_ghost(index_jm1)
   block(g)%p(i,j-1,k)   = block(g)%pt_ghost(index_jm1)
ENDIF
IF (block(g)%cell2(i, j, k+1)==2) THEN
   block(g)%u(i,j,k+1)   = block(g)%u2t_ghost(index_kp1)
   block(g)%u(i-1,j,k+1) = block(g)%u1t_ghost(index_kp1)
   block(g)%v(i,j,k+1)   = block(g)%v2t_ghost(index_kp1)
   block(g)%v(i,j-1,k+1) = block(g)%v1t_ghost(index_kp1)
   block(g)%w(i,j,k+1)   = block(g)%w2t_ghost(index_kp1)
   block(g)%p(i,j,k+1)   = block(g)%pt_ghost(index_kp1)
ENDIF
IF (block(g)%cell2(i, j, k-1)==2) THEN
   block(g)%u(i,j,k-1)   = block(g)%u2t_ghost(index_km1)
   block(g)%u(i-1,j,k-1) = block(g)%u1t_ghost(index_km1)
   block(g)%v(i,j,k-1)   = block(g)%v2t_ghost(index_km1)
   block(g)%v(i,j-1,k-1) = block(g)%v1t_ghost(index_km1)
   block(g)%w(i,j,k-2)   = block(g)%w1t_ghost(index_km1)
   block(g)%p(i,j,k-1)   = block(g)%pt_ghost(index_km1)
ENDIF
       ENDDO
!c***********************************************************************
     !$acc end parallel loop
        !$acc wait

       ENDDO
      end subroutine nsMomentum2order


    pure function compute_f(theta) result(f)
        real(dp), intent(in) :: theta(3)
        real(dp) :: f(8)

        f(1) = theta(3)
        f(2) = 2.0_dp*theta(3)+theta(3)**2.0_dp
        f(3) = 3.0_dp*theta(3)+3.0_dp*theta(3)**2.0_dp+theta(3)**3.0_dp
        f(4) = 4.0_dp*theta(3)+6.0_dp*theta(3)**2.0_dp+4.0_dp*theta(3)**3.0_dp+theta(3)**4.0_dp
        f(5) = 1.0_dp/(theta(1)*theta(2))
        f(6) = (2.0_dp*theta(1)+1.0_dp)/((theta(1)**2.0_dp)*(theta(2)**2.0_dp))
        f(7) = (3.0_dp*theta(1)**2.0_dp+3.0_dp*theta(1)+1.0_dp)/ &
               ((theta(1)**3.0_dp)*(theta(2)**3.0_dp))
        f(8) = (4.0_dp*theta(1)**3.0_dp+6.0_dp*theta(1)**2.0_dp+4.0_dp*theta(1)+1.0_dp)/ &
               ((theta(1)**4.0_dp)*(theta(2)**4.0_dp))
    end function compute_f


    pure function compute_s(theta_2, f) result(s)
    !> Only need the 2nd element of theta
       real(dp), intent(in) :: theta_2
       real(dp), intent(in) :: f(8)
       real(dp) :: s(5:10, 1:7)

       s(5, 1) = -1.0_dp/(theta_2**4.0_dp)
       s(5, 2) = f(4)
       s(5, 3) = f(4)
       s(5, 4) = s(5, 1)
       s(5, 5) = -(f(1)/(theta_2**4.0_dp)+f(4)/theta_2)
       s(5, 6) = (f(4)/(theta_2**2.0_dp)-f(2)/(theta_2**4.0_dp))
       s(5, 7) = -(f(3)/(theta_2**4.0_dp)+f(4)/(theta_2**3.0_dp))

       s(6, 1) = -1.0_dp
       s(6, 2) = f(8)
       s(6, 3) = f(8)
       s(6, 4) = -1.0_dp
       s(6, 5) = (f(5)+f(8))
       s(6, 6) = (f(8)-f(6))
       s(6, 7) = (f(7)+f(8))

       s(7, 1) = -1.0_dp
       s(7, 2) = (1.0_dp+f(4))
       s(7, 3) = f(4)
       s(7, 4) = (f(4)-f(1))
       s(7, 5) = (f(4)-f(2))
       s(7, 6) = (f(4)-f(3))

       s(8, 1) =-1.0_dp/(theta_2**4.0_dp)
       s(8, 2) =(f(8)+(1.0_dp/theta_2**4.0_dp))
       s(8, 3) =f(8)
       s(8, 4) =(f(5)/(theta_2**4.0_dp)-f(8)/theta_2)
       s(8, 5) =(f(8)/(theta_2**2.0_dp)-f(6)/(theta_2**4.0_dp))
       s(8, 6) =(f(7)/(theta_2**4.0_dp)-f(8)/(theta_2**3.0_dp))

       s(9, 1) = -s(6, 6) * s(5, 1)
       s(9, 2) = s(5, 6) * s(6, 1)
       s(9, 3) = -(s(5, 4) * s(6, 6) + s(5, 6) * s(6, 2))
       s(9, 4) = (s(5, 6) * s(6, 4)+ s(6, 6) * s(5, 2))
       s(9, 5) = (s(5, 6) * s(6, 3)- s(6, 6) * s(5, 3))
       s(9, 6) = (s(5, 6) * s(6, 5)- s(6, 6) * s(5, 5))
       s(9, 7) = (s(5, 6) * s(6, 7)- s(6, 6) * s(5, 7))

       s(10, 1) = -s(8, 5) * s(7, 1)
       s(10, 2) = -s(8, 5) * s(7, 2)
       s(10, 3) = s(7, 5) * s(8, 1)
       s(10, 4) = s(7, 5) * s(8, 2)
       s(10, 5) = (s(7, 5) * s(8, 3) - s(8, 5) * s(7, 3))
       s(10, 6) = (s(7, 5) * s(8, 4) - s(8, 5) * s(7, 4))
       s(10, 7) = (s(7, 5) * s(8, 6) - s(8, 5) * s(7, 6))
    end function compute_s

    pure function compute_ak(theta) result(ak)
       real(dp), intent(in) :: theta(3)
       real(dp) :: ak(7)

       ak(1) = (1.0_dp+2.0_dp*theta(1))*(theta(3)+theta(3)**2.0_dp)*theta(2)
       ak(2) = (1.0_dp+theta(1))*(2.0_dp*theta(3)+theta(3)**2.0_dp)
       ak(3) = (1.0_dp+theta(1))
       ak(4) = (1.0_dp+theta(1))*((1.0_dp+theta(3))**2.0_dp)
       ak(5) = ((1.0_dp+theta(1))**2.0_dp)*(theta(3)+theta(3)**2.0_dp)*theta(2)
       ak(6) = (theta(1)**2.0_dp)*theta(2)*(theta(3)+theta(3)**2.0_dp)
       ak(7) = (1.0_dp+theta(1))*(theta(3)+theta(3)**2.0_dp)*theta(2)
    end function compute_ak
end module biocfd_navier_stokes

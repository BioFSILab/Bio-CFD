module biocfd_navier_stokes
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use global
  implicit none
  private

  public :: non_uni_coeff, nsmomentum2order

contains
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	SUBROUTINE non_uni_coeff
	!write(*,*)'entered non_uni_coeff'
	INTEGER  (dp) :: i, j, k, g, nx_var, ny_var, nz_var
	REAL (dp)   :: f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8,                             &
                         s51, s52, s53, s54, s55, s56, s57, s61, s62, s63, s64, s65, s66, s67,  &
                         s71, s72, s73, s74, s75, s76, s81, s82, s83, s84, s85, s86, s91, s92,  &
                         s93, s94, s95, s96, s97, s101, s102, s103, s104, s105, s106, s107,     &
                         ak_1, ak_2, ak_3, ak_4, ak_5, ak_6, ak_7, theta_1, theta_2, theta_3,   &
                         tmp_dx1, tmp_dx2, tmp_dx3, tmp_dx4, tmp_dy1, tmp_dy2, tmp_dy3,         &
                         tmp_dy4, tmp_dz1, tmp_dz2, tmp_dz3, tmp_dz4


        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz

        ALLOCATE(block(g)%ca1_uu(nx_var+2), block(g)%ca2_uu(nx_var+2), block(g)%ca3_uu(nx_var+2), &
             block(g)%ca4_uu(nx_var+2), block(g)%ca5_uu(nx_var+2), block(g)%ca6_uu(nx_var+2), &
             block(g)%ck1_uu(nx_var+2), block(g)%ck2_uu(nx_var+2), block(g)%ck3_uu(nx_var+2), &
             block(g)%ck4_uu(nx_var+2), block(g)%ck5_uu(nx_var+2), block(g)%ck6_uu(nx_var+2), &
             block(g)%ca1_vv(ny_var+2), block(g)%ca2_vv(ny_var+2), block(g)%ca3_vv(ny_var+2), &
             block(g)%ca4_vv(ny_var+2), block(g)%ca5_vv(ny_var+2), block(g)%ca6_vv(ny_var+2), &
             block(g)%ck1_vv(ny_var+2), block(g)%ck2_vv(ny_var+2), block(g)%ck3_vv(ny_var+2), &
             block(g)%ck4_vv(ny_var+2), block(g)%ck5_vv(ny_var+2), block(g)%ck6_vv(ny_var+2), &
             block(g)%ca1_ww(nz_var+2), block(g)%ca2_ww(nz_var+2), block(g)%ca3_ww(nz_var+2), &
             block(g)%ca4_ww(nz_var+2), block(g)%ca5_ww(nz_var+2), block(g)%ca6_ww(nz_var+2), &
             block(g)%ck1_ww(nz_var+2), block(g)%ck2_ww(nz_var+2), block(g)%ck3_ww(nz_var+2), &
             block(g)%ck4_ww(nz_var+2), block(g)%ck5_ww(nz_var+2), block(g)%ck6_ww(nz_var+2), &
             block(g)%ca1_uv(nx_var+2), block(g)%ca2_uv(nx_var+2), block(g)%ca3_uv(nx_var+2), &
             block(g)%ca4_uv(nx_var+2), block(g)%ca5_uv(nx_var+2), block(g)%ca6_uv(nx_var+2), &
             block(g)%ck1_uv(nx_var+2), block(g)%ck2_uv(nx_var+2), block(g)%ck3_uv(nx_var+2), &
             block(g)%ck4_uv(nx_var+2), block(g)%ck5_uv(nx_var+2), block(g)%ck6_uv(nx_var+2), &
             block(g)%ca1_uw(nx_var+2), block(g)%ca2_uw(nx_var+2), block(g)%ca3_uw(nx_var+2), &
             block(g)%ca4_uw(nx_var+2), block(g)%ca5_uw(nx_var+2), block(g)%ca6_uw(nx_var+2))
        ALLOCATE(block(g)%ck1_uw(nx_var+2), block(g)%ck2_uw(nx_var+2), block(g)%ck3_uw(nx_var+2), &
             block(g)%ck4_uw(nx_var+2), block(g)%ck5_uw(nx_var+2), block(g)%ck6_uw(nx_var+2), &
             block(g)%ca1_vu(ny_var+2), block(g)%ca2_vu(ny_var+2), block(g)%ca3_vu(ny_var+2), &
             block(g)%ca4_vu(ny_var+2), block(g)%ca5_vu(ny_var+2), block(g)%ca6_vu(ny_var+2), &
             block(g)%ck1_vu(ny_var+2), block(g)%ck2_vu(ny_var+2), block(g)%ck3_vu(ny_var+2), &
             block(g)%ck4_vu(ny_var+2), block(g)%ck5_vu(ny_var+2), block(g)%ck6_vu(ny_var+2), &
             block(g)%ca1_vw(ny_var+2), block(g)%ca2_vw(ny_var+2), block(g)%ca3_vw(ny_var+2), &
             block(g)%ca4_vw(ny_var+2), block(g)%ca5_vw(ny_var+2), block(g)%ca6_vw(ny_var+2), &
             block(g)%ck1_vw(ny_var+2), block(g)%ck2_vw(ny_var+2), block(g)%ck3_vw(ny_var+2), &
             block(g)%ck4_vw(ny_var+2), block(g)%ck5_vw(ny_var+2), block(g)%ck6_vw(ny_var+2), &
             block(g)%ca1_wu(nz_var+2), block(g)%ca2_wu(nz_var+2), block(g)%ca3_wu(nz_var+2), &
             block(g)%ca4_wu(nz_var+2), block(g)%ca5_wu(nz_var+2), block(g)%ca6_wu(nz_var+2), &
             block(g)%ck1_wu(nz_var+2), block(g)%ck2_wu(nz_var+2), block(g)%ck3_wu(nz_var+2), &
             block(g)%ck4_wu(nz_var+2), block(g)%ck5_wu(nz_var+2), block(g)%ck6_wu(nz_var+2), &
             block(g)%ca1_wv(nz_var+2), block(g)%ca2_wv(nz_var+2), block(g)%ca3_wv(nz_var+2), &
             block(g)%ca4_wv(nz_var+2), block(g)%ca5_wv(nz_var+2), block(g)%ca6_wv(nz_var+2), &
             block(g)%ck1_wv(nz_var+2), block(g)%ck2_wv(nz_var+2), block(g)%ck3_wv(nz_var+2), &
             block(g)%ck4_wv(nz_var+2), block(g)%ck5_wv(nz_var+2), block(g)%ck6_wv(nz_var+2))
	END DO

        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz
	do i=2,nx_var
	theta_1=block(g)%deltax(i)/block(g)%deltax(i-1)
	theta_2=block(g)%deltax(i+1)/block(g)%deltax(i)
	theta_3=block(g)%deltax(i+2)/block(g)%deltax(i+1)

        f_1=theta_3
        f_2=2.0_dp*theta_3+theta_3**2.0_dp
	f_3=3.0_dp*theta_3+3.0_dp*theta_3**2.0_dp+theta_3**3.0_dp
	f_4=4.0_dp*theta_3+6.0_dp*theta_3**2.0_dp+4.0_dp*theta_3**3.0_dp+theta_3**4.0_dp
	f_5=1.0_dp/(theta_1*theta_2)
	f_6=(2.0_dp*theta_1+1.0_dp)/((theta_1**2.0_dp)*(theta_2**2.0_dp))
	f_7=(3.0_dp*theta_1**2.0_dp+3.0_dp*theta_1+1.0_dp)/((theta_1**3.0_dp)*(theta_2**3.0_dp))
        f_8=(4.0_dp*theta_1**3.0_dp+6.0_dp*theta_1**2.0_dp+4.0_dp*theta_1+1.0_dp)/&
            ((theta_1**4.0_dp)*(theta_2**4.0_dp))

	s51=-1.0_dp/(theta_2**4.0_dp)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0_dp)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0_dp)-f_2/(theta_2**4.0_dp))
	s57=-(f_3/(theta_2**4.0_dp)+f_4/(theta_2**3.0_dp))

	s61=-1.0_dp
	s62=f_8
	s63=f_8
	s64=-1.0
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0_dp
	s72=(1.0_dp+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0_dp/(theta_2**4.0_dp)
	s82=(f_8+(1.0_dp/theta_2**4.0_dp))
	s83=f_8
	s84=(f_5/(theta_2**4.0_dp)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0_dp)-f_6/(theta_2**4.0_dp))
	s86=(f_7/(theta_2**4.0_dp)-f_8/(theta_2**3.0_dp))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	block(g)%ca1_uu(i)=(s97*s101-s107*s91)
	block(g)%ca2_uu(i)=(s97*s102+s107*s93)
	block(g)%ca3_uu(i)=(s95*s107-s105*s97)
	block(g)%ca4_uu(i)=(s94*s107+s104*s97)
	block(g)%ca5_uu(i)=(s97*s103-s107*s92)
	block(g)%ca6_uu(i)=(s97*s106-s107*s96)

	ak_1=(1.0_dp+2.0_dp*theta_1)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_2=(1.0_dp+theta_1)*(2.0_dp*theta_3+theta_3**2.0_dp)
	ak_3=(1.0_dp+theta_1)
	ak_4=(1.0_dp+theta_1)*((1.0_dp+theta_3)**2.0_dp)
	ak_5=((1.0_dp+theta_1)**2.0_dp)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_6=(theta_1**2.0_dp)*theta_2*(theta_3+theta_3**2.0_dp)
	ak_7=(1.0_dp+theta_1)*(theta_3+theta_3**2.0_dp)*theta_2

	block(g)%ck1_uu(i)=ak_3
	block(g)%ck2_uu(i)=-ak_4
	block(g)%ck3_uu(i)=(ak_1+ak_2)
	block(g)%ck4_uu(i)=-ak_5
	block(g)%ck5_uu(i)=ak_6
	block(g)%ck6_uu(i)=ak_7
	enddo
        ENDDO
        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz
	do j=2,ny_var
	theta_1=block(g)%deltay(j)/block(g)%deltay(j-1)
	theta_2=block(g)%deltay(j+1)/block(g)%deltay(j)
	theta_3=block(g)%deltay(j+2)/block(g)%deltay(j+1)

	f_1=theta_3
	f_2=2.0_dp*theta_3+theta_3**2.0_dp
	f_3=3.0_dp*theta_3+3.0_dp*theta_3**2.0_dp+theta_3**3.0_dp
	f_4=4.0_dp*theta_3+6.0_dp*theta_3**2.0_dp+4.0_dp*theta_3**3.0_dp+theta_3**4.0_dp
	f_5=1.0_dp/(theta_1*theta_2)
	f_6=(2.0_dp*theta_1+1.0_dp)/((theta_1**2.0_dp)*(theta_2**2.0_dp))
	f_7=(3.0_dp*theta_1**2.0_dp+3.0_dp*theta_1+1.0_dp)/((theta_1**3.0_dp)*(theta_2**3.0_dp))
        f_8=(4.0_dp*theta_1**3.0_dp+6.0_dp*theta_1**2.0_dp+4.0_dp*theta_1+1.0_dp)/&
            ((theta_1**4.0_dp)*(theta_2**4.0_dp))

	s51=-1.0_dp/(theta_2**4.0_dp)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0_dp)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0_dp)-f_2/(theta_2**4.0_dp))
	s57=-(f_3/(theta_2**4.0_dp)+f_4/(theta_2**3.0_dp))

	s61=-1.0_dp
	s62=f_8
	s63=f_8
	s64=-1.0_dp
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0_dp
	s72=(1.0_dp+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0_dp/(theta_2**4.0_dp)
	s82=(f_8+(1.0_dp/theta_2**4.0_dp))
	s83=f_8
	s84=(f_5/(theta_2**4.0_dp)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0_dp)-f_6/(theta_2**4.0_dp))
	s86=(f_7/(theta_2**4.0_dp)-f_8/(theta_2**3.0_dp))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	block(g)%ca1_vv(j)=(s97*s101-s107*s91)
	block(g)%ca2_vv(j)=(s97*s102+s107*s93)
	block(g)%ca3_vv(j)=(s95*s107-s105*s97)
	block(g)%ca4_vv(j)=(s94*s107+s104*s97)
	block(g)%ca5_vv(j)=(s97*s103-s107*s92)
	block(g)%ca6_vv(j)=(s97*s106-s107*s96)

	ak_1=(1.0_dp+2.0_dp*theta_1)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_2=(1.0_dp+theta_1)*(2.0_dp*theta_3+theta_3**2.0_dp)
	ak_3=(1.0_dp+theta_1)
	ak_4=(1.0_dp+theta_1)*((1.0_dp+theta_3)**2.0_dp)
	ak_5=((1.0_dp+theta_1)**2.0_dp)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_6=(theta_1**2.0_dp)*theta_2*(theta_3+theta_3**2.0_dp)
	ak_7=(1.0_dp+theta_1)*(theta_3+theta_3**2.0_dp)*theta_2

	block(g)%ck1_vv(j)=ak_3
	block(g)%ck2_vv(j)=-ak_4
	block(g)%ck3_vv(j)=(ak_1+ak_2)
	block(g)%ck4_vv(j)=-ak_5
	block(g)%ck5_vv(j)=ak_6
	block(g)%ck6_vv(j)=ak_7
	enddo
        ENDDO

        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz
	do k=2,nz_var
	theta_1=block(g)%deltaz(k)/block(g)%deltaz(k-1)
	theta_2=block(g)%deltaz(k+1)/block(g)%deltaz(k)
	theta_3=block(g)%deltaz(k+2)/block(g)%deltaz(k+1)

	f_1=theta_3
	f_2=2.0_dp*theta_3+theta_3**2.0_dp
	f_3=3.0_dp*theta_3+3.0_dp*theta_3**2.0_dp+theta_3**3.0_dp
	f_4=4.0_dp*theta_3+6.0_dp*theta_3**2.0_dp+4.0_dp*theta_3**3.0_dp+theta_3**4.0_dp
	f_5=1.0_dp/(theta_1*theta_2)
	f_6=(2.0_dp*theta_1+1.0_dp)/((theta_1**2.0_dp)*(theta_2**2.0_dp))
	f_7=(3.0_dp*theta_1**2.0_dp+3.0_dp*theta_1+1.0_dp)/((theta_1**3.0_dp)*(theta_2**3.0_dp))
        f_8=(4.0_dp*theta_1**3.0_dp+6.0_dp*theta_1**2.0_dp+4.0_dp*theta_1+1.0_dp)/&
            ((theta_1**4.0_dp)*(theta_2**4.0_dp))

	s51=-1.0_dp/(theta_2**4.0_dp)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0_dp)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0_dp)-f_2/(theta_2**4.0_dp))
	s57=-(f_3/(theta_2**4.0_dp)+f_4/(theta_2**3.0_dp))

	s61=-1.0_dp
	s62=f_8
	s63=f_8
	s64=-1.0_dp
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0_dp
	s72=(1.0_dp+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0_dp/(theta_2**4.0_dp)
	s82=(f_8+(1.0_dp/theta_2**4.0_dp))
	s83=f_8
	s84=(f_5/(theta_2**4.0_dp)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0_dp)-f_6/(theta_2**4.0_dp))
	s86=(f_7/(theta_2**4.0_dp)-f_8/(theta_2**3.0_dp))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	block(g)%ca1_ww(k)=(s97*s101-s107*s91)
	block(g)%ca2_ww(k)=(s97*s102+s107*s93)
	block(g)%ca3_ww(k)=(s95*s107-s105*s97)
	block(g)%ca4_ww(k)=(s94*s107+s104*s97)
	block(g)%ca5_ww(k)=(s97*s103-s107*s92)
	block(g)%ca6_ww(k)=(s97*s106-s107*s96)

	ak_1=(1.0_dp+2.0_dp*theta_1)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_2=(1.0_dp+theta_1)*(2.0_dp*theta_3+theta_3**2.0_dp)
	ak_3=(1.0_dp+theta_1)
	ak_4=(1.0_dp+theta_1)*((1.0_dp+theta_3)**2.0_dp)
	ak_5=((1.0_dp+theta_1)**2.0_dp)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_6=(theta_1**2.0_dp)*theta_2*(theta_3+theta_3**2.0_dp)
	ak_7=(1.0_dp+theta_1)*(theta_3+theta_3**2.0_dp)*theta_2

	block(g)%ck1_ww(k)=ak_3
	block(g)%ck2_ww(k)=-ak_4
	block(g)%ck3_ww(k)=(ak_1+ak_2)
	block(g)%ck4_ww(k)=-ak_5
	block(g)%ck5_ww(k)=ak_6
	block(g)%ck6_ww(k)=ak_7
	enddo
        END DO

        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz
	do i=2,nx_var
	if(i==2)then
	tmp_dx1=block(g)%deltax(i-1)
	else
	tmp_dx1=0.5_dp*(block(g)%deltax(i-1)+block(g)%deltax(i-2))
	endif
	tmp_dx2=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i-1))
	tmp_dx3=0.5_dp*(block(g)%deltax(i)+block(g)%deltax(i+1))
	tmp_dx4=0.5_dp*(block(g)%deltax(i+1)+block(g)%deltax(i+2))

	theta_1=tmp_dx2/tmp_dx1
	theta_2=tmp_dx3/tmp_dx2
	theta_3=tmp_dx4/tmp_dx3

	f_1=theta_3
	f_2=2.0_dp*theta_3+theta_3**2.0_dp
	f_3=3.0_dp*theta_3+3.0_dp*theta_3**2.0_dp+theta_3**3.0_dp
	f_4=4.0_dp*theta_3+6.0_dp*theta_3**2.0_dp+4.0_dp*theta_3**3.0_dp+theta_3**4.0_dp
	f_5=1.0_dp/(theta_1*theta_2)
	f_6=(2.0_dp*theta_1+1.0_dp)/((theta_1**2.0_dp)*(theta_2**2.0_dp))
	f_7=(3.0_dp*theta_1**2.0_dp+3.0_dp*theta_1+1.0_dp)/((theta_1**3.0_dp)*(theta_2**3.0_dp))
        f_8=(4.0_dp*theta_1**3.0_dp+6.0_dp*theta_1**2.0_dp+4.0_dp*theta_1+1.0_dp)/&
            ((theta_1**4.0_dp)*(theta_2**4.0_dp))

	s51=-1.0_dp/(theta_2**4.0_dp)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0_dp)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0_dp)-f_2/(theta_2**4.0_dp))
	s57=-(f_3/(theta_2**4.0_dp)+f_4/(theta_2**3.0_dp))

	s61=-1.0_dp
	s62=f_8
	s63=f_8
	s64=-1.0_dp
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0_dp
	s72=(1.0_dp+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0_dp/(theta_2**4.0_dp)
	s82=(f_8+(1.0_dp/theta_2**4.0_dp))
	s83=f_8
	s84=(f_5/(theta_2**4.0_dp)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0_dp)-f_6/(theta_2**4.0_dp))
	s86=(f_7/(theta_2**4.0_dp)-f_8/(theta_2**3.0_dp))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	block(g)%ca1_uv(i)=(s97*s101-s107*s91)
	block(g)%ca1_uw(i)=block(g)%ca1_uv(i)
	block(g)%ca2_uv(i)=(s97*s102+s107*s93)
	block(g)%ca2_uw(i)=block(g)%ca2_uv(i)
	block(g)%ca3_uv(i)=(s95*s107-s105*s97)
	block(g)%ca3_uw(i)=block(g)%ca3_uv(i)
	block(g)%ca4_uv(i)=(s94*s107+s104*s97)
	block(g)%ca4_uw(i)=block(g)%ca4_uv(i)
	block(g)%ca5_uv(i)=(s97*s103-s107*s92)
	block(g)%ca5_uw(i)=block(g)%ca5_uv(i)
	block(g)%ca6_uv(i)=(s97*s106-s107*s96)
	block(g)%ca6_uw(i)=block(g)%ca6_uv(i)

	ak_1=(1.0_dp+2.0_dp*theta_1)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_2=(1.0_dp+theta_1)*(2.0_dp*theta_3+theta_3**2.0_dp)
	ak_3=(1.0_dp+theta_1)
	ak_4=(1.0_dp+theta_1)*((1.0_dp+theta_3)**2.0_dp)
	ak_5=((1.0_dp+theta_1)**2.0_dp)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_6=(theta_1**2.0_dp)*theta_2*(theta_3+theta_3**2.0_dp)
	ak_7=(1.0_dp+theta_1)*(theta_3+theta_3**2.0_dp)*theta_2

	block(g)%ck1_uv(i)=ak_3
	block(g)%ck1_uw(i)=block(g)%ck1_uv(i)
	block(g)%ck2_uv(i)=-ak_4
	block(g)%ck2_uw(i)=block(g)%ck2_uv(i)
	block(g)%ck3_uv(i)=(ak_1+ak_2)
	block(g)%ck3_uw(i)=block(g)%ck3_uv(i)
	block(g)%ck4_uv(i)=-ak_5
	block(g)%ck4_uw(i)=block(g)%ck4_uv(i)
	block(g)%ck5_uv(i)=ak_6
	block(g)%ck5_uw(i)=block(g)%ck5_uv(i)
	block(g)%ck6_uv(i)=ak_7
	block(g)%ck6_uw(i)=block(g)%ck6_uv(i)
	enddo
        END DO

        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz
	do j=2,ny_var
	if(j==2)then
	tmp_dy1=block(g)%deltay(j-1)
	else
	tmp_dy1=0.5_dp*(block(g)%deltay(j-1)+block(g)%deltay(j-2))
	endif
	tmp_dy2=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j-1))
	tmp_dy3=0.5_dp*(block(g)%deltay(j)+block(g)%deltay(j+1))
	tmp_dy4=0.5_dp*(block(g)%deltay(j+1)+block(g)%deltay(j+2))

	theta_1=tmp_dy2/tmp_dy1
	theta_2=tmp_dy3/tmp_dy2
	theta_3=tmp_dy4/tmp_dy3

	f_1=theta_3
	f_2=2.0_dp*theta_3+theta_3**2.0_dp
	f_3=3.0_dp*theta_3+3.0_dp*theta_3**2.0_dp+theta_3**3.0_dp
	f_4=4.0_dp*theta_3+6.0_dp*theta_3**2.0_dp+4.0_dp*theta_3**3.0_dp+theta_3**4.0_dp
	f_5=1.0_dp/(theta_1*theta_2)
	f_6=(2.0_dp*theta_1+1.0_dp)/((theta_1**2.0_dp)*(theta_2**2.0_dp))
	f_7=(3.0_dp*theta_1**2.0_dp+3.0_dp*theta_1+1.0_dp)/((theta_1**3.0_dp)*(theta_2**3.0_dp))
        f_8=(4.0_dp*theta_1**3.0_dp+6.0_dp*theta_1**2.0_dp+4.0_dp*theta_1+1.0_dp)/&
            ((theta_1**4.0_dp)*(theta_2**4.0_dp))

	s51=-1.0_dp/(theta_2**4.0_dp)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0_dp)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0_dp)-f_2/(theta_2**4.0_dp))
	s57=-(f_3/(theta_2**4.0_dp)+f_4/(theta_2**3.0_dp))

	s61=-1.0_dp
	s62=f_8
	s63=f_8
	s64=-1.0_dp
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0_dp
	s72=(1.0_dp+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0_dp/(theta_2**4.0_dp)
	s82=(f_8+(1.0_dp/theta_2**4.0_dp))
	s83=f_8
	s84=(f_5/(theta_2**4.0_dp)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0_dp)-f_6/(theta_2**4.0_dp))
	s86=(f_7/(theta_2**4.0_dp)-f_8/(theta_2**3.0_dp))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	block(g)%ca1_vu(j)=(s97*s101-s107*s91)
	block(g)%ca1_vw(j)=block(g)%ca1_vu(j)
	block(g)%ca2_vu(j)=(s97*s102+s107*s93)
	block(g)%ca2_vw(j)=block(g)%ca2_vu(j)
	block(g)%ca3_vu(j)=(s95*s107-s105*s97)
	block(g)%ca3_vw(j)=block(g)%ca3_vu(j)
	block(g)%ca4_vu(j)=(s94*s107+s104*s97)
	block(g)%ca4_vw(j)=block(g)%ca4_vu(j)
	block(g)%ca5_vu(j)=(s97*s103-s107*s92)
	block(g)%ca5_vw(j)=block(g)%ca5_vu(j)
	block(g)%ca6_vu(j)=(s97*s106-s107*s96)
	block(g)%ca6_vw(j)=block(g)%ca6_vu(j)

	ak_1=(1.0_dp+2.0_dp*theta_1)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_2=(1.0_dp+theta_1)*(2.0_dp*theta_3+theta_3**2.0_dp)
	ak_3=(1.0_dp+theta_1)
	ak_4=(1.0_dp+theta_1)*((1.0_dp+theta_3)**2.0_dp)
	ak_5=((1.0_dp+theta_1)**2.0_dp)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_6=(theta_1**2.0_dp)*theta_2*(theta_3+theta_3**2.0_dp)
	ak_7=(1.0_dp+theta_1)*(theta_3+theta_3**2.0_dp)*theta_2

	block(g)%ck1_vu(j)=ak_3
	block(g)%ck1_vw(j)=block(g)%ck1_vu(j)
	block(g)%ck2_vu(j)=-ak_4
	block(g)%ck2_vw(j)=block(g)%ck2_vu(j)
	block(g)%ck3_vu(j)=(ak_1+ak_2)
	block(g)%ck3_vw(j)=block(g)%ck3_vu(j)
	block(g)%ck4_vu(j)=-ak_5
	block(g)%ck4_vw(j)=block(g)%ck4_vu(j)
	block(g)%ck5_vu(j)=ak_6
	block(g)%ck5_vw(j)=block(g)%ck5_vu(j)
	block(g)%ck6_vu(j)=ak_7
	block(g)%ck6_vw(j)=block(g)%ck6_vu(j)
	enddo
        ENDDO

        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz
	do k=2,nz_var
	if(k==2)then
	tmp_dz1=block(g)%deltaz(k-1)
	else
	tmp_dz1=0.5_dp*(block(g)%deltaz(k-1)+block(g)%deltaz(k-2))
	endif
	tmp_dz2=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
	tmp_dz3=0.5_dp*(block(g)%deltaz(k)+block(g)%deltaz(k+1))
	tmp_dz4=0.5_dp*(block(g)%deltaz(k+1)+block(g)%deltaz(k+2))

	theta_1=tmp_dz2/tmp_dz1
	theta_2=tmp_dz3/tmp_dz2
	theta_3=tmp_dz4/tmp_dz3

	f_1=theta_3
	f_2=2.0_dp*theta_3+theta_3**2.0_dp
	f_3=3.0_dp*theta_3+3.0_dp*theta_3**2.0_dp+theta_3**3.0_dp
	f_4=4.0_dp*theta_3+6.0_dp*theta_3**2.0_dp+4.0_dp*theta_3**3.0_dp+theta_3**4.0_dp
	f_5=1.0_dp/(theta_1*theta_2)
	f_6=(2.0_dp*theta_1+1.0_dp)/((theta_1**2.0_dp)*(theta_2**2.0_dp))
	f_7=(3.0_dp*theta_1**2.0_dp+3.0_dp*theta_1+1.0_dp)/((theta_1**3.0_dp)*(theta_2**3.0_dp))
        f_8=(4.0_dp*theta_1**3.0_dp+6.0_dp*theta_1**2.0_dp+4.0_dp*theta_1+1.0_dp)/&
            ((theta_1**4.0_dp)*(theta_2**4.0_dp))

	s51=-1.0_dp/(theta_2**4.0_dp)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0_dp)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0_dp)-f_2/(theta_2**4.0_dp))
	s57=-(f_3/(theta_2**4.0_dp)+f_4/(theta_2**3.0_dp))

	s61=-1.0_dp
	s62=f_8
	s63=f_8
	s64=-1.0_dp
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0_dp
	s72=(1.0_dp+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0_dp/(theta_2**4.0_dp)
	s82=(f_8+(1.0_dp/theta_2**4.0_dp))
	s83=f_8
	s84=(f_5/(theta_2**4.0_dp)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0_dp)-f_6/(theta_2**4.0_dp))
	s86=(f_7/(theta_2**4.0_dp)-f_8/(theta_2**3.0_dp))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	block(g)%ca1_wu(k)=(s97*s101-s107*s91)
	block(g)%ca1_wv(k)=block(g)%ca1_wu(k)
	block(g)%ca2_wu(k)=(s97*s102+s107*s93)
	block(g)%ca2_wv(k)=block(g)%ca2_wu(k)
	block(g)%ca3_wu(k)=(s95*s107-s105*s97)
	block(g)%ca3_wv(k)=block(g)%ca3_wu(k)
	block(g)%ca4_wu(k)=(s94*s107+s104*s97)
	block(g)%ca4_wv(k)=block(g)%ca4_wu(k)
	block(g)%ca5_wu(k)=(s97*s103-s107*s92)
	block(g)%ca5_wv(k)=block(g)%ca5_wu(k)
	block(g)%ca6_wu(k)=(s97*s106-s107*s96)
	block(g)%ca6_wv(k)=block(g)%ca6_wu(k)

	ak_1=(1.0_dp+2.0_dp*theta_1)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_2=(1.0_dp+theta_1)*(2.0_dp*theta_3+theta_3**2.0_dp)
	ak_3=(1.0_dp+theta_1)
	ak_4=(1.0_dp+theta_1)*((1.0_dp+theta_3)**2.0_dp)
	ak_5=((1.0_dp+theta_1)**2.0_dp)*(theta_3+theta_3**2.0_dp)*theta_2
	ak_6=(theta_1**2.0_dp)*theta_2*(theta_3+theta_3**2.0_dp)
	ak_7=(1.0_dp+theta_1)*(theta_3+theta_3**2.0_dp)*theta_2

	block(g)%ck1_wu(k)=ak_3
	block(g)%ck1_wv(k)=block(g)%ck1_wu(k)
	block(g)%ck2_wu(k)=-ak_4
	block(g)%ck2_wv(k)=block(g)%ck2_wu(k)
	block(g)%ck3_wu(k)=(ak_1+ak_2)
	block(g)%ck3_wv(k)=block(g)%ck3_wu(k)
	block(g)%ck4_wu(k)=-ak_5
	block(g)%ck4_wv(k)=block(g)%ck4_wu(k)
	block(g)%ck5_wu(k)=ak_6
	block(g)%ck5_wv(k)=block(g)%ck5_wu(k)
	block(g)%ck6_wu(k)=ak_7
	block(g)%ck6_wv(k)=block(g)%ck6_wu(k)
	enddo
        END DO
       !!$acc update device (                                 &
       !!$acc ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
	!!$acc ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
	!!$acc ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
	!!$acc ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
	!!$acc ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
	!!$acc ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww, &
	!!$acc ca1_uv, ca2_uv, ca3_uv, ca4_uv, ca5_uv, ca6_uv, &
	!!$acc ck1_uv, ck2_uv, ck3_uv, ck4_uv, ck5_uv, ck6_uv, &
	!!$acc ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
	!!$acc ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
	!!$acc ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
	!!$acc ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
	!!$acc ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
	!!$acc ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
	!!$acc ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
	!!$acc ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
	!!$acc ca1_wv, ca2_wv, ca3_wv, ca4_wv, ca5_wv, ca6_wv, &
	!!$acc ck1_wv, ck2_wv, ck3_wv, ck4_wv, ck5_wv, ck6_wv)

	write(*,*)'leaving non_uni_coeff'

	END SUBROUTINE non_uni_coeff

!csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine nsMomentum2order
!c***********************************************************************
!c     navier-stokes equations for constant properties
!c
!c***********************************************************************
!c***********************************************************************
      !write(6,*)'has entered nseqcp'
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



!!$omp parallel do &
!!$omp private (i, j, k, n1, dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,   &
!!$omp	        u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
!!$omp		 w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
!!$omp	        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
!!$omp	         dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          &
!!$omp		 r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
!!$omp	        dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
!!$omp	        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,            &
!!$omp	        dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
!!$omp	        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
!!$omp		 duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw, &
!!$omp           index_ip1,index_im1,index_jp1,index_jm1,index_kp1,index_km1) &
!!$omp           i11,j11,k11,temp_u2dotn,temp_u1dotn) &
!!$omp           temp_v2dotn,temp_v1dotn) &
!!$omp           temp_w2dotn,temp_w1dotn, temp_pdotn,g) &
!!$omp firstprivate(nx_var, ny_var, nz_var, rev, deltat, al) num_threads(3)
        DO g=1,nblocks
        nx_var=block(g)%nx
        ny_var=block(g)%ny
        nz_var=block(g)%nz
!!$acc parallel loop present(cell, cell2,  block(g)%TSIndexPtr, block(g)%index_ts, &
!!$acc        block(g)%fluidIndexPtr, block(g)%xcent, block(g)%ycent, block(g)%zcent, &
!!$acc        zp, z1, y1, xp, x1, yp) &
!!$acc present(block(g)%cosAlpha, block(g)%cosBeta, block(g)%cosGamma, block(g)%nelu2,&
!!$acc block(g)%nelu1, block(g)%nelv2, block(g)%nelv1, block(g)%nelw2, block(g)%nelw1,
!!$accc block(g)%nelp) &
!!$acc present(block(g)%p_ghost, block(g)%pt_ghost, block(g)%u2_ghost, &
!!$acc block(g)%u2t_ghost, block(g)%v2_ghost, block(g)%v2t_ghost, block(g)%w2_ghost, &
!!$acc block(g)%w2t_ghost, block(g)%u1_ghost, block(g)%u1t_ghost, block(g)%v1_ghost, &
!!$acc block(g)%v1t_ghost, block(g)%w1_ghost, block(g)%w1t_ghost, ut, vt, wt, u, v, w, p, &
!!$acc deltax, deltay, deltaz, resi_u, resi_v, resi_w) &
!!$acc present(ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, ck1_uu, ck2_uu, &
!!$acc ck3_uu, ck4_uu, ck5_uu, ck6_uu, ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
!!$acc ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, ca1_ww, ca2_ww, ca3_ww, ca4_ww, &
!!$acc ca5_ww, ca6_ww, ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww, ca1_uv, ca2_uv, ca3_uv, &
!!$acc ca4_uv, ca5_uv, ca6_uv, ck1_uv, ck2_uv, ck3_uv, ck4_uv, ck5_uv, ck6_uv, ca1_uw, ca2_uw, &
!!$acc ca3_uw, ca4_uw, ca5_uw, ca6_uw, ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, ca1_vu, &
!!$acc ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
!!$acc ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, &
!!$acc ck6_vw, ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, ck1_wu, ck2_wu, ck3_wu, ck4_wu, &
!!$acc ck5_wu, ck6_wu, ca1_wv, ca2_wv, ca3_wv, ca4_wv, ca5_wv, ca6_wv, ck1_wv, ck2_wv, ck3_wv, &
!!$acc ck4_wv, ck5_wv, ck6_wv)



!$acc parallel loop gang vector  &
!$acc private (i, j, k, n1, dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,   &
!$acc	        u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
!$acc		 w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
!$acc	        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
!$acc		 dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          &
!$acc		 r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
!$acc	        dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
!$acc	        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,            &
!$acc	        dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
!$acc	        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
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
      !GOTO 11
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
                   !temp_u1dotn = (block(g)%x1(i+1) - &
                   !block(g)%xcent(block(g)%nelu1(block(g)%index_ts(index_ip1))))*&
                   !block(g)%cosAlpha(block(g)%nelu1(block(g)%index_ts(index_ip1))) +  &
                !(block(g)%yp(j)   - block(g)%ycent(block(g)%nelu1(block(g)%index_ts(index_ip1))))*&
                   !block(g)%cosBeta(block(g)%nelu1(block(g)%index_ts(index_ip1)))  +  &
                 !(block(g)%zp(k)   - block(g)%zcent(block(g)%nelu1(block(g)%index_ts(index_ip1))))&
                   !*block(g)%cosGamma(block(g)%nelu1(block(g)%index_ts(index_ip1)))
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
                   !ELSE
                   !   block(g)%u(i+1,j,k) = block(g)%u2t_ghost(index_ip1)
                   ENDIF
                   !IF (temp_u1dotn.LT.0) THEN
                      !block(g)%u(i,j,k) = block(g)%u1_ghost(index_ip1)
                   !ELSE
                      !block(g)%u(i,j,k) = block(g)%u1t_ghost(index_ip1)
                   !ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i+1,j,k) = block(g)%v2_ghost(index_ip1)
                   !ELSE
                   !   block(g)%v(i+1,j,k) = block(g)%v2t_ghost(index_ip1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i+1,j-1,k) = block(g)%v1_ghost(index_ip1)
                   !ELSE
                      !block(g)%v(i+1,j-1,k) = block(g)%v1t_ghost(index_ip1)
                   ENDIF
                   IF (temp_w2dotn<0) THEN
                      block(g)%w(i+1,j,k) = block(g)%w2_ghost(index_ip1)
                   !ELSE
                      !block(g)%w(i+1,j,k) = block(g)%w2t_ghost(index_ip1)
                   ENDIF
                   IF (temp_w1dotn<0) THEN
                      block(g)%w(i+1,j,k-1) = block(g)%w1_ghost(index_ip1)
                   !ELSE
                      !block(g)%w(i+1,j,k-1) = block(g)%w1t_ghost(index_ip1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i+1,j,k) = block(g)%p_ghost(index_ip1)
                   !ELSE
                      !block(g)%p(i+1,j,k) = block(g)%pt_ghost(index_ip1)
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
                   !temp_u2dotn = (block(g)%x1(i)   - &
                   !block(g)%xcent(block(g)%nelu2(block(g)%index_ts(index_im1))))*&
                   !block(g)%cosAlpha(block(g)%nelu2(block(g)%index_ts(index_im1))) +  &
                   !(block(g)%yp(j)   - &
                   !block(g)%ycent(block(g)%nelu2(block(g)%index_ts(index_im1))))*&
                   !block(g)%cosBeta(block(g)%nelu2(block(g)%index_ts(index_im1)))  +  &
		   !(block(g)%zp(k)   - &
                   !block(g)%zcent(block(g)%nelu2(block(g)%index_ts(index_im1))))*&
                   !block(g)%cosGamma(block(g)%nelu2(block(g)%index_ts(index_im1)))
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
                   !IF (temp_u2dotn.LT.0) THEN
                      !block(g)%u(i-1,j,k) = block(g)%u2_ghost(index_im1)
                   !ELSE
                      !block(g)%u(i-1,j,k) = block(g)%u2t_ghost(index_im1)
                   !ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-2,j,k) = block(g)%u1_ghost(index_im1)
                   !ELSE
                      !block(g)%u(i-2,j,k) = block(g)%u1t_ghost(index_im1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i-1,j,k) = block(g)%v2_ghost(index_im1)
                   !ELSE
                      !block(g)%v(i-1,j,k) = block(g)%v2t_ghost(index_im1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i-1,j-1,k) = block(g)%v1_ghost(index_im1)
                   !ELSE
                      !block(g)%v(i-1,j-1,k) = block(g)%v1t_ghost(index_im1)
                   ENDIF
		             IF (temp_w2dotn<0) THEN
                      block(g)%w(i-1,j,k) = block(g)%w2_ghost(index_im1)
                   !ELSE
                      !block(g)%w(i-1,j,k) = block(g)%w2t_ghost(index_im1)
                   ENDIF
		             IF (temp_w1dotn<0) THEN
                      block(g)%w(i-1,j,k-1) = block(g)%w1_ghost(index_im1)
                   !ELSE
                      !block(g)%w(i-1,j,k-1) = block(g)%w1t_ghost(index_im1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i-1,j,k) = block(g)%p_ghost(index_im1)
                   !ELSE
                      !block(g)%p(i-1,j,k) = block(g)%pt_ghost(index_im1)
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
                   !temp_v1dotn = (block(g)%xp(i)   - &
                   !block(g)%xcent(block(g)%nelv1(block(g)%index_ts(index_jp1))))*&
                   !block(g)%cosAlpha(block(g)%nelv1(block(g)%index_ts(index_jp1))) +  &
                !(block(g)%y1(j+1) - block(g)%ycent(block(g)%nelv1(block(g)%index_ts(index_jp1))))*&
                   !block(g)%cosBeta(block(g)%nelv1(block(g)%index_ts(index_jp1)))  +  &
                !(block(g)%zp(k)   - block(g)%zcent(block(g)%nelv1(block(g)%index_ts(index_jp1))))*&
                   !block(g)%cosGamma(block(g)%nelv1(block(g)%index_ts(index_jp1)))
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
                   !ELSE
                      !block(g)%u(i,j+1,k) = block(g)%u2t_ghost(index_jp1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j+1,k) = block(g)%u1_ghost(index_jp1)
                   !ELSE
                      !block(g)%u(i-1,j+1,k) = block(g)%u1t_ghost(index_jp1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i,j+1,k) = block(g)%v2_ghost(index_jp1)
                   !ELSE
                      !block(g)%v(i,j+1,k) = block(g)%v2t_ghost(index_jp1)
                   ENDIF
                   !IF (temp_v1dotn.LT.0) THEN
                   !   block(g)%v(i,j,k) = block(g)%v1_ghost(index_jp1)
                   !ELSE
                   !   block(g)%v(i,j,k) = block(g)%v1t_ghost(index_jp1)
                   !ENDIF
		             IF (temp_w2dotn<0) THEN
                      block(g)%w(i,j+1,k) = block(g)%w2_ghost(index_jp1)
                   !ELSE
                      !block(g)%w(i,j+1,k) = block(g)%w2t_ghost(index_jp1)
                   ENDIF
		             IF (temp_w1dotn<0) THEN
                      block(g)%w(i,j+1,k-1) = block(g)%w1_ghost(index_jp1)
                   !ELSE
                      !block(g)%w(i,j+1,k-1) = block(g)%w1t_ghost(index_jp1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j+1,k) = block(g)%p_ghost(index_jp1)
                   !ELSE
                      !block(g)%p(i,j+1,k) = block(g)%pt_ghost(index_jp1)
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
                   !temp_v2dotn = (block(g)%xp(i)   - &
                   !block(g)%xcent(block(g)%nelv2(block(g)%index_ts(index_jm1))))*&
                   !block(g)%cosAlpha(block(g)%nelv2(block(g)%index_ts(index_jm1))) +  &
               ! (block(g)%y1(j)   - block(g)%ycent(block(g)%nelv2(block(g)%index_ts(index_jm1))))*&
                   !block(g)%cosBeta(block(g)%nelv2(block(g)%index_ts(index_jm1)))  +  &
                !(block(g)%zp(k)   - block(g)%zcent(block(g)%nelv2(block(g)%index_ts(index_jm1))))*&
                   !block(g)%cosGamma(block(g)%nelv2(block(g)%index_ts(index_jm1)))
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
                   !ELSE
                      !block(g)%u(i,j-1,k) = block(g)%u2t_ghost(index_jm1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j-1,k) = block(g)%u1_ghost(index_jm1)
                   !ELSE
                      !block(g)%u(i-1,j-1,k) = block(g)%u1t_ghost(index_jm1)
                   ENDIF
                   !IF (temp_v2dotn.LT.0) THEN
                   !   block(g)%v(i,j-1,k) = block(g)%v2_ghost(index_jm1)
                   !ELSE
                   !   block(g)%v(i,j-1,k) = block(g)%v2t_ghost(index_jm1)
                   !ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i,j-2,k) = block(g)%v1_ghost(index_jm1)
                   !ELSE
                      !block(g)%v(i,j-2,k) = block(g)%v1t_ghost(index_jm1)
                   ENDIF
                   IF (temp_w2dotn<0) THEN
                      block(g)%w(i,j-1,k) = block(g)%w2_ghost(index_jm1)
                   !ELSE
                      !block(g)%w(i,j-1,k) = block(g)%w2t_ghost(index_jm1)
                   ENDIF
                   IF (temp_w1dotn<0) THEN
                      block(g)%w(i,j-1,k-1) = block(g)%w1_ghost(index_jm1)
                   !ELSE
                      !block(g)%w(i,j-1,k-1) = block(g)%w1t_ghost(index_jm1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j-1,k) = block(g)%p_ghost(index_jm1)
                   !ELSE
                      !block(g)%p(i,j-1,k) = block(g)%pt_ghost(index_jm1)
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
                     !temp_w1dotn = (block(g)%xp(i)   - &
                     !block(g)%xcent(block(g)%nelw1(block(g)%index_ts(index_kp1))))*&
                     !block(g)%cosAlpha(block(g)%nelw1(block(g)%index_ts(index_kp1))) +  &
                !(block(g)%yp(j)   - block(g)%ycent(block(g)%nelw1(block(g)%index_ts(index_kp1))))*&
                     !block(g)%cosBeta(block(g)%nelw1(block(g)%index_ts(index_kp1)))  +  &
                !(block(g)%z1(k+1) - block(g)%zcent(block(g)%nelw1(block(g)%index_ts(index_kp1))))*&
                     !block(g)%cosGamma(block(g)%nelw1(block(g)%index_ts(index_kp1)))
                     temp_pdotn  = (block(g)%xp(i)   - &
                          block(g)%xcent(block(g)%nelp(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosAlpha(block(g)%nelp(block(g)%index_ts(index_kp1)))  +  &
                  (block(g)%yp(j)   - block(g)%ycent(block(g)%nelp(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosBeta(block(g)%nelp(block(g)%index_ts(index_kp1)))   +  &
                  (block(g)%zp(k)   - block(g)%zcent(block(g)%nelp(block(g)%index_ts(index_kp1))))*&
                          block(g)%cosGamma(block(g)%nelp(block(g)%index_ts(index_kp1)))
                   IF (temp_u2dotn<0) THEN
                      block(g)%u(i,j,k+1) = block(g)%u2_ghost(index_kp1)
                   !ELSE
                      !block(g)%u(i,j,k+1) = block(g)%u2t_ghost(index_kp1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j,k+1) = block(g)%u1_ghost(index_kp1)
                   !ELSE
                      !block(g)%u(i-1,j,k+1) = block(g)%u1t_ghost(index_kp1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i,j,k+1) = block(g)%v2_ghost(index_kp1)
                   !ELSE
                      !block(g)%v(i,j,k+1) = block(g)%v2t_ghost(index_kp1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i,j-1,k+1) = block(g)%v1_ghost(index_kp1)
                   !ELSE
                      !block(g)%v(i,j-1,k+1) = block(g)%v1t_ghost(index_kp1)
                   ENDIF
                   IF (temp_w2dotn<0) THEN
                      block(g)%w(i,j,k+1) = block(g)%w2_ghost(index_kp1)
                   !ELSE
                      !block(g)%w(i,j,k+1) = block(g)%w2t_ghost(index_kp1)
                   ENDIF
                   !IF (temp_w1dotn.LT.0) THEN
                      !block(g)%w(i,j,k) = block(g)%w1_ghost(index_kp1)
                   !ELSE
                      !block(g)%w(i,j,k) = block(g)%w1t_ghost(index_kp1)
                   !ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j,k+1) = block(g)%p_ghost(index_kp1)
                   !ELSE
                      !block(g)%p(i,j,k+1) = block(g)%pt_ghost(index_kp1)
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
                   !temp_w2dotn = (block(g)%xp(i)   - &
                   !block(g)%xcent(block(g)%nelw2(block(g)%index_ts(index_km1))))*&
                   !block(g)%cosAlpha(block(g)%nelw2(block(g)%index_ts(index_km1))) +  &
                !(block(g)%yp(j)   - block(g)%ycent(block(g)%nelw2(block(g)%index_ts(index_km1))))*&
                   !block(g)%cosBeta(block(g)%nelw2(block(g)%index_ts(index_km1)))  +  &
               ! (block(g)%z1(k)   - block(g)%zcent(block(g)%nelw2(block(g)%index_ts(index_km1))))*&
                   !block(g)%cosGamma(block(g)%nelw2(block(g)%index_ts(index_km1)))
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
                   !ELSE
                      !block(g)%u(i,j,k-1) = block(g)%u2t_ghost(index_km1)
                   ENDIF
                   IF (temp_u1dotn<0) THEN
                      block(g)%u(i-1,j,k-1) = block(g)%u1_ghost(index_km1)
                   !ELSE
                      !block(g)%u(i-1,j,k-1) = block(g)%u1t_ghost(index_km1)
                   ENDIF
                   IF (temp_v2dotn<0) THEN
                      block(g)%v(i,j,k-1) = block(g)%v2_ghost(index_km1)
                   !ELSE
                      !block(g)%v(i,j,k-1) = block(g)%v2t_ghost(index_km1)
                   ENDIF
                   IF (temp_v1dotn<0) THEN
                      block(g)%v(i,j-1,k-1) = block(g)%v1_ghost(index_km1)
                   !ELSE
                      !block(g)%v(i,j-1,k-1) = block(g)%v1t_ghost(index_km1)
                   ENDIF
                   !IF (temp_w2dotn.LT.0) THEN
               !       block(g)%w(i,j,k-1) = block(g)%w2_ghost(index_km1)
               !    ELSE
               !       block(g)%w(i,j,k-1) = block(g)%w2t_ghost(index_km1)
                !   ENDIF
                   IF (temp_w1dotn<0) THEN
                      block(g)%w(i,j,k-2) = block(g)%w1_ghost(index_km1)
                   !ELSE
                      !block(g)%w(i,j,k-2) = block(g)%w1t_ghost(index_km1)
                   ENDIF
                   IF (temp_pdotn<0) THEN
                      block(g)%p(i,j,k-1) = block(g)%p_ghost(index_km1)
                   !ELSE
                      !block(g)%p(i,j,k-1) = block(g)%pt_ghost(index_km1)
                   ENDIF
                ENDIF

             ENDDO
      ENDIF
! 11  CONTINUE
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

       duutdx=block(g)%u(i,j,k)*(block(g)%ca1_uu(i)*block(g)%u(i+2,j,k)+&
            block(g)%ca2_uu(i)*block(g)%u(i+1,j,k) &
   +block(g)%ca3_uu(i)*block(g)%u(i,j,k)+block(g)%ca4_uu(i)*block(g)%u(i-1,j,k)+block(g)%ca5_uu(i)*&
      block(g)%u(i-2,j,k))/(block(g)%ca6_uu(i)*block(g)%deltax(i+1))+dabs(block(g)%u(i,j,k))* &
      (block(g)%ck1_uu(i)*block(g)%u(i+2,j,k)+block(g)%ck2_uu(i)*block(g)%u(i+1,j,k) &
  +block(g)%ck3_uu(i)*block(g)%u(i,j,k)+block(g)%ck4_uu(i)*block(g)%u(i-1,j,k)+block(g)%ck5_uu(i)*&
      block(g)%u(i-2,j,k))/(2.0_dp*block(g)%ck6_uu(i)*block(g)%deltax(i))

       dvutdy=v_in_um*(block(g)%ca1_vu(j)*block(g)%u(i,j+2,k)+&
            block(g)%ca2_vu(j)*block(g)%u(i,j+1,k) &
            +block(g)%ca3_vu(j)*block(g)%u(i,j,k)+block(g)%ca4_vu(j)*&
            block(g)%u(i,j-1,k)+block(g)%ca5_vu(j)* &
      block(g)%u(i,j-2,k))/(block(g)%ca6_vu(j)*ddye)+dabs(v_in_um)* &
      (block(g)%ck1_vu(j)*block(g)%u(i,j+2,k)+block(g)%ck2_vu(j)*block(g)%u(i,j+1,k) &
      +block(g)%ck3_vu(j)*block(g)%u(i,j,k)+block(g)%ck4_vu(j)*&
      block(g)%u(i,j-1,k)+block(g)%ck5_vu(j)* &
      block(g)%u(i,j-2,k))/(2.0_dp*block(g)%ck6_vu(j)*ddy)

       dwutdz=w_in_um*(block(g)%ca1_wu(k)*block(g)%u(i,j,k+2)+&
            block(g)%ca2_wu(k)*block(g)%u(i,j,k+1) &
            +block(g)%ca3_wu(k)*block(g)%u(i,j,k)+block(g)%ca4_wu(k)*&
            block(g)%u(i,j,k-1)+block(g)%ca5_wu(k)* &
      block(g)%u(i,j,k-2))/(block(g)%ca6_wu(k)*ddzr)+dabs(w_in_um)* &
      (block(g)%ck1_wu(k)*block(g)%u(i,j,k+2)+block(g)%ck2_wu(k)*block(g)%u(i,j,k+1) &
      +block(g)%ck3_wu(k)*block(g)%u(i,j,k)+block(g)%ck4_wu(k)*block(g)%u(i,j,k-1)+&
      block(g)%ck5_wu(k)* &
      block(g)%u(i,j,k-2))/(2.0_dp*block(g)%ck6_wu(k)*ddz)

     	duudx=duutdx
       dvudy=dvutdy
       dwudz=dwutdz

	else
 !cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
   duudx = -0.25_dp*( u1a*u1a + alpha*dabs(u1a)*u22-u3*u3-alpha*dabs(u3)*u4)/block(g)%deltax(i)

   dvudy = -0.25_dp*( v1a*u5 + alpha*dabs(v1a)*u6 - v22*u7 - alpha*dabs(v22)*u8)/block(g)%deltay(j)

   dwudz = -0.25_dp*(w1a*u9 + alpha*dabs(w1a)*u10-w22*u11 - alpha*dabs(w22)*u12) /block(g)%deltaz(k)

       !duudx=-duuwdx
       !dvudy=-dvuwdy
       !dwudz=-dwuwdz

       endif

!ccccccccccccccccccccccccc  grad of u part cccccccccccccccccccccccccccccc
	!d2udx2=(2.0/dx2xr)*((-u4/dxr)+(u22/dx))
       !d2udy2=(2.0/dy)*((-u8/dy2ye)+(u6/dy2yw))
	!d2udz2=(2.0/dz)*((-u12/dz2zt)+(u10/dz2zb))

	d2udx2=(2.0_dp/dx2xr)*((-u4/dxr)+(u22/dx))

       d2udy2=(2.0_dp/(0.5_dp*(dy2ye+dy2yw)))*((-u8/(0.5_dp*dy2ye))+ &
      (u6/(0.5_dp*dy2yw)))

	d2udz2=(2.0_dp/(0.5_dp*(dz2zt+dz2zb)))*((-u12/(0.5_dp*dz2zt))+ &
      (u10/(0.5_dp*dz2zb)))

       xtt2=(d2udx2+d2udy2+d2udz2)/re

       residu=(-duudx-dvudy-dwudz+xtt2)

       !if(ita.eq.1)then
        block(g)%ut(i,j,k)=block(g)%u(i,j,k)+deltat*(residu+dpdx)
      ! else
       ! block(g)%ut(i,j,k)=block(g)%u(i,j,k)+deltat*(0.5*(3.0*residu-block(g)%resi_u(i,j,k))+dpdx)
      ! endif

       block(g)%resi_u(i,j,k)=residu
!c***********************************************************************

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

     duvtdx=u_in_vm*(block(g)%ca1_uw(i)*block(g)%v(i+2,j,k)+block(g)%ca2_uw(i)*block(g)%v(i+1,j,k) &
  +block(g)%ca3_uw(i)*block(g)%v(i,j,k)+block(g)%ca4_uw(i)*block(g)%v(i-1,j,k)+block(g)%ca5_uw(i)* &
      block(g)%v(i-2,j,k))/(block(g)%ca6_uw(i)*ddxr)+dabs(u_in_vm)* &
      (block(g)%ck1_uw(i)*block(g)%v(i+2,j,k)+block(g)%ck2_uw(i)*block(g)%v(i+1,j,k) &
  +block(g)%ck3_uw(i)*block(g)%v(i,j,k)+block(g)%ck4_uw(i)*block(g)%v(i-1,j,k)+block(g)%ck5_uw(i)* &
      block(g)%v(i-2,j,k))/(2.0_dp*block(g)%ck6_uw(i)*ddx)

     dvvtdy=block(g)%v(i,j,k)*(block(g)%ca1_vv(j)*block(g)%v(i,j+2,k)+&
          block(g)%ca2_vv(j)*block(g)%v(i,j+1,k) &
          +block(g)%ca3_vv(j)*block(g)%v(i,j,k)+block(g)%ca4_vv(j)*block(g)%v(i,j-1,k)+&
          block(g)%ca5_vv(j)* &
      block(g)%v(i,j-2,k))/(block(g)%ca6_vv(j)*block(g)%deltay(j+1))+dabs(block(g)%v(i,j,k))* &
      (block(g)%ck1_vv(j)*block(g)%v(i,j+2,k)+block(g)%ck2_vv(j)*block(g)%v(i,j+1,k) &
      +block(g)%ck3_vv(j)*block(g)%v(i,j,k)+block(g)%ck4_vv(j)*block(g)%v(i,j-1,k)+&
      block(g)%ck5_vv(j)* &
      block(g)%v(i,j-2,k))/(2.0_dp*block(g)%ck6_vv(j)*block(g)%deltay(j))

     dwvtdz=w_in_vm*(block(g)%ca1_wu(k)*block(g)%v(i,j,k+2)+block(g)%ca2_wu(k)*block(g)%v(i,j,k+1) &
  +block(g)%ca3_wu(k)*block(g)%v(i,j,k)+block(g)%ca4_wu(k)*block(g)%v(i,j,k-1)+block(g)%ca5_wu(k)* &
      block(g)%v(i,j,k-2))/(block(g)%ca6_wu(k)*ddzr)+dabs(w_in_vm)* &
      (block(g)%ck1_wu(k)*block(g)%v(i,j,k+2)+block(g)%ck2_wu(k)*block(g)%v(i,j,k+1) &
      +block(g)%ck3_wu(k)*block(g)%v(i,j,k)+block(g)%ck4_wu(k)*block(g)%v(i,j,k-1)+&
      block(g)%ck5_wu(k)* &
      block(g)%v(i,j,k-2))/(2.0_dp*block(g)%ck6_wu(k)*ddz)

       duvdx=duvtdx
       dvvdy=dvvtdy
       dwvdz=dwvtdz

       else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
	  duvdx = -0.25_dp*(u13*v3 + alpha*dabs(u13)*v4- u14*v5 - alpha*dabs(u14)*v6)/block(g)%deltax(i)

      dvvdy = -0.25_dp*(v7*v7 + alpha*dabs(v7)*v8 -v9*v9 - alpha*dabs(v9)*v10)/block(g)%deltay(j)

      dwvdz = -0.25_dp*(w3*v11 + alpha*dabs(w3)*v12 -w4*v13 - alpha*dabs(w4)*v14)/block(g)%deltaz(k)

       !duvdx=duvwdx
       !dvvdy=dvvwdy
       !dwvdz=dwvwdz

       endif

!ccccccccccccccccccccccccc  grad of v part cccccccccccccccccccccccccccccc
	!d2vdx2=(2.0/dx)*((-v6/dx2xr)+(v4/dx2xl))
       !d2vdy2=(2.0/dy2ye)*((-v10/dye)+(v8/dy))
	!d2vdz2=(2.0/dz)*((-v14/dz2zt)+(v12/dz2zb))

	d2vdx2=(2.0_dp/(0.5_dp*(dx2xr+dx2xl)))*((-v6/(0.5_dp*dx2xr))+ &
	(v4/(0.5_dp*dx2xl)))

       d2vdy2=(2.0_dp/dy2ye)*((-v10/dye)+(v8/dy))

	d2vdz2=(2.0_dp/(0.5_dp*(dz2zt+dz2zb)))*((-v14/(0.5_dp*dz2zt))+ &
	(v12/(0.5_dp*dz2zb)))

       ytt2=(d2vdx2+d2vdy2+d2vdz2)/re

       residv=(-duvdx-dvvdy-dwvdz+ytt2)

      ! if(ita.eq.1 .or. irest.eq.1)then
       block(g)%vt(i,j,k)=block(g)%v(i,j,k)+deltat*(residv+dpdy)
      ! else
       !block(g)%vt(i,j,k)=block(g)%v(i,j,k)+deltat*(0.5*(3.0*residv-block(g)%resi_v(i,j,k))+dpdy)
      ! endif

       block(g)%resi_v(i,j,k)=residv
!c***********************************************************************

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

       duwtdx=u_in_wm*(block(g)%ca1_uw(i)*block(g)%w(i+2,j,k)+block(g)%ca2_uw(i)*&
            block(g)%w(i+1,j,k) &
            +block(g)%ca3_uw(i)*block(g)%w(i,j,k)+block(g)%ca4_uw(i)*block(g)%w(i-1,j,k)+&
            block(g)%ca5_uw(i)* &
      block(g)%w(i-2,j,k))/(block(g)%ca6_uw(i)*ddxr)+dabs(u_in_wm)* &
      (block(g)%ck1_uw(i)*block(g)%w(i+2,j,k)+block(g)%ck2_uw(i)*block(g)%w(i+1,j,k) &
      +block(g)%ck3_uw(i)*block(g)%w(i,j,k)+block(g)%ck4_uw(i)*block(g)%w(i-1,j,k)+&
      block(g)%ck5_uw(i)* &
      block(g)%w(i-2,j,k))/(2.0_dp*block(g)%ck6_uw(i)*ddx)

       dvwtdy=v_in_wm*(block(g)%ca1_vw(j)*block(g)%w(i,j+2,k)+block(g)%ca2_vw(j)*&
            block(g)%w(i,j+1,k) &
            +block(g)%ca3_vw(j)*block(g)%w(i,j,k)+block(g)%ca4_vw(j)*block(g)%w(i,j-1,k)+&
            block(g)%ca5_vw(j)* &
      block(g)%w(i,j-2,k))/(block(g)%ca6_vw(j)*ddye)+dabs(v_in_wm)* &
      (block(g)%ck1_vw(j)*block(g)%w(i,j+2,k)+block(g)%ck2_vw(j)*block(g)%w(i,j+1,k) &
      +block(g)%ck3_vw(j)*block(g)%w(i,j,k)+block(g)%ck4_vw(j)*block(g)%w(i,j-1,k)+&
       block(g)%ck5_vw(j)* &
      block(g)%w(i,j-2,k))/(2.0_dp*block(g)%ck6_vw(j)*ddy)

       dwwtdz=block(g)%w(i,j,k)*(block(g)%ca1_ww(k)*block(g)%w(i,j,k+2)+block(g)%ca2_ww(k)*&
            block(g)%w(i,j,k+1) &
            +block(g)%ca3_ww(k)*block(g)%w(i,j,k)+block(g)%ca4_ww(k)*block(g)%w(i,j,k-1)+&
            block(g)%ca5_ww(k)* &
      block(g)%w(i,j,k-2))/(block(g)%ca6_ww(k)*block(g)%deltaz(k+1))+dabs(block(g)%w(i,j,k))* &
      (block(g)%ck1_ww(k)*block(g)%w(i,j,k+2)+block(g)%ck2_ww(k)*block(g)%w(i,j,k+1) &
      +block(g)%ck3_ww(k)*block(g)%w(i,j,k)+block(g)%ck4_ww(k)*block(g)%w(i,j,k-1)+&
      block(g)%ck5_ww(k)* &
      block(g)%w(i,j,k-2))/(2.0_dp*block(g)%ck6_ww(k)*block(g)%deltaz(k))

     	duwdx=duwtdx
       dvwdy=dvwtdy
       dwwdz=dwwtdz

	else
!cccccccccccccc---First Order Upwinding ----cccccccccccccccccccccccccccc

   duwdx=-0.5_dp*(u15*w5 + alpha*dabs(u15)*w6 -  u16*w7 - alpha*dabs(u16)*w8)/block(g)%deltax(i)

   dvwdy=-0.5_dp*(v15*w9 + alpha*dabs(v15)*w10 - v16*w11 -alpha*dabs(v16)*w12)/block(g)%deltay(j)

   dwwdz=-0.5_dp*(w13*w13 + alpha*dabs(w13)*w14 - w15*w15 - alpha*dabs(w15)*w16)/block(g)%deltaz(k)

       !duwdx=duwwdx
       !dvwdy=dvwwdy
       !dwwdz=dwwwdz

       endif

!ccccccccccccccccccccccccc  grad of w part cccccccccccccccccccccccccccccc
	!d2wdx2=(2.0/dx)*((-w8/dx2xr)+(w6/dx2xl))
	!d2wdy2=(2.0/dy)*((-w12/dy2ye)+(w10/dy2yw))
	!d2wdz2=(2.0/dz2zt)*((-w16/dzt)+(w14/dz))

	d2wdx2=(2.0_dp/(0.5_dp*(dx2xr+dx2xl)))*((-W8/(0.5_dp*dx2xr))+ &
      (W6/(0.5_dp*dx2xl)))

	d2wdy2=(2.0_dp/(0.5_dp*(dy2ye+dy2yw)))*((-W12/(0.5_dp*dy2ye))+ &
      (W10/(0.5_dp*dy2yw)))

	d2wdz2=(2.0_dp/dz2zt)*((-w16/dzt)+(w14/dz))

       ztt2=(d2wdx2+d2wdy2+d2wdz2)/re

       residw =(-duwdx-dvwdy-dwwdz+ztt2)

       !if(ita.eq.1 .or. irest.eq.1)then
       block(g)%wt(i,j,k)=block(g)%w(i,j,k)+deltat*(residw+dpdz)
       !else
       !block(g)%wt(i,j,k)=block(g)%w(i,j,k)+deltat*(0.5*(3.0*residw-block(g)%resi_w(i,j,k))+dpdz)
       !endif

       block(g)%resi_w(i,j,k)=residw
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
!c
	ENDDO
       !write(6,*) 'leaving nseqcp '
!c***********************************************************************
     !$acc end parallel
        !$acc wait
	ENDDO
      !  !$omp end parallel do
    ! return
      end subroutine nsMomentum2order
!csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!     SUBROUTINE nsMomentumAB
!***********************************************************************
!     navier-stokes equations for constant properties
!***********************************************************************
!        USE global
!***********************************************************************
!     !write(*,*)'has entered nseqcp'
!         IMPLICIT NONE
!         INTEGER          :: i, j, k, n
!     REAL (KIND = 8)  :: dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,         &
!                         u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,       &
!                w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
!       	   dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
!       	  dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          &
!       	   r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
!       	  dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
!       		   uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,        &
!                    dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w, &
!       	  u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
!      	        duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw
!
!       al = 1.

!!$acc parallel loop gang vector  &
!!$acc private (i, j, k, dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,   &
!!$acc	        u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
!!$acc		 w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
!!$acc	        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
!!$acc		 dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          &
!!$acc		 r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
!!$acc	        dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
!!$acc	        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,            &
!!$acc	        dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
!!$acc	        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
!!$acc		 duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw) &
!!$acc present (block(g)%fluidIndexPtr, block(g)%deltax, block(g)%deltay, block(g)%deltaz, u, v, w,&
!!$acc          ut, vt, wt, p, cell, x1, y1, z1, xp, yp, zp,                &
!!$acc          xu, yu, zu, xv, yv, zv, xw, yw, zw, resi_u, resi_v, resi_w, &
!!$acc          ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
!!$acc          ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
!!$acc          ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
!!$acc          ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
!!$acc          ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
!!$acc          ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
!!$acc          ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
!!$acc          ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
!!$acc          ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
!!$acc          ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
!!$acc          ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
!!$acc          ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
!!$acc          ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
!!$acc          ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
!!$acc          ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
!!$acc          ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
!!$acc          ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
!!$acc          ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww) &
!!$acc firstprivate(nx_var, ny_var, nz_var, rev, deltat, al)
!       DO n = 1, block(g)%fluidCellCount
!
!      i = block(g)%fluidIndexPtr(n, 1)
!      j = block(g)%fluidIndexPtr(n, 2)
!      k = block(g)%fluidIndexPtr(n, 3)
!
!       dxr=block(g)%deltax(i+1)
!       dx=block(g)%deltax(i)
!       dxl=block(g)%deltax(i-1)
!       dye=block(g)%deltay(j+1)
!       dy=block(g)%deltay(j)
!       dyw=block(g)%deltay(j-1)
!       dzt=block(g)%deltaz(k+1)
!       dz=block(g)%deltaz(k)
!       dzb=block(g)%deltaz(k-1)
!cccccccccccccccccccccccccc  diff-u     ccccccccccccccccccccccccccccccccc
!     duu / dx
!      u1a = block(g)%u(i-1,j,k) + block(g)%u(i,j,k)
!      u22 = block(g)%u(i-1,j,k) - block(g)%u(i,j,k)
!      u3  = block(g)%u(i,j,k)   + block(g)%u(i+1,j,k)
!      u4  = block(g)%u(i,j,k)   - block(g)%u(i+1,j,k)

!     duv  / dy
!      u5 = block(g)%u(i,j-1,k)  + block(g)%u(i,j,k)
!      u6 = block(g)%u(i,j-1,k)  - block(g)%u(i,j,k)
!      u7 = block(g)%u(i,j,k)    + block(g)%u(i,j+1,k)
!      u8 = block(g)%u(i,j,k)    - block(g)%u(i,j+1,k)

!     dwu  /  dz
!      u9 = block(g)%u(i,j,k-1)  + block(g)%u(i,j,k)
!      u10= block(g)%u(i,j,k-1)  - block(g)%u(i,j,k)
!      u11= block(g)%u(i,j,k)    + block(g)%u(i,j,k+1)
!      u12= block(g)%u(i,j,k)    - block(g)%u(i,j,k+1)

!     duv/dx
!      u13 = block(g)%u(i-1,j,k) + block(g)%u(i-1,j+1,k)
!      u14 = u7

!     duw / dx
!      u15 = block(g)%u(i-1,j,k) + block(g)%u(i-1,j,k+1)
!      u16 = u11
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccc    diff -v    ccccccccccccccccccccccccccccccc
!     dvu / dx
!      v1a = block(g)%v(i,j-1,k)  + block(g)%v(i+1,j-1,k)
!      v22 = block(g)%v(i,j,k)    + block(g)%v(i+1,j,k)

!     duv / dx
!      v3 = block(g)%v(i-1,j,k)   + block(g)%v(i,j,k)
!      v4 = block(g)%v(i-1,j,k)   - block(g)%v(i,j,k)
!      v5 = block(g)%v(i,j,k)     + block(g)%v(i+1,j,k)
!      v6 = block(g)%v(i,j,k)     - block(g)%v(i+1,j,k)

!     dvv / dy
!      v7 = block(g)%v(i,j-1,k)   + block(g)%v(i,j,k)
!      v8 = block(g)%v(i,j-1,k)   - block(g)%v(i,j,k)
!      v9 = block(g)%v(i,j,k)     + block(g)%v(i,j+1,k)
!      v10= block(g)%v(i,j,k)     - block(g)%v(i,j+1,k)

!     dwu / dz
!      v11 = block(g)%v(i,j,k-1)  + block(g)%v(i,j,k)
!      v12 = block(g)%v(i,j,k-1)  - block(g)%v(i,j,k)
!      v13 = block(g)%v(i,j,k)    + block(g)%v(i,j,k+1)
!      v14 = block(g)%v(i,j,k)    - block(g)%v(i,j,k+1)

!     dvw / dy
!      v15 = block(g)%v(i,j-1,k)  + block(g)%v(i,j-1,k+1)
!      v16 = v13
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccc  diff - w  cccccccccccccccccccccccccccccccc
!     dwu / dz
!      w1a = block(g)%w(i,j,k-1)  + block(g)%w(i+1,j,k-1)
!      w22 = block(g)%w(i,j,k)    + block(g)%w(i+1,j,k)

!     dwv / dz
!      w3 = block(g)%w(i,j,k-1)   + block(g)%w(i,j+1,k-1)
!      w4 = block(g)%w(i,j,k)     + block(g)%w(i,j+1,k)

!     duw / dx
!      w5 = block(g)%w(i-1,j,k)   + block(g)%w(i,j,k)
!      w6 = block(g)%w(i-1,j,k)   - block(g)%w(i,j,k)
!      w7 = block(g)%w(i,j,k)     + block(g)%w(i+1,j,k)
!      w8 = block(g)%w(i,j,k)     - block(g)%w(i+1,j,k)

!     dvw / dy
!      w9 = block(g)%w(i,j-1,k)   + block(g)%w(i,j,k)
!      w10 = block(g)%w(i,j-1,k)  - block(g)%w(i,j,k)
!      w11 = block(g)%w(i,j,k)    + block(g)%w(i,j+1,k)
!      w12 = block(g)%w(i,j,k)    - block(g)%w(i,j+1,k)

!     dww / dz
!      w13 = block(g)%w(i,j,k-1)  + block(g)%w(i,j,k)
!      w14 = block(g)%w(i,j,k-1)  - block(g)%w(i,j,k)
!      w15 = block(g)%w(i,j,k)    + block(g)%w(i,j,k+1)
!      w16 = block(g)%w(i,j,k)    - block(g)%w(i,j,k+1)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      dpdx = (block(g)%p(i,j,k) - block(g)%p(i+1,j,k))/(0.5*(dxr+dx))
!      dpdy = (block(g)%p(i,j,k) - block(g)%p(i,j+1,k))/(0.5*(dye+dy))
!      dpdz = (block(g)%p(i,j,k) - block(g)%p(i,j,k+1))/(0.5*(dzt+dz))
!
!   	dx2xr = dx + dxr
!      dx2xl = dx + dxl
!       dy2ye = dy + dye
!       dy2yw = dy + dyw
!       dz2zt = dz + dzt
!       dz2zb = dz + dzb
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*********************** U - Momentum **********************************
!       vu_e=block(g)%v(i+1,j,k)+(dxr/(dx2xr))*(block(g)%v(i,j,k)-block(g)%v(i+1,j,k))
!       vu_w=block(g)%v(i+1,j-1,k)+(dxr/(dx2xr))*(block(g)%v(i,j-1,k)-block(g)%v(i+1,j-1,k))
!       v_in_um=0.5*(vu_e+vu_w)
!
!       wu_n=block(g)%w(i+1,j,k)+(dxr/(dx2xr))*(block(g)%w(i,j,k)-block(g)%w(i+1,j,k))
!       wu_s=block(g)%w(i+1,j,k-1)+(dxr/(dx2xr))*(block(g)%w(i,j,k-1)-block(g)%w(i+1,j,k-1))
!       w_in_um=0.5*(wu_n+wu_s)
!
!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
!      if(i.ne.2.and.i.lt.nx_var.and.j.ne.2.and.j.lt.ny_var+1.and.    &
!      k.ne.2.and.k.lt.nz_var+1.and.block(g)%cell(i+1,j,k).ne.2.and.       &
!      block(g)%cell(i-1,j,k).ne.2.and.block(g)%cell(i,j+1,k).ne.2.and.         &
!      block(g)%cell(i,j-1,k).ne.2.and.block(g)%cell(i,j,k+1).ne.2.and.         &
!      block(g)%cell(i,j,k-1).ne.2) then
!
!       ddy=0.5*(block(g)%deltay(j)+block(g)%deltay(j-1))
!      ddye=0.5*(block(g)%deltay(j)+block(g)%deltay(j+1))
!      ddz=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
!      ddzr=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k+1))

      !      duutdx=block(g)%u(i,j,k)*(block(g)%ca1_uu(i)*block(g)%u(i+2,j,k)+&
      !block(g)%ca2_uu(i)*block(g)%u(i+1,j,k) &
      !      +block(g)%ca3_uu(i)*block(g)%u(i,j,k)+block(g)%ca4_uu(i)*block(g)%u(i-1,j,k)+&
      !block(g)%ca5_uu(i)* &
!      block(g)%u(i-2,j,k))/(block(g)%ca6_uu(i)*block(g)%deltax(i+1))+dabs(block(g)%u(i,j,k))* &
!      (block(g)%ck1_uu(i)*block(g)%u(i+2,j,k)+block(g)%ck2_uu(i)*block(g)%u(i+1,j,k) &
      !      +block(g)%ck3_uu(i)*block(g)%u(i,j,k)+block(g)%ck4_uu(i)*block(g)%u(i-1,j,k)+&
      !block(g)%ck5_uu(i)* &
!      block(g)%u(i-2,j,k))/(2.0*block(g)%ck6_uu(i)*block(g)%deltax(i))

      !      dvutdy=v_in_um*(block(g)%ca1_vu(j)*block(g)%u(i,j+2,k)+block(g)%ca2_vu(j)*&
      !block(g)%u(i,j+1,k) &
      !      +block(g)%ca3_vu(j)*block(g)%u(i,j,k)+block(g)%ca4_vu(j)*block(g)%u(i,j-1,k)+&
      !block(g)%ca5_vu(j)* &
!      block(g)%u(i,j-2,k))/(block(g)%ca6_vu(j)*ddye)+dabs(v_in_um)* &
!      (block(g)%ck1_vu(j)*block(g)%u(i,j+2,k)+block(g)%ck2_vu(j)*block(g)%u(i,j+1,k) &
      !      +block(g)%ck3_vu(j)*block(g)%u(i,j,k)+block(g)%ck4_vu(j)*block(g)%u(i,j-1,k)+&
      !block(g)%ck5_vu(j)* &
!      block(g)%u(i,j-2,k))/(2.0*block(g)%ck6_vu(j)*ddy)

      !      dwutdz=w_in_um*(block(g)%ca1_wu(k)*block(g)%u(i,j,k+2)+block(g)%ca2_wu(k)*
      !block(g)%u(i,j,k+1) &
      !      +block(g)%ca3_wu(k)*block(g)%u(i,j,k)+block(g)%ca4_wu(k)*block(g)%u(i,j,k-1)+&
      !block(g)%ca5_wu(k)* &
!      block(g)%u(i,j,k-2))/(block(g)%ca6_wu(k)*ddzr)+dabs(w_in_um)* &
!      (block(g)%ck1_wu(k)*block(g)%u(i,j,k+2)+block(g)%ck2_wu(k)*block(g)%u(i,j,k+1) &
      !      +block(g)%ck3_wu(k)*block(g)%u(i,j,k)+block(g)%ck4_wu(k)*block(g)%u(i,j,k-1)+&
      !block(g)%ck5_wu(k)* &
!      block(g)%u(i,j,k-2))/(2.0*block(g)%ck6_wu(k)*ddz)
!
!    	duudx=duutdx
!      dvudy=dvutdy
!      dwudz=dwutdz
!
!       else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
!       r1x=block(g)%deltax(i)/block(g)%deltax(i+1)
!      r1yn=0.5*(block(g)%deltay(j)+block(g)%deltay(j-1))
!      r1yd=0.5*(block(g)%deltay(j)+block(g)%deltay(j+1))
!      r1y=r1yn/r1yd
!       r1zn=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
!       r1zd=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k+1))
!       r1z=r1zn/r1zd

!       duuwdx=(block(g)%u(i,j,k)/block(g)%deltax(i+1))*((-1.0/(r1x*(r1x+1.0)))* &
!    	block(g)%u(i-1,j,k)-(1.0-1.0/r1x)*block(g)%u(i,j,k)+   &
!    	(r1x/(r1x+1.0))*block(g)%u(i+1,j,k))-al*dabs(block(g)%u(i,j,k))* &
!    	(1.0/(2.0*block(g)%deltax(i+1)))*(2.0/(r1x*(r1x+1.0))* &
!    	block(g)%u(i-1,j,k)-2.0/r1x*block(g)%u(i,j,k)+2.0/(r1x+1.0)* &
!      block(g)%u(i+1,j,k))

!      dvuwdy=(v_in_um/r1yd)*((-1.0/(r1y*(r1y+1.0)))* &
!      block(g)%u(i,j-1,k)-(1.0-1.0/r1y)*block(g)%u(i,j,k)+ &
!      (r1y/(r1y+1.0))*block(g)%u(i,j+1,k))-al*dabs(v_in_um)* &
!      (1.0/(2.0*r1yd))*(2.0/(r1y*(r1y+1.0))* &
!      block(g)%u(i,j-1,k)-2.0/r1y*block(g)%u(i,j,k)+2.0/(r1y+1.0)* &
!      block(g)%u(i,j+1,k))

!       dwuwdz=(w_in_um/r1zd)*((-1.0/(r1z*(r1z+1.0)))* &
!    	block(g)%u(i,j,k-1)-(1.0-1.0/r1z)*block(g)%u(i,j,k)+ &
!   	(r1z/(r1z+1.0))*block(g)%u(i,j,k+1))-al*dabs(w_in_um)* &
!    	(1.0/(2.0*r1zd))*(2.0/(r1z*(r1z+1.0))* &
!    	block(g)%u(i,j,k-1)-2.0/r1z*block(g)%u(i,j,k)+2.0/(r1z+1.0)* &
!     	block(g)%u(i,j,k+1))
!
!      duudx=duuwdx
!      dvudy=dvuwdy
!      dwudz=dwuwdz
!
!      endif
!
!ccccccccccccccccccccccccc  grad of u part cccccccccccccccccccccccccccccc
!       !d2udx2=(2.0/dx2xr)*((-u4/dxr)+(u22/dx))
!      !d2udy2=(2.0/dy)*((-u8/dy2ye)+(u6/dy2yw))
!       !d2udz2=(2.0/dz)*((-u12/dz2zt)+(u10/dz2zb))
!
!       d2udx2=(2.0/dx2xr)*((-u4/dxr)+(u22/dx))
!
!      d2udy2=(2.0/(0.5*(dy2ye+dy2yw)))*((-u8/(0.5*dy2ye))+ &
!      (u6/(0.5*dy2yw)))
!
!       d2udz2=(2.0/(0.5*(dz2zt+dz2zb)))*((-u12/(0.5*dz2zt))+ &
!      (u10/(0.5*dz2zb)))

!      xtt2=rev*(d2udx2+d2udy2+d2udz2)
!
!      residu=(-duudx-dvudy-dwudz+xtt2)

!      !if(ita.eq.1)then
!      !block(g)%ut(i,j,k)=block(g)%u(i,j,k)+deltat*(residu+dpdx)
!      !else
!      block(g)%ut(i,j,k)=block(g)%u(i,j,k)+deltat*(0.5*(3.0*residu-block(g)%resi_u(i,j,k))+dpdx)
!      !endif

!      block(g)%resi_u(i,j,k)=residu
!c***********************************************************************

!c*********************** V - Momentum *********************************
!       uv_e=block(g)%u(i,j+1,k)+(dye/(dy2ye))*(block(g)%u(i,j,k)-block(g)%u(i,j+1,k))
!       uv_w=block(g)%u(i-1,j+1,k)+(dye/(dy2ye))*(block(g)%u(i-1,j,k)-block(g)%u(i-1,j+1,k))
!       u_in_vm=0.5*(uv_e+uv_w)
!
!       wv_n=block(g)%w(i,j+1,k)+(dye/(dy2ye))*(block(g)%w(i,j,k)-block(g)%w(i,j+1,k))
!       wv_s=block(g)%w(i,j+1,k-1)+(dye/(dy2ye))*(block(g)%w(i,j,k-1)-block(g)%w(i,j+1,k-1))
!       w_in_vm=0.5*(wv_n+wv_s)

!cccccccccccccc---Third Order Upwinding ----cccccccccccccccccccccccccccc
!     if(i.ne.2.and.i.lt.nx_var+1.and.j.ne.2.and.j.lt.ny_var.and. &
!     k.ne.2.and.k.lt.nz_var+1.and.block(g)%cell(i+1,j,k).ne.2.and.    &
!     block(g)%cell(i-1,j,k).ne.2.and.block(g)%cell(i,j+1,k).ne.2.and.      &
!     block(g)%cell(i,j-1,k).ne.2.and.block(g)%cell(i,j,k+1).ne.2.and.      &
!     block(g)%cell(i,j,k-1).ne.2) then
!
!       ddx=0.5*(block(g)%deltax(i)+block(g)%deltax(i-1))
!      ddxr=0.5*(block(g)%deltax(i)+block(g)%deltax(i+1))
!      ddz=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
!      ddzr=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k+1))

      !      duvtdx=u_in_vm*(block(g)%ca1_uw(i)*block(g)%v(i+2,j,k)+block(g)%ca2_uw(i)*&
      !block(g)%v(i+1,j,k) &
      !      +block(g)%ca3_uw(i)*block(g)%v(i,j,k)+block(g)%ca4_uw(i)*block(g)%v(i-1,j,k)+&
      !block(g)%ca5_uw(i)* &
!      block(g)%v(i-2,j,k))/(block(g)%ca6_uw(i)*ddxr)+dabs(u_in_vm)* &
!      (block(g)%ck1_uw(i)*block(g)%v(i+2,j,k)+block(g)%ck2_uw(i)*block(g)%v(i+1,j,k) &
      !      +block(g)%ck3_uw(i)*block(g)%v(i,j,k)+block(g)%ck4_uw(i)*block(g)%v(i-1,j,k)+&
      !block(g)%ck5_uw(i)* &
!      block(g)%v(i-2,j,k))/(2.0*block(g)%ck6_uw(i)*ddx)

      !      dvvtdy=block(g)%v(i,j,k)*(block(g)%ca1_vv(j)*block(g)%v(i,j+2,k)+&
      !block(g)%ca2_vv(j)*block(g)%v(i,j+1,k) &
      !      +block(g)%ca3_vv(j)*block(g)%v(i,j,k)+block(g)%ca4_vv(j)*block(g)%v(i,j-1,k)+&
      !block(g)%ca5_vv(j)* &
!      block(g)%v(i,j-2,k))/(block(g)%ca6_vv(j)*block(g)%deltay(j+1))+dabs(block(g)%v(i,j,k))* &
!      (block(g)%ck1_vvj)*block(g)%v(i,j+2,k)+block(g)%ck2_vvj)*block(g)%v(i,j+1,k) &
      !      +block(g)%ck3_vvj)*block(g)%v(i,j,k)+block(g)%ck4_vvj)*block(g)%v(i,j-1,k)+&
      !block(g)%ck5_vv(j)* &
!      block(g)%v(i,j-2,k))/(2.0*block(g)%ck6_vv(j)*block(g)%deltay(j))

      !      dwvtdz=w_in_vm*(block(g)%ca1_wu(k)*block(g)%v(i,j,k+2)+block(g)%ca2_wu(k)*&
      !block(g)%v(i,j,k+1) &
      !      +block(g)%ca3_wu(k)*block(g)%v(i,j,k)+block(g)%ca4_wu(k)*block(g)%v(i,j,k-1)+&
      !block(g)%ca5_wu(k)* &
!      block(g)%v(i,j,k-2))/(block(g)%ca6_wu(k)*ddzr)+dabs(w_in_vm)* &
!      (block(g)%ck1_wu(k)*block(g)%v(i,j,k+2)+block(g)%ck2_wu(k)*block(g)%v(i,j,k+1) &
      !      +block(g)%ck3_wu(k)*block(g)%v(i,j,k)+block(g)%ck4_wu(k)*block(g)%v(i,j,k-1)+&
      !block(g)%ck5_wu(k)* &
!      block(g)%v(i,j,k-2))/(2.0*block(g)%ck6_wu(k)*ddz)
!
!      duvdx=duvtdx
!      dvvdy=dvvtdy
!      dwvdz=dwvtdz
!
!      else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
!       r1xn=0.5*(block(g)%deltax(i)+block(g)%deltax(i-1))
!       r1xd=0.5*(block(g)%deltax(i)+block(g)%deltax(i+1))
!       r1x=r1xn/r1xd
!      r1y=block(g)%deltay(j)/block(g)%deltay(j+1)
!       r1zn=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k-1))
!       r1zd=0.5*(block(g)%deltaz(k)+block(g)%deltaz(k+1))
!       r1z=r1zn/r1zd

!       duvwdx=(u_in_vm/r1xd)*((-1.0/(r1x*(r1x+1.0)))* &
!    	block(g)%v(i-1,j,k)-(1.0-1.0/r1x)*block(g)%v(i,j,k)+ &
!    	(r1x/(r1x+1.0))*block(g)%v(i+1,j,k))-al*dabs(u_in_vm)* &
!    	(1.0/(2.0*r1xd))*(2.0/(r1x*(r1x+1.0))* &
!    	block(g)%v(i-1,j,k)-2.0/r1x*block(g)%v(i,j,k)+2.0/(r1x+1.0)* &
!    	block(g)%v(i+1,j,k))

!      dvvwdy=(block(g)%v(i,j,k)/block(g)%deltay(j+1))*((-1.0/(r1y*(r1y+1.0)))* &
!      block(g)%v(i,j-1,k)-(1.0-1.0/r1y)*block(g)%v(i,j,k)+ &
!      (r1y/(r1y+1.0))*block(g)%v(i,j+1,k))-al*dabs(block(g)%v(i,j,k))* &
!      (1.0/(2.0*block(g)%deltay(j+1)))*(2.0/(r1y*(r1y+1.0))* &
!      block(g)%v(i,j-1,k)-2.0/r1y*block(g)%v(i,j,k)+2.0/(r1y+1.0)* &
!      block(g)%v(i,j+1,k))

!       dwvwdz=(w_in_vm/r1zd)*((-1.0/(r1z*(r1z+1.0)))* &
!    	block(g)%v(i,j,k-1)-(1.0-1.0/r1z)*block(g)%v(i,j,k)+ &
!    	(r1z/(r1z+1.0))*block(g)%v(i,j,k+1))-al*dabs(w_in_vm)* &
!    	(1.0/(2.0*r1zd))*(2.0/(r1z*(r1z+1.0))* &
!    	block(g)%v(i,j,k-1)-2.0/r1z*block(g)%v(i,j,k)+2.0/(r1z+1.0)* &
!    	block(g)%v(i,j,k+1))

!      duvdx=duvwdx
!      dvvdy=dvvwdy
!      dwvdz=dwvwdz
!
!      endif

!ccccccccccccccccccccccccc  grad of v part cccccccccccccccccccccccccccccc
!       !d2vdx2=(2.0/dx)*((-v6/dx2xr)+(v4/dx2xl))
!      !d2vdy2=(2.0/dy2ye)*((-v10/dye)+(v8/dy))
!       !d2vdz2=(2.0/dz)*((-v14/dz2zt)+(v12/dz2zb))
!
!       d2vdx2=(2.0/(0.5*(dx2xr+dx2xl)))*((-v6/(0.5*dx2xr))+ &
!       (v4/(0.5*dx2xl)))
!
!      d2vdy2=(2.0/dy2ye)*((-v10/dye)+(v8/dy))
!
!       d2vdz2=(2.0/(0.5*(dz2zt+dz2zb)))*((-v14/(0.5*dz2zt))+ &
!       (v12/(0.5*dz2zb)))
!
!      ytt2=rev*(d2vdx2+d2vdy2+d2vdz2)

!      residv=(-duvdx-dvvdy-dwvdz+ytt2)

!      !if(ita.eq.1 .or. irest.eq.1)then
!      !block(g)%vt(i,j,k)=block(g)%v(i,j,k)+deltat*(residv+dpdy)
!      !else
!      block(g)%vt(i,j,k)=block(g)%v(i,j,k)+deltat*(0.5*(3.0*residv-block(g)%resi_v(i,j,k))+dpdy)
!      !endif

!      block(g)%resi_v(i,j,k)=residv
!c***********************************************************************

!c*********************** W - Momentum **********************************
!       uw_e=block(g)%u(i,j,k+1)+(dzt/(dz2zt))*(block(g)%u(i,j,k)-block(g)%u(i,j,k+1))
!       uw_w=block(g)%u(i-1,j,k+1)+(dzt/(dz2zt))*(block(g)%u(i-1,j,k)-block(g)%u(i-1,j,k+1))
!       u_in_wm=0.5*(uw_e+uw_w)
!
!       vw_n=block(g)%v(i,j,k+1)+(dzt/(dz2zt))*(block(g)%v(i,j,k)-block(g)%v(i,j,k+1))
!       vw_s=block(g)%v(i,j-1,k+1)+(dzt/(dz2zt))*(block(g)%v(i,j-1,k)-block(g)%v(i,j-1,k+1))
!       v_in_wm=0.5*(vw_n+vw_s)

!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
!     if(i.ne.2.and.i.lt.nx_var+1.and.j.ne.2.and.j.lt.ny_var+1.and. &
!     k.ne.2.and.k.lt.nz_var.and.block(g)%cell(i+1,j,k).ne.2.and.        &
!     block(g)%cell(i-1,j,k).ne.2.and.block(g)%cell(i,j+1,k).ne.2.and.        &
!     block(g)%cell(i,j-1,k).ne.2.and.block(g)%cell(i,j,k+1).ne.2.and.        &
!     block(g)%cell(i,j,k-1).ne.2) then

!       ddx=0.5*(block(g)%deltax(i)+block(g)%deltax(i-1))
!      ddxr=0.5*(block(g)%deltax(i)+block(g)%deltax(i+1))
!      ddy=0.5*(block(g)%deltay(j)+block(g)%deltay(j-1))
!      ddye=0.5*(block(g)%deltay(j)+block(g)%deltay(j+1))

      !      duwtdx=u_in_wm*(block(g)%ca1_uw(i)*block(g)%w(i+2,j,k)+block(g)%ca2_uw(i)*&
      !block(g)%w(i+1,j,k) &
      !      +block(g)%ca3_uw(i)*block(g)%w(i,j,k)+block(g)%ca4_uw(i)*block(g)%w(i-1,j,k)+&
      !block(g)%ca5_uw(i)* &
!      block(g)%w(i-2,j,k))/(block(g)%ca6_uw(i)*ddxr)+dabs(u_in_wm)* &
!      (block(g)%ck1_uw(i)*block(g)%w(i+2,j,k)+block(g)%ck2_uw(i)*block(g)%w(i+1,j,k) &
      !      +block(g)%ck3_uw(i)*block(g)%w(i,j,k)+block(g)%ck4_uw(i)*block(g)%w(i-1,j,k)+&
      !block(g)%ck5_uw(i)* &
!      block(g)%w(i-2,j,k))/(2.0*block(g)%ck6_uw(i)*ddx)

      !      dvwtdy=v_in_wm*(block(g)%ca1_vw(j)*block(g)%w(i,j+2,k)+block(g)%ca2_vw(j)*&
      !block(g)%w(i,j+1,k) &
      !      +block(g)%ca3_vw(j)*block(g)%w(i,j,k)+block(g)%ca4_vw(j)*block(g)%w(i,j-1,k)+!
      !block(g)%ca5_vw(j)* &
!      block(g)%w(i,j-2,k))/(block(g)%ca6_vw(j)*ddye)+dabs(v_in_wm)* &
!      (block(g)%ck1_vw(j)*block(g)%w(i,j+2,k)+block(g)%ck2_vw(j)*block(g)%w(i,j+1,k) &
      !      +block(g)%ck3_vw(j)*block(g)%w(i,j,k)+block(g)%ck4_vw(j)*block(g)%w(i,j-1,k)+&
      !block(g)%ck5_vw(j)* &
!      block(g)%w(i,j-2,k))/(2.0*block(g)%ck6_vw(j)*ddy)

      !      dwwtdz=block(g)%w(i,j,k)*(block(g)%ca1_ww(k)*block(g)%w(i,j,k+2)+block(g)%ca2_ww(k)*&
      !block(g)%w(i,j,k+1) &
      !      +block(g)%ca3_ww(k)*block(g)%w(i,j,k)+block(g)%ca4_ww(k)*block(g)%w(i,j,k-1)+&
      !block(g)%ca5_ww(k)* &
!      block(g)%w(i,j,k-2))/(block(g)%ca6_ww(k)*block(g)%deltaz(k+1))+dabs(block(g)%w(i,j,k))* &
!      (block(g)%ck1_ww(k)*block(g)%w(i,j,k+2)+block(g)%ck2_ww(k)*block(g)%w(i,j,k+1) &
      !      +block(g)%ck3_ww(k)*block(g)%w(i,j,k)+block(g)%ck4_ww(k)*block(g)%w(i,j,k-1)+&
      !block(g)%ck5_ww(k)* &
!      block(g)%w(i,j,k-2))/(2.0*block(g)%ck6_ww(k)*block(g)%deltaz(k))
!
!    	duwdx=duwtdx
!      dvwdy=dvwtdy
!      dwwdz=dwwtdz

!       else
!cccccccccccccc---First Order Upwinding ----cccccccccccccccccccccccccccc
!       r1xn=0.5*(block(g)%deltax(i)+block(g)%deltax(i-1))
!       r1xd=0.5*(block(g)%deltax(i)+block(g)%deltax(i+1))
!       r1x=r1xn/r1xd
!      r1yn=0.5*(block(g)%deltay(j)+block(g)%deltay(j-1))
!      r1yd=0.5*(block(g)%deltay(j)+block(g)%deltay(j+1))
!      r1y=r1yn/r1yd
!       r1z=block(g)%deltaz(k)/block(g)%deltaz(k+1)

!       duwwdx=(u_in_wm/r1xd)*((-1.0/(r1x*(r1x+1.0)))* &
!    	block(g)%w(i-1,j,k)-(1.0-1.0/r1x)*block(g)%w(i,j,k)+ &
!    	(r1x/(r1x+1.0))*block(g)%w(i+1,j,k))-al*dabs(u_in_wm)* &
!    	(1.0/(2.0*r1xd))*(2.0/(r1x*(r1x+1.0))* &
!    	block(g)%w(i-1,j,k)-2.0/r1x*block(g)%w(i,j,k)+2.0/(r1x+1.0)* &
!    	block(g)%w(i+1,j,k))

!      dvwwdy=(v_in_wm/r1yd)*((-1.0/(r1y*(r1y+1.0)))* &
!      block(g)%w(i,j-1,k)-(1.0-1.0/r1y)*block(g)%w(i,j,k)+ &
!      (r1y/(r1y+1.0))*block(g)%w(i,j+1,k))-al*dabs(v_in_wm)* &
!      (1.0/(2.0*r1yd))*(2.0/(r1y*(r1y+1.0))* &
!      block(g)%w(i,j-1,k)-2.0/r1y*block(g)%w(i,j,k)+2.0/(r1y+1.0)* &
!      block(g)%w(i,j+1,k))

!       dwwwdz=(block(g)%w(i,j,k)/block(g)%deltaz(k+1))*((-1.0/(r1z*(r1z+1.0)))* &
!    	block(g)%w(i,j,k-1)-(1.0-1.0/r1z)*block(g)%w(i,j,k)+ &
!    	(r1z/(r1z+1.0))*block(g)%w(i,j,k+1))-al*dabs(block(g)%w(i,j,k))* &
!    	(1.0/(2.0*block(g)%deltaz(k+1)))*(2.0/(r1z*(r1z+1.0))* &
!    	block(g)%w(i,j,k-1)-2.0/r1z*block(g)%w(i,j,k)+2.0/(r1z+1.0)* &
!    	block(g)%w(i,j,k+1))
!
!
!      duwdx=duwwdx
!      dvwdy=dvwwdy
!      dwwdz=dwwwdz
!
!      endif

!ccccccccccccccccccccccccc  grad of w part cccccccccccccccccccccccccccccc
!       !d2wdx2=(2.0/dx)*((-w8/dx2xr)+(w6/dx2xl))
!       !d2wdy2=(2.0/dy)*((-w12/dy2ye)+(w10/dy2yw))
!       !d2wdz2=(2.0/dz2zt)*((-w16/dzt)+(w14/dz))
!
!       d2wdx2=(2.0/(0.5*(dx2xr+dx2xl)))*((-W8/(0.5*dx2xr))+ &
!      (W6/(0.5*dx2xl)))
!
!       d2wdy2=(2.0/(0.5*(dy2ye+dy2yw)))*((-W12/(0.5*dy2ye))+ &
!      (W10/(0.5*dy2yw)))
!
!       d2wdz2=(2.0/dz2zt)*((-w16/dzt)+(w14/dz))

!      ztt2=rev*(d2wdx2+d2wdy2+d2wdz2)

!      residw =(-duwdx-dvwdy-dwwdz+ztt2)

!      !if(ita.eq.1 .or. irest.eq.1)then
!      !block(g)%wt(i,j,k)=w(i,j,k)+deltat*(residw+dpdz)
!      !else
!      block(g)%wt(i,j,k)=block(g)%w(i,j,k)+deltat*(0.5*(3.0*residw-block(g)%resi_w(i,j,k))+dpdz)
!      !endif

!      block(g)%resi_w(i,j,k)=residw
!***********************************************************************
!       ENDDO
!!$acc end parallel
!      !write(6,*) 'leaving nseqcp '
!***********************************************************************
!     END SUBROUTINE nsMomentumAB

end module biocfd_navier_stokes

module biocfd_coarse_update
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use global, only : block, intfr, intflines
  use biocfd_interpolation, only: bilinear_interpolation, linear_interpolation
#ifdef BIOCFD_MPI
  use mpi_f08
#endif
  implicit none
  private

  public ::  coarseUpdate_newv, coarseUpdate_pc, coarseUpdate

  contains
subroutine coarseUpdate
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

#ifdef BIOCFD_MPI
        integer :: rank, num_proc, ierror, token
        call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

        ! TN: I *think* that this interpolation needs to essentially
        ! happen serially across ranks. The idea is that rank 0 will
        ! run first, then update block(1) on rank 1 and tell rank 1 it
        ! can start. Rank 1 would then update the block on rank 2 and
        ! tell rank 2 it can start. This process will repeat up to
        ! rank n. Once rank n finishes we broadcast the block(1)
        ! arrays to all ranks so that everything is up to date.
        !
        ! We don't wait for anything if we are on rank 0 and just get
        ! started
        if (rank /= 0) then
          ! This essentially acts as a barrier in that rank n wont
          ! start until rank n-1 has finished. Assumes that a_blk number is always 1
          call MPI_Recv(block(1)%ut, size(block(1)%ut), MPI_DOUBLE_PRECISION, rank-1, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
          call MPI_Recv(block(1)%vt, size(block(1)%vt), MPI_DOUBLE_PRECISION, rank-1, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
          call MPI_Recv(block(1)%wt, size(block(1)%wt), MPI_DOUBLE_PRECISION, rank-1, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        end if

#endif

        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

#ifdef BIOCFD_MPI
        if (.not. allocated(block(b_blk_no)%ut)) then
            ! If this array isn't allocated we aren't on the right
            ! rank to deal with this so keep going until we find
            ! one that is on this rank
            cycle
        end if
#endif
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

                ! The following block can replace the above once we have an answer on GH issue #116
                ! block(a_blk_no)%wt(tar_x, tar_y, tar_z-1) = trilinear_interpolation( &
                !      block(a_blk_no)%xw(tar_x), block(a_blk_no)%yw(tar_y), block(a_blk_no)%zw(tar_z), &
                !      loc_x, loc_y, loc_z, &
                !      block(b_blk_no)%xw, block(b_blk_no)%yw, block(b_blk_no)%zw, &
                !      3, block(b_blk_no)%wt)



        enddo
        enddo
        enddo
        ENDDO

#ifdef BIOCFD_MPI
        ! We have no finished processing the ut, vt, and wt arrays on
        ! rank n - send these arrays to rank n+1 to continue, unless
        ! we are the last rank, in which case we need will broadcast
        ! the updated array to all ranks (this is maybe something that
        ! can be optimised as we might not need everything to be up to
        ! date here)
        if (rank /= (num_proc-1)) then
          call MPI_Send(block(1)%ut, size(block(1)%ut), MPI_DOUBLE_PRECISION, rank+1, 0, &
                        MPI_COMM_WORLD, ierror)
          call MPI_Send(block(1)%vt, size(block(1)%vt), MPI_DOUBLE_PRECISION, rank+1, 0, &
                        MPI_COMM_WORLD, ierror)
          call MPI_Send(block(1)%wt, size(block(1)%wt), MPI_DOUBLE_PRECISION, rank+1, 0, &
                        MPI_COMM_WORLD, ierror)
        end if

        ! This will block everything until the last rank is done
        call MPI_Bcast(block(1)%ut, size(block(1)%ut), MPI_DOUBLE_PRECISION, num_proc-1, &
                       MPI_COMM_WORLD)
        call MPI_Bcast(block(1)%vt, size(block(1)%vt), MPI_DOUBLE_PRECISION, num_proc-1, &
                       MPI_COMM_WORLD)
        call MPI_Bcast(block(1)%wt, size(block(1)%wt), MPI_DOUBLE_PRECISION, num_proc-1, &
                       MPI_COMM_WORLD)
        ! Now block(1) ut, vt, and wt, should be the same on all ranks
        print *, "Rank = ", rank, "block(1)%ut(1, 1, 1) = ", block(1)%ut(1, 1, 1)
        print *, "Rank = ", rank, "block(1)%vt(1, 1, 1) = ", block(1)%vt(1, 1, 1)
        print *, "Rank = ", rank, "block(1)%wt(1, 1, 1) = ", block(1)%wt(1, 1, 1)
#endif

        end subroutine coarseUpdate_newv
        subroutine coarseUpdate_pc

        REAL (dp) :: bl_interp_ans

        INTEGER :: i,j,k,g, tar_x, tar_y, tar_z, loc_x, loc_y, loc_z, a_blk_no,b_blk_no
        INTEGER :: st_idx, en_idx
        INTEGER :: st_idy, en_idy
        INTEGER :: st_idz, en_idz

#ifdef BIOCFD_MPI
        integer :: rank, num_proc, ierror, token
        call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

        ! See comment in coarseUpdate_newv
        !
        ! We don't wait for anything if we are on rank 0 and just get
        ! started
        if (rank /= 0) then
          ! This essentially acts as a barrier in that rank n wont
          ! start until rank n-1 has finished. Assumes that a_blk number is always 1
          call MPI_Recv(block(1)%pc, size(block(1)%pc), MPI_DOUBLE_PRECISION, rank-1, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
          ! TN: We could consider just doing the copy locally on each rank here
          call MPI_Recv(block(1)%pco, size(block(1)%pco), MPI_DOUBLE_PRECISION, rank-1, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        end if

#endif


        DO g=1,intflines
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk

#ifdef BIOCFD_MPI
        if (.not. allocated(block(b_blk_no)%pc)) then
            ! If this array isn't allocated we aren't on the right
            ! rank to deal with this so keep going until we find
            ! one that is on this rank
            cycle
        end if
#endif

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

#ifdef BIOCFD_MPI
        ! See comment in coarseUpdate_newv
        if (rank /= (num_proc-1)) then
          call MPI_Send(block(1)%pc, size(block(1)%pc), MPI_DOUBLE_PRECISION, rank+1, 0, &
                        MPI_COMM_WORLD, ierror)
          call MPI_Send(block(1)%pco, size(block(1)%pco), MPI_DOUBLE_PRECISION, rank+1, 0, &
                        MPI_COMM_WORLD, ierror)
        end if

        ! This will block everything until the last rank is done
        call MPI_Bcast(block(1)%pc, size(block(1)%pc), MPI_DOUBLE_PRECISION, num_proc-1, &
                       MPI_COMM_WORLD)
        call MPI_Bcast(block(1)%pco, size(block(1)%pco), MPI_DOUBLE_PRECISION, num_proc-1, &
                       MPI_COMM_WORLD)
        ! Now block(1) pc and pco should be the same on all ranks
        print *, "Rank = ", rank, "block(1)%pc(1, 1, 1) = ", block(1)%pc(1, 1, 1)
        print *, "Rank = ", rank, "block(1)%pco(1, 1, 1) = ", block(1)%pco(1, 1, 1)
#endif
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

module biocfd_write_output_corner1
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, ita, totime, re, &
       totime, ita1
  use biocfd_block_type, only: Block_t
#if USE_HDF5 == 1
  use biocfd_hdf5_io, only: hdf5_write_real, hdf5_write_int
#endif
  implicit none

  private

#if USE_HDF5 == 1
  public :: write_output_hdf5, body_plot, writeresult
#else
  public :: write_output_ascii, body_plot, writeresult
#endif

contains

#if USE_HDF5 == 1
      SUBROUTINE write_output_hdf5(blk,blk_no)
       type(Block_t), intent(in) :: blk
       integer (int64), intent(in) :: blk_no
       CHARACTER(len=150)  :: filename
       INTEGER  :: k, i, j
       REAL (dp), allocatable :: u1(:,:,:), v1(:,:,:), w1(:,:,:)

         if (mod(ita,200_int64) /=0 .and. ita > 2) return
            allocate(u1(2:blk%nx+1,2:blk%ny+1,2:blk%nz+1),&
                     v1(2:blk%nx+1,2:blk%ny+1,2:blk%nz+1),&
                     w1(2:blk%nx+1,2:blk%ny+1,2:blk%nz+1))
            DO k = 2, blk%nz+1
            DO j = 2, blk%ny+1
            DO i = 2, blk%nx+1
               u1(i,j,k) = 0.5_dp*(blk%u(i,j,k)+blk%u(i-1,j,k))
               v1(i,j,k) = 0.5_dp*(blk%v(i,j,k)+blk%v(i,j-1,k))
               w1(i,j,k) = 0.5_dp*(blk%w(i,j,k)+blk%w(i,j,k-1))
            END DO
            END DO
            END DO
            write(filename, "('out/timestep_', i5.5, '_block_', i5.5, '.h5')") ita, blk_no
            call hdf5_write_real(filename=filename,&
                                 array_input_3d=u1,key="u1",group="/")
            call hdf5_write_real(filename=filename,&
                                 array_input_3d=v1,key="v1",group="/")
            call hdf5_write_real(filename=filename,&
                                 array_input_3d=w1,key="w1",group="/")
            call hdf5_write_real(filename=filename,&
                                 array_input_1d=blk%xp(2:blk%nx+1), key="xp", group="/")
            call hdf5_write_real(filename=filename,&
                                 array_input_1d=blk%yp(2:blk%ny+1), key="yp",group="/")
            call hdf5_write_real(filename=filename,&
                                 array_input_1d=blk%zp(2:blk%nz+1), key="zp",group="/")
            call hdf5_write_int(filename=filename,&
                                scalar_input=blk%nx,key="zonei",group="/")
            call hdf5_write_int(filename=filename,&
                                scalar_input=blk%nx,key="zonej",group="/")
            call hdf5_write_int(filename=filename,&
                                scalar_input=blk%nz,key="zonek",group="/")
            call hdf5_write_real(filename=filename,&
                                 array_input_3d=blk%p(2:blk%nx+1, 2:blk%ny+1, 2:blk%nz+1), &
                                 key="p",group="/")
            call hdf5_write_real(filename=filename,&
                                 scalar_input=totime,key="totime",group="/")
            call hdf5_write_int(filename=filename,&
                                array_input_3d=blk%cell(2:blk%nx+1, 2:blk%ny+1, 2:blk%nz+1), &
                                key="cell",group="/")
            call hdf5_write_int(filename=filename,&
                                array_input_3d=blk%cell_n(2:blk%nx+1, 2:blk%ny+1, 2:blk%nz+1), &
                                key="cell_n",group="/")
            call hdf5_write_int(filename=filename,&
                                array_input_3d=blk%cell_pr(2:blk%nx+1, 2:blk%ny+1, 2:blk%nz+1), &
                                key="cell_pr",group="/")
            deallocate(u1,v1,w1)

       END SUBROUTINE write_output_hdf5
#else
      SUBROUTINE write_output_ascii(blk,blk_no,char_f)
       type(Block_t), intent(in) :: blk
       integer (int64), intent(in) :: blk_no
       CHARACTER(len=150)  :: filename1
       CHARACTER (LEN = 3),INTENT(IN)   :: char_f
       INTEGER  :: k, i, j
       REAL (dp) :: u1, v1, w1

         IF((mod(ita,200_int64) ==0 .or. ita <= 2))then

       WRITE(filename1,1)char_f,ita,blk_no,re,blk%dx,size(block)
1     FORMAT("out/",A3,"_butter_fielddata.",i9.9,".",i3.3,".",f7.1,".",f8.6,".",i3.3,".dat")
           OPEN(UNIT = 786, FILE = filename1, STATUS = "unknown")
            WRITE(786,*)'variables="x","y","z","u","v","w","p","totime","cellid","cell_n","cell_pr"'
            WRITE(786,*) "zone, ", "i = ", blk%nx," j = ", blk%ny, " k = ", blk%nz

            DO k = 2, blk%nz+1
            DO j = 2, blk%ny+1
            DO i = 2, blk%nx+1
               u1 = 0.5_dp*(blk%u(i,j,k)+blk%u(i-1,j,k))
               v1 = 0.5_dp*(blk%v(i,j,k)+blk%v(i,j-1,k))
               w1 = 0.5_dp*(blk%w(i,j,k)+blk%w(i,j,k-1))
               WRITE(786,*) blk%xp(i),blk%yp(j),blk%zp(k), u1, v1, w1, &
                            blk%p(i,j,k),totime,blk%cell(i,j,k) , &
                            blk%cell_n(i,j,k),blk%cell_pr(i,j,k)
            END DO
            END DO
            END DO
            CLOSE(786)
        !$acc wait
         END IF
      END SUBROUTINE write_output_ascii
#endif
      SUBROUTINE writeResult(blk,blk_no,char_f)
        type(Block_t), intent(in) :: blk
        integer (int64), intent(in) :: blk_no
        INTEGER ::  i, j, k
        CHARACTER(len=70)  :: filename1
        CHARACTER (LEN = 3),INTENT(IN)   :: char_f
        IF(mod(ita,500_int64)/=0) return
           WRITE(filename1,22)char_f,blk_no,re,blk%dx
 22          FORMAT("out/Chkpt/",A3,"_butter_chkpt.",i3.3,".",f6.1,".",f8.6,".dat")
        OPEN (1,FILE=filename1,FORM="formatted")
        DO k = 1, blk%nz+2
        DO j = 1, blk%ny+2
        DO i = 1, blk%nx+2
          WRITE(1,*) blk%u(i,j,k), blk%v(i,j,k), blk%w(i,j,k), &
        blk%p(i,j,k), totime, ita, ita1
       END DO
       END DO
       END DO
        CLOSE(1)

      END SUBROUTINE writeResult

         SUBROUTINE body_plot(blk, blk_id)
         type(Block_t), intent(in) :: blk
         integer(int64), intent(in) :: blk_id
         INTEGER(int64) :: inode, ielem
         CHARACTER(len=150) :: filename1

          if (ita /= 1) return
          WRITE(filename1, "('out/butterfly_', i3.3, '.dat')") blk_id
          OPEN(UNIT=857,FILE=filename1,STATUS="unknown")
          WRITE(857,*) 'TITLE = "FEstressplot"'
          WRITE(857,*) 'VARIABLES= "x", "y", "z","bd_n"'
          WRITE(857,*) "ZONE NODES= ",blk%ibNodes,",ELEMENTS=",blk%ibElems,",DATAPACKING=POINT, ZONETYPE=FETRIANGLE"
          DO inode = 1, blk%ibNodes
            WRITE(857,*) blk%xnode1(inode),blk%ynode1(inode),blk%znode1(inode), 99
          END DO
          DO ielem = 1, blk%ibElems
            WRITE(857,*) blk%ibElP1(ielem),blk%ibElP2(ielem),blk%ibElP3(ielem)
          END DO
          CLOSE(857)
        end subroutine body_plot

end module biocfd_write_output_corner1

module biocfd_write_output_corner1
  use, intrinsic :: iso_fortran_env, only: dp => real64, int64
  use global, only : block, ita, totime, re, nblocks, &
       totime, ita1
  use biocfd_blocks, only: Blocks
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
      SUBROUTINE write_output_hdf5
       CHARACTER(len=150)  :: filename1
       INTEGER  :: k, i, j, g
       REAL (dp), allocatable :: u1(:,:,:), v1(:,:,:), w1(:,:,:)
       character (len=11) :: dummy_1
       character (len=5) ::dummy_2

         IF((mod(ita,200_int64) ==0 .or. ita <= 2 ))then
         do g=1,nblocks
            write(dummy_1,'(A6,I5.5)') 'block_',g
            write(dummy_2,'(I5.5)') ita
            allocate(u1(2:block(g)%nx+1,2:block(g)%ny+1,2:block(g)%nz+1),&
                     v1(2:block(g)%nx+1,2:block(g)%ny+1,2:block(g)%nz+1),&
                     w1(2:block(g)%nx+1,2:block(g)%ny+1,2:block(g)%nz+1))
            DO k = 2, block(g)%nz+1
            DO j = 2, block(g)%ny+1
            DO i = 2, block(g)%nx+1
               u1(i,j,k) = 0.5_dp*(block(g)%u(i,j,k)+block(g)%u(i-1,j,k))
               v1(i,j,k) = 0.5_dp*(block(g)%v(i,j,k)+block(g)%v(i,j-1,k))
               w1(i,j,k) = 0.5_dp*(block(g)%w(i,j,k)+block(g)%w(i,j,k-1))
            END DO
            END DO
            END DO
            filename1="out/output_"//trim(dummy_2)//".h5"
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=u1,key='u1',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=v1,key='v1',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=w1,key='w1',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_1d=block(g)%xp,key='xp',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_1d=block(g)%yp,key='yp',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_1d=block(g)%zp,key='zp',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                scalar_input=block(g)%nx,key='zonei',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                scalar_input=block(g)%nx,key='zonej',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                scalar_input=block(g)%nz,key='zonek',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 array_input_3d=block(g)%p,key='p',group=dummy_1)
            call hdf5_write_real(filename=filename1,&
                                 scalar_input=totime,key='totime',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                array_input_3d=block(g)%cell,key='cell',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                array_input_3d=block(g)%cell_n,key='cell_n',group=dummy_1)
            call hdf5_write_int(filename=filename1,&
                                array_input_3d=block(g)%cell_pr,key='cell_pr',group=dummy_1)
            deallocate(u1,v1,w1)
         end do
        ENDIF
       END SUBROUTINE write_output_hdf5
#else
      SUBROUTINE write_output_ascii(blk,blk_no,char_f)
       type(Blocks), intent(in) :: blk
       integer (int64), intent(in) :: blk_no
       CHARACTER(len=150)  :: filename1
       CHARACTER (LEN = 3),INTENT(IN)   :: char_f
       INTEGER  :: k, i, j
       REAL (dp) :: u1, v1, w1

         IF((mod(ita,200_int64) ==0 .or. ita <= 2 ))then

       WRITE(filename1,1)char_f,ita,blk_no,re,blk%dx,nblocks
1     FORMAT('out/',A3,'_butter_fielddata.',i9.9,'.',i3.3,'.',f7.1,'.',f8.6,'.',i3.3,".dat")
           OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables="x","y","z","u","v","w","p","totime","cellid","cell_n","cell_pr"'
            WRITE(786,*) 'zone, ', 'i = ', blk%nx,' j = ', blk%ny, ' k = ', blk%nz

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
         ENDIF
      END SUBROUTINE write_output_ascii
#endif
      SUBROUTINE writeResult(char_f)
        INTEGER::  i, j, k,g
        CHARACTER(len=70)  :: filename1
        CHARACTER (LEN = 3),INTENT(IN)   :: char_f
        IF(mod(ita,500_int64)==0)THEN
           Do g=1,nblocks
           WRITE(filename1,22)char_f,g,re,block(2)%dx
 22          FORMAT('out/Chkpt/',A3,'_butter_chkpt.',i3.3,'.',f6.1,'.',f8.6,".dat")
        OPEN (1,FILE=filename1,FORM='formatted')
        DO k = 1, block(g)%nz+2
        DO j = 1, block(g)%ny+2
        DO i = 1, block(g)%nx+2
          WRITE(1,*) block(g)%u(i,j,k), block(g)%v(i,j,k), block(g)%w(i,j,k), &
        block(g)%p(i,j,k), totime, ita, ita1
       END DO
       END DO
       END DO
        CLOSE(22)
        END DO
        END IF
      END SUBROUTINE writeResult

         SUBROUTINE body_plot
         INTEGER(int64) :: inode, ielem, g
         CHARACTER(len=150) :: filename1

          DO g=1,nblocks
          if (ita == 1 )then
          WRITE(filename1,108)
  108     FORMAT("out/butterfly.dat")
          OPEN(UNIT=857,FILE=filename1,STATUS='unknown')
          WRITE(857,*) 'TITLE = "FEstressplot"'
          WRITE(857,*) 'VARIABLES= "x", "y", "z","bd_n"'
          WRITE(857,*) 'ZONE NODES= ',block(g)%ibNodes,',ELEMENTS=',block(g)%ibElems,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
          DO inode = 1, block(g)%ibNodes
            WRITE(857,*)block(g)%xnode1(inode),block(g)%ynode1(inode),block(g)%znode1(inode), 99
          END DO
          DO ielem = 1, block(g)%ibElems
            WRITE(857,*) block(g)%ibElP1(ielem), block(g)%ibElP2(ielem), block(g)%ibElP3(ielem)
          END DO
          CLOSE(857)
         END IF
          END DO
        end subroutine body_plot

end module biocfd_write_output_corner1

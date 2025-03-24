      SUBROUTINE allocateArrays
        USE global
        IMPLICIT NONE

        integer(kind=4) :: i
        integer(kind=4) :: nx_var, ny_var, nz_var

        DO i=1,nblocks
        nx_var=block(i)%nx
        ny_var=block(i)%ny
        nz_var=block(i)%nz

        ALLOCATE ( block(i)%u(nx_var+2,ny_var+2,nz_var+2), block(i)%ut(nx_var+2,ny_var+2,nz_var+2), block(i)%v(nx_var+2,ny_var+2,nz_var+2), &
                block(i)%vt(nx_var+2,ny_var+2,nz_var+2), block(i)%w(nx_var+2,ny_var+2,nz_var+2), block(i)%wt(nx_var+2,ny_var+2,nz_var+2),  &
                block(i)%p(nx_var+2,ny_var+2,nz_var+2), block(i)%u_sum(nx_var+2,ny_var+2,nz_var+2), block(i)%v_sum(nx_var+2,ny_var+2,nz_var+2), &
                 block(i)%w_sum(nx_var+2,ny_var+2,nz_var+2), block(i)%p_sum(nx_var+2,ny_var+2,nz_var+2), block(i)%u_avg(nx_var+2,ny_var+2,nz_var+2), &
                block(i)%v_avg(nx_var+2,ny_var+2,nz_var+2),  block(i)%w_avg(nx_var+2,ny_var+2,nz_var+2), block(i)%p_avg(nx_var+2,ny_var+2, nz_var+2), &
	block(i)%resi_u(nx_var+2,ny_var+2,nz_var+2), block(i)%resi_v(nx_var+2,ny_var+2,nz_var+2), block(i)%resi_w(nx_var+2,ny_var+2,nz_var+2) )
        ALLOCATE(block(i)%u_dum( nx_var+2,ny_var+2,nz_var+2), block(i)%v_dum( nx_var+2,ny_var+2,nz_var+2),block(i)%w_dum(nx_var+2,ny_var+2,nz_var+2), block(i)%p_dum(nx_var+2,ny_var+2,nz_var+2))
        ALLOCATE ( block(i)%cell(nx_var+2,ny_var+2,nz_var+2) )
        ALLOCATE ( block(i)%cell2(nx_var+3,ny_var+3,nz_var+3) )
        ALLOCATE ( block(i)%cell_pr(nx_var+3,ny_var+3,nz_var+3) )
        ALLOCATE ( block(i)%nodeIdTag(nx_var+3,ny_var+3,nz_var+3) )
        ALLOCATE ( block(i)%b(nx_var+2,ny_var+2,nz_var+2) )
        ALLOCATE ( block(i)%cell_n(nx_var+2,ny_var+2,nz_var+2) )
        ALLOCATE ( block(i)%Acx(nx_var,3), block(i)%Acy(ny_var,3), block(i)%Acz(nz_var,3) )
        ALLOCATE ( block(i)%pc(nx_var+2,ny_var+2, nz_var+2), block(i)%pco(nx_var+2,ny_var+2, nz_var+2) )
       ALLOCATE ( block(i)%Ac(nx_var*ny_var*nz_var,7) )
       ALLOCATE ( block(i)%An(nx_var*ny_var*nz_var,7) )

        ALLOCATE( block(i)%uv_sum(nx_var+2, ny_var+2, nz_var+2), block(i)%vw_sum(nx_var+2, ny_var+2, nz_var+2), block(i)%uw_sum(nx_var+2,ny_var+2,nz_var +2))
        ALLOCATE( block(i)%uv_avg(nx_var+2, ny_var+2, nz_var+2), block(i)%vw_avg(nx_var+2, ny_var+2, nz_var+2), block(i)%uw_avg(nx_var+2,ny_var+2,nz_var +2))

        ALLOCATE( block(i)%uflu_avg(nx_var+2, ny_var+2, nz_var+2), block(i)%vflu_avg(nx_var+2, ny_var+2, nz_var+2), block(i)%wflu_avg(nx_var+2,ny_var+2,nz_var +2), block(i)%pflu_avg(nx_var+2, ny_var+2, nz_var+2))
        ALLOCATE( block(i)%uflu_rms(nx_var+2, ny_var+2, nz_var+2), block(i)%vflu_rms(nx_var+2, ny_var+2, nz_var+2), block(i)%wflu_rms(nx_var+2,ny_var+2,nz_var +2), block(i)%pflu_rms(nx_var+2, ny_var+2, nz_var+2))
        ALLOCATE( block(i)%u2_sum(nx_var+2, ny_var+2, nz_var+2), block(i)%v2_sum(nx_var+2, ny_var+2, nz_var+2), block(i)%w2_sum(nx_var+2,ny_var+2,nz_var +2),block(i)%p2_sum(nx_var+2,ny_var+2,nz_var +2))
        ALLOCATE( block(i)%u2_avg(nx_var+2, ny_var+2, nz_var+2), block(i)%v2_avg(nx_var+2, ny_var+2, nz_var+2), block(i)%w2_avg(nx_var+2,ny_var+2,nz_var +2),block(i)%p2_avg(nx_var+2,ny_var+2,nz_var +2))
        ALLOCATE( block(i)%ufl(nx_var+2, ny_var+2, nz_var+2), block(i)%vfl(nx_var+2, ny_var+2, nz_var+2), block(i)%wfl(nx_var+2,ny_var+2,nz_var +2))


         END DO
      ! ALLOCATE ( SUMWSS(ibNodes), SIGNWSS(ibNodes,3) )
      ! ALLOCATE ( INSTWSS(ibNodes))
      ! ALLOCATE ( SQSUMWSS(ibNodes))
      ! ALLOCATE ( WSSRMS(ibNodes))
      ! ALLOCATE( TAWSS(ibNodes), OSI(ibNodes), RRT(ibNodes))
      ! ALLOCATE( Afnode(ibNodes,3) ,stressNode(ibNodes,3),Anode(ibNodes,3))



      END SUBROUTINE allocateArrays

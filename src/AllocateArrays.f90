module biocfd_allocate_arrays
    !* Allocate arrays for a number of block variable that are used in
    !  the program. Note that not all of the arrays are allocated
    !  here, several are allocated in other subroutines
    use, intrinsic :: iso_fortran_env, only: int64
    use biocfd_block_type, only: Block_t
    implicit none
    private

    public :: allocateArrays

    contains

    !> Allocate arrays required for various calculations for a block
    SUBROUTINE allocateArrays(blk)

        !> The block that arrays are being allocated for
        type(Block_t), intent(inout) :: blk
        integer(int64) :: nx, ny, nz

        nx = blk%nx
        ny = blk%ny
        nz = blk%nz

        ALLOCATE(&
            blk%u(nx+2, ny+2, nz+2), &
            blk%ut(nx+2, ny+2, nz+2), &
            blk%v(nx+2, ny+2, nz+2), &
            blk%vt(nx+2, ny+2, nz+2), &
            blk%w(nx+2, ny+2, nz+2), &
            blk%wt(nx+2, ny+2, nz+2),  &
            blk%p(nx+2, ny+2, nz+2))

        ALLOCATE(&
            blk%u_dum(nx+2, ny+2, nz+2), &
            blk%v_dum(nx+2, ny+2, nz+2), &
            blk%w_dum(nx+2, ny+2, nz+2), &
            blk%p_dum(nx+2, ny+2, nz+2))

        ALLOCATE(blk%cell(nx+2,ny+2,nz+2))
        ALLOCATE(blk%cell2(nx+3,ny+3,nz+3))
        ALLOCATE(blk%cell_pr(nx+3,ny+3,nz+3))
        ALLOCATE(blk%nodeIdTag(nx+3,ny+3,nz+3))
        ALLOCATE(blk%b(nx+2,ny+2,nz+2))
        ALLOCATE(blk%cell_n(nx+2,ny+2,nz+2))
        ALLOCATE(blk%Acx(nx,3), blk%Acy(ny,3), blk%Acz(nz,3))
        ALLOCATE(blk%pc(nx+2,ny+2, nz+2), blk%pco(nx+2,ny+2, nz+2))

        ALLOCATE(blk%xcent(blk%ibElems), &
                 blk%ycent(blk%ibElems), &
                 blk%zcent(blk%ibElems), &
                 blk%cosAlpha(blk%ibElems), &
                 blk%cosBeta(blk%ibElems), &
                 blk%cosGamma(blk%ibElems), &
                 blk%element_length(blk%ibElems))

      END SUBROUTINE allocateArrays
end module biocfd_allocate_arrays

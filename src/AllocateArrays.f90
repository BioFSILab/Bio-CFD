module biocfd_allocate_arrays
    !* Allocate arrays for a number of block variable that are used in
    !  the program. Note that not all of the arrays are allocated
    !  here, several are allocated in other subroutines
    use, intrinsic :: iso_fortran_env, only: int64
    use global, only : block
    implicit none
    private

    public :: allocateArrays

    contains

    SUBROUTINE allocateArrays

        integer(int64) :: i
        integer(int64) :: nx, ny, nz

        DO i=1, size(block)
            nx = block(i)%nx
            ny = block(i)%ny
            nz = block(i)%nz

            ALLOCATE(&
                block(i)%u(nx+2, ny+2, nz+2), &
                block(i)%ut(nx+2, ny+2, nz+2), &
                block(i)%v(nx+2, ny+2, nz+2), &
                block(i)%vt(nx+2, ny+2, nz+2), &
                block(i)%w(nx+2, ny+2, nz+2), &
                block(i)%wt(nx+2, ny+2, nz+2),  &
                block(i)%p(nx+2, ny+2, nz+2))

            ALLOCATE(&
                block(i)%u_dum(nx+2, ny+2, nz+2), &
                block(i)%v_dum(nx+2, ny+2, nz+2), &
                block(i)%w_dum(nx+2, ny+2, nz+2), &
                block(i)%p_dum(nx+2, ny+2, nz+2))

            ALLOCATE(block(i)%cell(nx+2,ny+2,nz+2))
            ALLOCATE(block(i)%cell2(nx+3,ny+3,nz+3))
            ALLOCATE(block(i)%cell_pr(nx+3,ny+3,nz+3))
            ALLOCATE(block(i)%nodeIdTag(nx+3,ny+3,nz+3))
            ALLOCATE(block(i)%b(nx+2,ny+2,nz+2))
            ALLOCATE(block(i)%cell_n(nx+2,ny+2,nz+2))
            ALLOCATE(block(i)%Acx(nx,3), block(i)%Acy(ny,3), block(i)%Acz(nz,3))
            ALLOCATE(block(i)%pc(nx+2,ny+2, nz+2), block(i)%pco(nx+2,ny+2, nz+2))

        END DO

      END SUBROUTINE allocateArrays
end module biocfd_allocate_arrays

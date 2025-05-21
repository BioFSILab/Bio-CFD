module biocfd_interpolation
        use, intrinsic :: iso_fortran_env, only: dp => real64, int64
        implicit none

        public :: linear_interpolation, bilinear_interpolation

        private

        contains

        ! acc routine directives are needed in the following functions
        ! to ensure that device versions are created

        !> Perform linear interpolation between two points
        !> See https://en.wikipedia.org/wiki/Linear_interpolation
        pure function linear_interpolation(x, x0, x1, y0, y1) result(y)
                !$acc routine
                real(dp), intent(in) :: x, x0, x1, y0, y1
                real(dp) :: y

                y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
        end function linear_interpolation


        !> Perform bilinear interpolation
        !> See https://en.wikipedia.org/wiki/Bilinear_interpolation
        pure function bilinear_interpolation(x, y, x1, x2, y1, y2, fvals) result(fout)
                !$acc routine
                real(dp), intent(in) :: x, y, x1, x2, y1, y2
                !> The values at [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]
                real(dp), intent(in) :: fvals(4)
                real(dp) :: fout

                real(dp) :: numer, denom

                numer = fvals(1) * (x2 - x) * (y2 - y) + fvals(3) * (x2 - x) * (y - y1) &
                      + fvals(2) * (x - x1) * (y2 - y) + fvals(4) * (x - x1) * (y - y1)

                denom = (x2 - x1) * (y2 - y1)

                fout = numer / denom
        end function bilinear_interpolation

end module biocfd_interpolation

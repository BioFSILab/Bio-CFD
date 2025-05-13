module biocfd_interpolation
        use, intrinsic :: iso_fortran_env, only: dp => real64, int64
        implicit none

        public :: linear_interpolation, bilinear_interpolation

        private

        contains

        !> Perform linear interpolation between two points
        !> See https://en.wikipedia.org/wiki/Linear_interpolation
        pure function linear_interpolation(x, x1, x0, y1, y0) result(y)
                real(dp), intent(in) :: x, x1, x0, y1, y0
                real(dp) :: y

                y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
        end function linear_interpolation


        !> Perform bilinear interpolation
        !> See https://en.wikipedia.org/wiki/Bilinear_interpolation
        pure function bilinear_interpolation(x, y, x2, x1, y2, y1, fvals) result(fout)
                real(dp), intent(in) :: x, y, x2, x1, y2, y1
                !> The values at [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]
                real(dp), intent(in) :: fvals(4)
                real(dp) :: fout

                real(dp) :: tmp(2, 1), tmp_fout(1, 1), denom

                tmp = matmul(reshape(fvals, [2, 2]), reshape([y2 - y, y - y1], [2, 1]))
                denom = (x2 - x1) * (y2 - y1)

                tmp_fout = (matmul(reshape([x2 - x, x - x1], [1, 2]), tmp) / denom)
                fout = tmp_fout(1, 1)
        end function bilinear_interpolation

end module biocfd_interpolation

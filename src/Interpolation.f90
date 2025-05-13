module biocfd_interpolation
        use, intrinsic :: iso_fortran_env, only: dp => real64, int64
        implicit none

        public :: linear_interpolation, bilinear_interpolation

        private

        contains

        SUBROUTINE billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1, &
                           bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4, &
                           bl_interp_ans)

        REAL (dp), INTENT(IN) :: bl_intp_valx, bl_intp_valy, bl_intp_x1, bl_intp_x2, &
                                 bl_intp_y1, bl_intp_y2, &
                                 bl_intp_f1, bl_intp_f2, bl_intp_f3, bl_intp_f4
        REAL (dp) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, &
                     bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (dp), INTENT(OUT) :: bl_interp_ans
        bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
        bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
        bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
        bl_intp_yty=(bl_intp_y2-bl_intp_valy)
        bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
        bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
        bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
        bl_intp_num= bl_intp_first_term + bl_intp_second_term

        bl_interp_ans=(bl_intp_num/bl_intp_deno)


        end subroutine billinearInterp


        SUBROUTINE distancePts(d_fine_x1, d_fine_y1, d_coarse_x1, d_coarse_y1, d_val_x1, d_val_y1, &
                               d_fine_frac, d_coarse_frac)


        REAL (dp), INTENT(IN) :: d_fine_x1, d_fine_y1, d_coarse_x1,d_coarse_y1, d_val_x1, d_val_y1
        REAL (dp), INTENT(OUT) :: d_fine_frac,d_coarse_frac
        REAL (dp) :: tot_d, fine_d

        tot_d=SQRT((d_fine_x1-d_coarse_x1)**2 +(d_fine_y1-d_coarse_y1)**2)
        fine_d=SQRT((d_fine_x1-d_val_x1)**2 +(d_fine_y1-d_val_y1)**2)
        d_coarse_frac=fine_d/tot_d
        d_fine_frac=1-d_fine_frac

        end subroutine distancePts

        SUBROUTINE linearInterp(l_intp_valx, l_intp_x1,l_intp_x2, l_intp_y1, l_intp_y2, l_intp_valy)

        REAL (dp) :: varr
        REAL (dp), INTENT(IN) :: l_intp_valx, l_intp_x1, l_intp_x2, l_intp_y1, l_intp_y2
        REAL (dp), INTENT(OUT) :: l_intp_valy
        varr=0.
        varr=(l_intp_valx-l_intp_x1)/(l_intp_x2-l_intp_x1)
        l_intp_valy=(l_intp_y1*(1._dp-varr))+(l_intp_y2*varr)

        end subroutine linearInterp

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

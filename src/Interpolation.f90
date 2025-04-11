module biocfd_interpolation
        use, intrinsic :: iso_fortran_env, only: dp => real64, int64
        implicit none
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

end module biocfd_interpolation

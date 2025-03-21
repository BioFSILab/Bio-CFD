        SUBROUTINE billinearInterp(bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1, &
                bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4, &
                bl_interp_ans)

        REAL (KIND=8), INTENT(IN) :: bl_intp_valx,bl_intp_valy,bl_intp_x1,bl_intp_x2,bl_intp_y1,bl_intp_y2,bl_intp_f1,bl_intp_f2,bl_intp_f3,bl_intp_f4
        REAL (KIND=8) :: bl_intp_deno, bl_intp_num, bl_intp_xtx, bl_intp_xxo, bl_intp_yty, bl_intp_yyo, bl_intp_first_term, bl_intp_second_term
        REAL (KIND=8), INTENT(OUT) :: bl_interp_ans
        bl_intp_deno= (bl_intp_x2-bl_intp_x1) * (bl_intp_y2-bl_intp_y1)
        bl_intp_xtx= (bl_intp_x2 -bl_intp_valx)
        bl_intp_xxo=(bl_intp_valx-bl_intp_x1)
        bl_intp_yty=(bl_intp_y2-bl_intp_valy)
        bl_intp_yyo=(bl_intp_valy-bl_intp_y1)
        bl_intp_first_term=(bl_intp_f1*bl_intp_xtx + bl_intp_f2*bl_intp_xxo)*bl_intp_yty
        bl_intp_second_term=(bl_intp_f4*bl_intp_xtx + bl_intp_f3*bl_intp_xxo)*bl_intp_yyo
        bl_intp_num= bl_intp_first_term + bl_intp_second_term

        bl_interp_ans=(bl_intp_num/bl_intp_deno)


        END SUBROUTINE


        SUBROUTINE distancePts(d_fine_x1, d_fine_y1,d_coarse_x1,d_coarse_y1, d_val_x1, d_val_y1,d_fine_frac,d_coarse_frac)


        REAL (KIND=8), INTENT(IN) :: d_fine_x1, d_fine_y1, d_coarse_x1,d_coarse_y1, d_val_x1, d_val_y1
        REAL (KIND=8), INTENT(OUT) :: d_fine_frac,d_coarse_frac
        REAL (KIND=8) :: tot_d, fine_d

        tot_d=SQRT((d_fine_x1-d_coarse_x1)**2 +(d_fine_y1-d_coarse_y1)**2)
        fine_d=SQRT((d_fine_x1-d_val_x1)**2 +(d_fine_y1-d_val_y1)**2)
        d_coarse_frac=fine_d/tot_d
        d_fine_frac=1-d_fine_frac

        END SUBROUTINE

        SUBROUTINE linearInterp(l_intp_valx, l_intp_x1,l_intp_x2, l_intp_y1, l_intp_y2, l_intp_valy)

        REAL (KIND=8) :: varr
        REAL (KIND=8), INTENT(IN) :: l_intp_valx, l_intp_x1, l_intp_x2, l_intp_y1, l_intp_y2
        REAL (KIND=8), INTENT(OUT) :: l_intp_valy
        varr=0.
        varr=(l_intp_valx-l_intp_x1)/(l_intp_x2-l_intp_x1)
        !print*,'valx',l_intp_valx,'x1',l_intp_x1,'x2',l_intp_x2,'varr',varr,'1-varr',1-varr
        l_intp_valy=(l_intp_y1*(1.-varr))+(l_intp_y2*varr)

        END SUBROUTINE

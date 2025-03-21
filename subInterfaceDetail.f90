        SUBROUTINE interfaceDetail
        USE global
        IMPLICIT NONE
        INTEGER (KIND=8) :: i, j, k, g, f, factor, a_maxx, a_maxy,a_maxz, a_mm,&
        b_maxx, b_maxy, b_maxz, b_mm, a_max_intf_length, b_max_intf_length,  &
        a_blk_no, b_blk_no, xx1_p, xx2_p, yy1_p, yy2_p, zz1_p, zz2_p,counter_coarse,checker,&
        max_yp,max_xp, max_zp,xx1_u, xx2_u, yy1_u, yy2_u, zz1_u,zz2_u,xx1_v, xx2_v, yy1_v, yy2_v,&
        zz1_v,zz2_v, xx1_w,xx2_w,yy1_w,yy2_w,zz1_w,zz2_w, a_mm_x, a_mm_y,a_mm_z
        INTEGER (KIND=8) :: starter, ender
        

       ! do g=1,nblocks
       !   !$ acc enter data copyin(block(g)%xp,block(g)%yp,block(g)%xu,block(g)%yu, block(g)%xv,block(g)%yv)
       !  end do
        
        print*, 'before allocation in interface'
        a_max_intf_length=0 
        b_max_intf_length=0 
        DO g=1, intflines
           a_blk_no=intfr(g)%a_blk
           b_blk_no=intfr(g)%b_blk
           factor=intfr(g)%b_msh/intfr(g)%a_msh
          !a_maxx=FLOOR((intfr(g)%xintf_end-intfr(g)%xintf_start)/block(a_blk_no)%dx) +1
          !a_maxy=FLOOR((intfr(g)%yintf_end-intfr(g)%yintf_start)/block(a_blk_no)%dy) +1
          !a_maxz=FLOOR((intfr(g)%zintf_end-intfr(g)%zintf_start)/block(a_blk_no)%dz) +1
           a_maxx=block(a_blk_no)%nx+3
           a_maxy=block(a_blk_no)%ny+3
           a_maxz=block(a_blk_no)%nz+3
           a_mm=max(a_maxx,a_maxy,a_maxz) 
!          if (a_max_intf_length .lt. a_mm)then 
!               a_max_intf_length=a_mm 
!          end if
           b_maxx=block(b_blk_no)%nx+3
           b_maxy=block(b_blk_no)%ny+3
           b_maxz=block(b_blk_no)%nz+3
          !b_maxx=FLOOR((intfr(g)%xintf_end-intfr(g)%xintf_start)/block(b_blk_no)%dx) +1
          !b_maxy=FLOOR((intfr(g)%yintf_end-intfr(g)%yintf_start)/block(b_blk_no)%dy) +1
          !b_maxz=FLOOR((intfr(g)%zintf_end-intfr(g)%zintf_start)/block(b_blk_no)%dz) +1
           b_mm=max(b_maxx,b_maxy,b_maxz) 
        a_mm_x=int(block(b_blk_no)%nx/factor) +3 
        a_mm_y=int(block(b_blk_no)%ny/factor) +3 
        a_mm_z=int(block(b_blk_no)%nz/factor) +3 
        print*,a_mm_x,a_mm_y,a_mm_z
        print*,block(b_blk_no)%nx, factor
!          if (b_max_intf_length .lt. b_mm)then 
!               b_max_intf_length=b_mm 
!          end if
        a_mm=max(a_mm,b_mm)
        ALLOCATE(intfr(g)%px_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%ux_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%vx_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%wx_interface_det(3,a_mm_x))
        ALLOCATE(intfr(g)%pz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%uz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%vz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%wz_interface_det(3,a_mm_z))
        ALLOCATE(intfr(g)%py_interface_det(3,a_mm_y))
        ALLOCATE(intfr(g)%uy_interface_det(3,a_mm_y))
        ALLOCATE(intfr(g)%vy_interface_det(3,a_mm_y))
        ALLOCATE(intfr(g)%wy_interface_det(3,a_mm_y))
        END DO                                
!       print*,'After first do loop in interface detail' 

!coarse mesh start and end index
        
        print*,'After allocation in interface detail' 
        DO g=1, intflines
        factor=intfr(g)%b_msh/intfr(g)%a_msh
        a_blk_no=intfr(g)%a_blk
        b_blk_no=intfr(g)%b_blk
        DO i=1,block(a_blk_no)%nx+3

        if(block(a_blk_no)%xp(i) .gt. intfr(g)%xintf_start)then

                xx1_p=i
                exit
        endif
        enddo

        !DO i=xx1_p,block(a_blk_no)%itn_en
        DO i=xx1_p,block(a_blk_no)%nx+2

        !if(block(a_blk_no)%xp(i) .gt. intfr(g)%xintf_end+(1.5*block(b_blk_no)%dx))then
        if(block(a_blk_no)%xp(i) .gt. intfr(g)%xintf_end)then

                xx2_p=i-1
                exit
        endif
        enddo
        intfr(g)%px_interface_det(1,1)=xx1_p-1
        intfr(g)%px_interface_det(2,1)=1
        intfr(g)%px_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counterxp=2
12      CONTINUE

        intfr(g)%px_interface_det(1,intfr(g)%counterxp)=xx1_p
        intfr(g)%px_interface_det(2,intfr(g)%counterxp)=starter
        intfr(g)%px_interface_det(3,intfr(g)%counterxp)=ender
        !write(*,91) 'px',starter,ender,block(a_blk_no)%xp(xx1_p),block(b_blk_no)%xp(starter),block(b_blk_no)%xp(ender)
!91      !format(A2,2I3,3F6.3)
        
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterxp=intfr(g)%counterxp+1 
        xx1_p=xx1_p+1
               
        if(xx1_p .le. xx2_p)then
                GOTO 12
        endif

        
        intfr(g)%px_interface_det(1,intfr(g)%counterxp)=xx2_p+1
       !intfr(g)%px_interface_det(2,intfr(g)%counterxp)=starter
       !intfr(g)%px_interface_det(3,intfr(g)%counterxp)=starter+1
        intfr(g)%px_interface_det(2,intfr(g)%counterxp)=starter
        intfr(g)%px_interface_det(3,intfr(g)%counterxp)=starter

        !DO j=1,block(a_blk_no)%jtn_en
        DO j=1,block(a_blk_no)%ny+2

       !!if(block(a_blk_no)%yp(j) .gt. intfr(g)%yintf_start-(1.5*block(b_blk_no)%dy))then
        if(block(a_blk_no)%yp(j) .gt. intfr(g)%yintf_start)then
         !       print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_start-(2*block(b_blk_no)%dy)
                yy1_p=j
                exit
        endif
        enddo

        !DO j=yy1_p,block(a_blk_no)%jtn_en
        DO j=yy1_p,block(a_blk_no)%ny+2

        !if(block(a_blk_no)%y1(j) .gt. intfr(g)%yintf_end+(1.5*block(b_blk_no)%dy))then
        if(block(a_blk_no)%yp(j) .gt. intfr(g)%yintf_end)then
          !      print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_end+(2*block(b_blk_no)%dy)

                yy2_p=j-1
                exit
        endif
        enddo
        intfr(g)%py_interface_det(1,1)=yy1_p-1
        intfr(g)%py_interface_det(2,1)=1
        intfr(g)%py_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counteryp=2
13      CONTINUE

        intfr(g)%py_interface_det(1,intfr(g)%counteryp)=yy1_p
        intfr(g)%py_interface_det(2,intfr(g)%counteryp)=starter
        intfr(g)%py_interface_det(3,intfr(g)%counteryp)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counteryp=intfr(g)%counteryp+1 
        yy1_p=yy1_p+1
               
        if(yy1_p .le. yy2_p)then
                GOTO 13
        endif

        
        intfr(g)%py_interface_det(1,intfr(g)%counteryp)=yy2_p+1
        intfr(g)%py_interface_det(2,intfr(g)%counteryp)=starter
        intfr(g)%py_interface_det(3,intfr(g)%counteryp)=starter

        DO j=1,block(a_blk_no)%nz+2

       !!if(block(a_blk_no)%yp(j) .gt. intfr(g)%yintf_start-(1.5*block(b_blk_no)%dy))then
        if(block(a_blk_no)%zp(j) .gt. intfr(g)%zintf_start)then
         !       print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_start-(2*block(b_blk_no)%dy)
                zz1_p=j
                exit
        endif
        enddo

        !DO j=yy1_p,block(a_blk_no)%jtn_en
        DO j=zz1_p,block(a_blk_no)%nz+2

        !if(block(a_blk_no)%y1(j) .gt. intfr(g)%yintf_end+(1.5*block(b_blk_no)%dy))then
        if(block(a_blk_no)%zp(j) .gt. intfr(g)%zintf_end)then
          !      print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_end+(2*block(b_blk_no)%dy)

                zz2_p=j-1
                exit
        endif
        enddo
        intfr(g)%pz_interface_det(1,1)=zz1_p-1
        intfr(g)%pz_interface_det(2,1)=1
        intfr(g)%pz_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counterzp=2
913      CONTINUE

        intfr(g)%pz_interface_det(1,intfr(g)%counterzp)=zz1_p
        intfr(g)%pz_interface_det(2,intfr(g)%counterzp)=starter
        intfr(g)%pz_interface_det(3,intfr(g)%counterzp)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterzp=intfr(g)%counterzp+1 
        zz1_p=zz1_p+1
               
        if(zz1_p .le. zz2_p)then
                GOTO 913
        endif

        
        intfr(g)%pz_interface_det(1,intfr(g)%counterzp)=zz2_p+1
        intfr(g)%pz_interface_det(2,intfr(g)%counterzp)=starter
        intfr(g)%pz_interface_det(3,intfr(g)%counterzp)=starter
!!!!!!!!!!!!!!!!!!!uuuuuuuuuuuuu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !DO i=1,block(a_blk_no)%itn_en
        DO i=1,block(a_blk_no)%nx+3

        if(block(a_blk_no)%xu(i) .ge. intfr(g)%xintf_start)then

                xx1_p=i
                exit
        endif
        enddo

        !DO i=xx1_p,block(a_blk_no)%itn_en
        DO i=xx1_p,block(a_blk_no)%nx+3

        !if(block(a_blk_no)%xu(i) .gt. intfr(g)%xintf_end+(2*block(b_blk_no)%dx))then
        if(block(a_blk_no)%xu(i) .gt. intfr(g)%xintf_end)then

                xx2_p=i-1
                exit
        endif
        enddo
        intfr(g)%ux_interface_det(1,1)=xx1_p
        intfr(g)%ux_interface_det(2,1)=2
        intfr(g)%ux_interface_det(3,1)=2
        starter=2+1
        ender=starter+factor-1
        intfr(g)%counterxu=2
        print*,block(a_blk_no)%xu(xx1_p),block(a_blk_no)%xu(xx2_p),block(b_blk_no)%xu(2)
121      CONTINUE

        xx1_p=xx1_p+1
        intfr(g)%ux_interface_det(1,intfr(g)%counterxu)=xx1_p
        intfr(g)%ux_interface_det(2,intfr(g)%counterxu)=starter
        intfr(g)%ux_interface_det(3,intfr(g)%counterxu)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterxu=intfr(g)%counterxu+1 
               
        if(xx1_p .lt. xx2_p)then
                GOTO 121
        endif

        
        intfr(g)%ux_interface_det(1,intfr(g)%counterxu)=xx2_p+1
        intfr(g)%ux_interface_det(2,intfr(g)%counterxu)=starter
        intfr(g)%ux_interface_det(3,intfr(g)%counterxu)=starter

        !DO j=1,block(a_blk_no)%jtn_en
        DO j=1,block(a_blk_no)%ny+2

        if(block(a_blk_no)%yu(j) .gt. intfr(g)%yintf_start)then
                !print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_start-(2*block(b_blk_no)%dy)
                yy1_p=j
                exit
        endif
        enddo

        !DO j=yy1_p,block(a_blk_no)%jtn_en
        DO j=yy1_p,block(a_blk_no)%ny+2

        !if(block(a_blk_no)%yu(j) .gt. intfr(g)%yintf_end+(1.5*block(b_blk_no)%dy))then
        if(block(a_blk_no)%yu(j) .gt. intfr(g)%yintf_end)then
        !        print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_end+(2*block(b_blk_no)%dy)

                yy2_p=j-1
                exit
        endif
        enddo
        intfr(g)%uy_interface_det(1,1)=yy1_p-1
        intfr(g)%uy_interface_det(2,1)=1
        intfr(g)%uy_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counteryu=2
131      CONTINUE

        intfr(g)%uy_interface_det(1,intfr(g)%counteryu)=yy1_p
        intfr(g)%uy_interface_det(2,intfr(g)%counteryu)=starter
        intfr(g)%uy_interface_det(3,intfr(g)%counteryu)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counteryu=intfr(g)%counteryu+1 
        yy1_p=yy1_p+1
               
        if(yy1_p .le. yy2_p)then
                GOTO 131
        endif

        
        intfr(g)%uy_interface_det(1,intfr(g)%counteryu)=yy2_p+1
        intfr(g)%uy_interface_det(2,intfr(g)%counteryu)=starter
        intfr(g)%uy_interface_det(3,intfr(g)%counteryu)=starter
        


        DO j=1,block(a_blk_no)%nz+2

        if(block(a_blk_no)%zu(j) .gt. intfr(g)%zintf_start)then
                !print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_start-(2*block(b_blk_no)%dy)
                zz1_p=j
                exit
        endif
        enddo

        !DO j=yy1_p,block(a_blk_no)%jtn_en
        DO j=zz1_p,block(a_blk_no)%nz+2

        !if(block(a_blk_no)%yu(j) .gt. intfr(g)%yintf_end+(1.5*block(b_blk_no)%dy))then
        if(block(a_blk_no)%zu(j) .gt. intfr(g)%zintf_end)then
        !        print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_end+(2*block(b_blk_no)%dy)

                zz2_p=j-1
                exit
        endif
        enddo
        intfr(g)%uz_interface_det(1,1)=zz1_p-1
        intfr(g)%uz_interface_det(2,1)=1
        intfr(g)%uz_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counterzu=2
9131      CONTINUE

        intfr(g)%uz_interface_det(1,intfr(g)%counterzu)=zz1_p
        intfr(g)%uz_interface_det(2,intfr(g)%counterzu)=starter
        intfr(g)%uz_interface_det(3,intfr(g)%counterzu)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterzu=intfr(g)%counterzu+1 
        zz1_p=zz1_p+1
               
        if(zz1_p .le. zz2_p)then
                GOTO 9131
        endif

        
        intfr(g)%uz_interface_det(1,intfr(g)%counterzu)=zz2_p+1
        intfr(g)%uz_interface_det(2,intfr(g)%counterzu)=starter
        intfr(g)%uz_interface_det(3,intfr(g)%counterzu)=starter
!!!!!!!!v
        !DO i=1,block(a_blk_no)%itn_en
        DO i=1,block(a_blk_no)%nx+2

        if(block(a_blk_no)%xv(i) .gt. intfr(g)%xintf_start)then

                xx1_p=i
                exit
        endif
        enddo

        !DO i=xx1_p,block(a_blk_no)%itn_en
        DO i=xx1_p,block(a_blk_no)%nx+2

        !if(block(a_blk_no)%xv(i) .gt. intfr(g)%xintf_end+(1.5*block(b_blk_no)%dx))then
        if(block(a_blk_no)%xv(i) .gt. intfr(g)%xintf_end)then

                xx2_p=i-1
                exit
        endif
        enddo
        intfr(g)%vx_interface_det(1,1)=xx1_p-1
        intfr(g)%vx_interface_det(2,1)=1
        intfr(g)%vx_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counterxv=2
122      CONTINUE

        intfr(g)%vx_interface_det(1,intfr(g)%counterxv)=xx1_p
        intfr(g)%vx_interface_det(2,intfr(g)%counterxv)=starter
        intfr(g)%vx_interface_det(3,intfr(g)%counterxv)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterxv=intfr(g)%counterxv+1 
        xx1_p=xx1_p+1
               
        if(xx1_p .le. xx2_p)then
                GOTO 122
        endif

        
        intfr(g)%vx_interface_det(1,intfr(g)%counterxv)=xx2_p+1
        intfr(g)%vx_interface_det(2,intfr(g)%counterxv)=starter
        intfr(g)%vx_interface_det(3,intfr(g)%counterxv)=starter

        !DO j=1,block(a_blk_no)%jtn_en
        DO j=1,block(a_blk_no)%ny+3

        if(block(a_blk_no)%yv(j) .ge. intfr(g)%yintf_start)then
                !print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_start-(2*block(b_blk_no)%dy)
                yy1_p=j
                exit
        endif
        enddo

        !DO j=yy1_p,block(a_blk_no)%jtn_en
        DO j=yy1_p,block(a_blk_no)%ny+3

        if(block(a_blk_no)%yv(j) .gt. intfr(g)%yintf_end)then
        !if(block(a_blk_no)%yv(j) .gt. intfr(g)%yintf_end+(2*block(b_blk_no)%dy))then
                !print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_end+(2*block(b_blk_no)%dy)

                yy2_p=j-1
                exit
        endif
        enddo
        intfr(g)%vy_interface_det(1,1)=yy1_p
        intfr(g)%vy_interface_det(2,1)=2
        intfr(g)%vy_interface_det(3,1)=2
        starter=2+1
        ender=starter+factor-1
        intfr(g)%counteryv=2
132      CONTINUE
        yy1_p=yy1_p+1

        intfr(g)%vy_interface_det(1,intfr(g)%counteryv)=yy1_p
        intfr(g)%vy_interface_det(2,intfr(g)%counteryv)=starter
        intfr(g)%vy_interface_det(3,intfr(g)%counteryv)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counteryv=intfr(g)%counteryv+1 
               
        if(yy1_p .lt. yy2_p)then
                GOTO 132
        endif

        
        intfr(g)%vy_interface_det(1,intfr(g)%counteryv)=yy2_p+1
        intfr(g)%vy_interface_det(2,intfr(g)%counteryv)=starter
        intfr(g)%vy_interface_det(3,intfr(g)%counteryv)=starter
        

        DO i=1,block(a_blk_no)%nz+2

        if(block(a_blk_no)%zv(i) .gt. intfr(g)%zintf_start)then

                zz1_p=i
                exit
        endif
        enddo

        !DO i=xx1_p,block(a_blk_no)%itn_en
        DO i=zz1_p,block(a_blk_no)%nz+2

        !if(block(a_blk_no)%xv(i) .gt. intfr(g)%xintf_end+(1.5*block(b_blk_no)%dx))then
        if(block(a_blk_no)%zv(i) .gt. intfr(g)%zintf_end)then

                zz2_p=i-1
                exit
        endif
        enddo
        intfr(g)%vz_interface_det(1,1)=zz1_p-1
        intfr(g)%vz_interface_det(2,1)=1
        intfr(g)%vz_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counterzv=2
9122      CONTINUE

        intfr(g)%vz_interface_det(1,intfr(g)%counterzv)=zz1_p
        intfr(g)%vz_interface_det(2,intfr(g)%counterzv)=starter
        intfr(g)%vz_interface_det(3,intfr(g)%counterzv)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterzv=intfr(g)%counterzv+1 
        zz1_p=zz1_p+1
               
        if(zz1_p .le. zz2_p)then
                GOTO 9122
        endif

        
        intfr(g)%vz_interface_det(1,intfr(g)%counterzv)=zz2_p+1
        intfr(g)%vz_interface_det(2,intfr(g)%counterzv)=starter
        intfr(g)%vz_interface_det(3,intfr(g)%counterzv)=starter

!!!!!!!!w!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !DO i=1,block(a_blk_no)%itn_en
        DO i=1,block(a_blk_no)%nx+2

        if(block(a_blk_no)%xw(i) .gt. intfr(g)%xintf_start)then

                xx1_p=i
                exit
        endif
        enddo

        !DO i=xx1_p,block(a_blk_no)%itn_en
        DO i=xx1_p,block(a_blk_no)%nx+2

        !if(block(a_blk_no)%xv(i) .gt. intfr(g)%xintf_end+(1.5*block(b_blk_no)%dx))then
        if(block(a_blk_no)%xw(i) .gt. intfr(g)%xintf_end)then

                xx2_p=i-1
                exit
        endif
        enddo
        intfr(g)%wx_interface_det(1,1)=xx1_p-1
        intfr(g)%wx_interface_det(2,1)=1
        intfr(g)%wx_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counterxw=2
1229      CONTINUE

        intfr(g)%wx_interface_det(1,intfr(g)%counterxw)=xx1_p
        intfr(g)%wx_interface_det(2,intfr(g)%counterxw)=starter
        intfr(g)%wx_interface_det(3,intfr(g)%counterxw)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterxw=intfr(g)%counterxw+1 
        xx1_p=xx1_p+1
               
        if(xx1_p .le. xx2_p)then
                GOTO 1229
        endif

        
        intfr(g)%wx_interface_det(1,intfr(g)%counterxw)=xx2_p+1
        intfr(g)%wx_interface_det(2,intfr(g)%counterxw)=starter
        intfr(g)%wx_interface_det(3,intfr(g)%counterxw)=starter

        !DO j=1,block(a_blk_no)%jtn_en
        DO j=1,block(a_blk_no)%ny+2

        if(block(a_blk_no)%yw(j) .gt. intfr(g)%yintf_start)then
                !print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_start-(2*block(b_blk_no)%dy)
                yy1_p=j
                exit
        endif
        enddo

        !DO j=yy1_p,block(a_blk_no)%jtn_en
        DO j=yy1_p,block(a_blk_no)%ny+2

        if(block(a_blk_no)%yw(j) .gt. intfr(g)%yintf_end)then
        !if(block(a_blk_no)%yv(j) .gt. intfr(g)%yintf_end+(2*block(b_blk_no)%dy))then
                !print*,block(a_blk_no)%y1(j),j,intfr(g)%yintf_end+(2*block(b_blk_no)%dy)

                yy2_p=j-1
                exit
        endif
        enddo
        intfr(g)%wy_interface_det(1,1)=yy1_p-1
        intfr(g)%wy_interface_det(2,1)=1
        intfr(g)%wy_interface_det(3,1)=1
        starter=1+1
        ender=starter+factor-1
        intfr(g)%counteryw=2
1329      CONTINUE

        intfr(g)%wy_interface_det(1,intfr(g)%counteryw)=yy1_p
        intfr(g)%wy_interface_det(2,intfr(g)%counteryw)=starter
        intfr(g)%wy_interface_det(3,intfr(g)%counteryw)=ender
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counteryw=intfr(g)%counteryw+1 
        yy1_p=yy1_p+1
               
        if(yy1_p .le. yy2_p)then
                GOTO 1329
        endif

        
        intfr(g)%wy_interface_det(1,intfr(g)%counteryw)=yy2_p+1
        intfr(g)%wy_interface_det(2,intfr(g)%counteryw)=starter
        intfr(g)%wy_interface_det(3,intfr(g)%counteryw)=starter
        

        DO i=1,block(a_blk_no)%nz+3

        if(block(a_blk_no)%zw(i) .ge. intfr(g)%zintf_start)then

                zz1_p=i
                exit
        endif
        enddo

        !DO i=xx1_p,block(a_blk_no)%itn_en
        DO i=zz1_p,block(a_blk_no)%nz+3

        !if(block(a_blk_no)%xv(i) .gt. intfr(g)%xintf_end+(1.5*block(b_blk_no)%dx))then
        if(block(a_blk_no)%zw(i) .gt. intfr(g)%zintf_end)then

                zz2_p=i-1
                exit
        endif
        enddo
        intfr(g)%wz_interface_det(1,1)=zz1_p
        intfr(g)%wz_interface_det(2,1)=2
        intfr(g)%wz_interface_det(3,1)=2
        starter=2+1
        ender=starter+factor-1
        intfr(g)%counterzw=2
91229      CONTINUE

        zz1_p=zz1_p+1
        intfr(g)%wz_interface_det(1,intfr(g)%counterzw)=zz1_p
        intfr(g)%wz_interface_det(2,intfr(g)%counterzw)=starter
        intfr(g)%wz_interface_det(3,intfr(g)%counterzw)=ender
!        write(*,91) 'wz',starter,ender,block(a_blk_no)%zw(xx1_p),block(b_blk_no)%zwp(starter),block(b_blk_no)%zw(ender)
!91     format(A2,2I3,3F6.3)
        starter=ender+1        
        ender=ender+factor       
        intfr(g)%counterzw=intfr(g)%counterzw+1 
               
        if(zz1_p .lt. zz2_p)then
                GOTO 91229
        endif

        
        intfr(g)%wz_interface_det(1,intfr(g)%counterzw)=zz2_p+1
        intfr(g)%wz_interface_det(2,intfr(g)%counterzw)=starter
        intfr(g)%wz_interface_det(3,intfr(g)%counterzw)=starter

        END DO

        DO j=1,intflines
        print*,'**********************px***************************'
        DO i=1,intfr(j)%counterxp
        WRITE(*,33)'px',i,intfr(j)%px_interface_det(1,i),intfr(j)%px_interface_det(2,i),intfr(j)%px_interface_det(3,i),block(a_blk_no)%xp(intfr(j)%px_interface_det(1,i)),block(b_blk_no)%xp(intfr(j)%px_interface_det(2,i)),block(b_blk_no)%xp(intfr(j)%px_interface_det(3,i))
 33       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************py***************************'
        DO i=1,intfr(j)%counteryp
        WRITE(*,133)'py',i,intfr(j)%py_interface_det(1,i),intfr(j)%py_interface_det(2,i),intfr(j)%py_interface_det(3,i),block(a_blk_no)%yp(intfr(j)%py_interface_det(1,i)),block(b_blk_no)%yp(intfr(j)%py_interface_det(2,i)),block(b_blk_no)%yp(intfr(j)%py_interface_det(3,i))
 133       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do
        DO j=1,intflines
        print*,'**********************pz***************************'
        DO i=1,intfr(j)%counterzp
        WRITE(*,933)'pz',i,intfr(j)%pz_interface_det(1,i),intfr(j)%pz_interface_det(2,i),intfr(j)%pz_interface_det(3,i),block(a_blk_no)%zp(intfr(j)%pz_interface_det(1,i)),block(b_blk_no)%zp(intfr(j)%pz_interface_det(2,i)),block(b_blk_no)%zp(intfr(j)%pz_interface_det(3,i))
 933       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do
        
        DO j=1,intflines
        print*,'**********************ux***************************'
        DO i=1,intfr(j)%counterxu
        WRITE(*,331)'ux',i,intfr(j)%ux_interface_det(1,i),intfr(j)%ux_interface_det(2,i),intfr(j)%ux_interface_det(3,i)  ,block(a_blk_no)%xu(intfr(j)%ux_interface_det(1,i)),block(b_blk_no)%xu(intfr(j)%ux_interface_det(2,i)),block(b_blk_no)%xu(intfr(j)%ux_interface_det(3,i))
 331       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************uy***************************'
        DO i=1,intfr(j)%counteryu
        WRITE(*,1331)'uy',i,intfr(j)%uy_interface_det(1,i),intfr(j)%uy_interface_det(2,i),intfr(j)%uy_interface_det(3,i) ,block(a_blk_no)%yu(intfr(j)%uy_interface_det(1,i)),block(b_blk_no)%yu(intfr(j)%uy_interface_det(2,i)),block(b_blk_no)%yu(intfr(j)%uy_interface_det(3,i))
 1331       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do
        DO j=1,intfLines
        print*,'**********************uz***************************'
        DO i=1,intfr(j)%counterzu
        WRITE(*,91331)'uz',i,intfr(j)%uz_interface_det(1,i),intfr(j)%uz_interface_det(2,i),intfr(j)%uz_interface_det(3,i) ,block(a_blk_no)%zu(intfr(j)%uz_interface_det(1,i)),block(b_blk_no)%zu(intfr(j)%uz_interface_det(2,i)),block(b_blk_no)%zu(intfr(j)%uz_interface_det(3,i))
91331       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do

        DO j=1,intflines
        print*,'**********************vx***************************'
        DO i=1,intfr(j)%counterxv
        WRITE(*,332)'vx',i,intfr(j)%vx_interface_det(1,i),intfr(j)%vx_interface_det(2,i),intfr(j)%vx_interface_det(3,i)  ,block(a_blk_no)%xv(intfr(j)%vx_interface_det(1,i)),block(b_blk_no)%xv(intfr(j)%vx_interface_det(2,i)),block(b_blk_no)%xv(intfr(j)%vx_interface_det(3,i))
 332       FORMAT(A3,I5,3I5,3F10.5)                                                                                                                                                                                                                                                      
        end do                                                                                                                                                                                                                                                                    
        end do                                                                                                                                                                                                                                                                    
        DO j=1,intfLines                                                                                                                                                                                                                                                          
        print*,'**********************vy***************************'                                                                                                                                                                                                              
        DO i=1,intfr(j)%counteryv                                                                                                                                                                                                                                                 
        WRITE(*,1332)'vy',i,intfr(j)%vy_interface_det(1,i),intfr(j)%vy_interface_det(2,i),intfr(j)%vy_interface_det(3,i) ,block(a_blk_no)%yv(intfr(j)%vy_interface_det(1,i)),block(b_blk_no)%yv(intfr(j)%vy_interface_det(2,i)),block(b_blk_no)%yv(intfr(j)%vy_interface_det(3,i))
 1332       FORMAT(A3,I5,3I5,3F10.5)                                                                                                                                                                                                                                                     
        end do                                                                                                                                                                                                                                                                    
        end do                                                                                                                                                                                                                                                                    
        DO j=1,intfLines                                                                                                                                                                                                                                                          
        print*,'**********************vz***************************'                                                                                                                                                                                                              
        DO i=1,intfr(j)%counterzv                                                                                                                                                                                                                                                 
        WRITE(*,91332)'vz',i,intfr(j)%vz_interface_det(1,i),intfr(j)%vz_interface_det(2,i),intfr(j)%vz_interface_det(3,i) ,block(a_blk_no)%zv(intfr(j)%vz_interface_det(1,i)),block(b_blk_no)%zv(intfr(j)%vz_interface_det(2,i)),block(b_blk_no)%zv(intfr(j)%vz_interface_det(3,i))
91332       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do


        DO j=1,intflines
        print*,'**********************wx***************************'
        DO i=1,intfr(j)%counterxw
        WRITE(*,3329)'wx',i,intfr(j)%wx_interface_det(1,i),intfr(j)%wx_interface_det(2,i),intfr(j)%wx_interface_det(3,i),block(a_blk_no)%xw(intfr(j)%wx_interface_det(1,i)),block(b_blk_no)%xw(intfr(j)%wx_interface_det(2,i)),block(b_blk_no)%xw(intfr(j)%wx_interface_det(3,i))
 3329       FORMAT(A3,I5,3I5,3F10.5)                                                                                                                                                                                                                                                    
        end do                                                                                                                                                                                                                                                                   
        end do                                                                                                                                                                                                                                                                   
        DO j=1,intfLines                                                                                                                                                                                                                                                         
        print*,'**********************wy***************************'                                                                                                                                                                                                             
        DO i=1,intfr(j)%counteryw                                                                                                                                                                                                                                                
        WRITE(*,13329)'wy',i,intfr(j)%wy_interface_det(1,i),intfr(j)%wy_interface_det(2,i),intfr(j)%wy_interface_det(3,i),block(a_blk_no)%yw(intfr(j)%wy_interface_det(1,i)),block(b_blk_no)%yw(intfr(j)%wy_interface_det(2,i)),block(b_blk_no)%yw(intfr(j)%wy_interface_det(3,i)) 
13329       FORMAT(A3,I5,3I5,3F10.5)                                                                                                                                                                                                                                                    
        end do                                                                                                                                                                                                                                                                   
        end do                                                                                                                                                                                                                                                                   
        DO j=1,intfLines                                                                                                                                                                                                                                                         
        print*,'**********************wz***************************'                                                                                                                                                                                                             
        DO i=1,intfr(j)%counterzw                                                                                                                                                                                                                                                
        WRITE(*,93329)'wz',i,intfr(j)%wz_interface_det(1,i),intfr(j)%wz_interface_det(2,i),intfr(j)%wz_interface_det(3,i),block(a_blk_no)%zw(intfr(j)%wz_interface_det(1,i)),block(b_blk_no)%zw(intfr(j)%wz_interface_det(2,i)),block(b_blk_no)%zw(intfr(j)%wz_interface_det(3,i))
93329       FORMAT(A3,I5,3I5,3F10.5)
        end do
        end do




        END SUBROUTINE


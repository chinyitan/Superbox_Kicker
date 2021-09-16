c**********************************************************************
c
c     P U S H E R 1 . F
c
c**********************************************************************
c!! Y ahora el mas dificil todavia, meter la formula de Binney (1977)
c!! (el halo esta fijo==>su self-potential)
c!! Que Dios nos pille confesaos

c
      subroutine pusher1 (gnr,step)
c
      include 'super.cb'
c
c     this routine updates the velocities of all stars of galaxy
c     gnr by differentiating the potential numerically in second
c     order and integrating the accelerations via Leap-frog-scheme
c
c     the variables used in this routine are changed to double
c     precision because of the small quantities computed here
c
      double precision rout1,rcore1,hepot
      double precision sf1,sf2,sf3,sf4,sf5
      double precision sf1q,sf2q,sf3q,sf4q,sf5q
      double precision xcms,ycms,zcms
      double precision ex,ey,ez,fxx,fxy,fxz,fyy,fyz,fzz
      double precision x,y,z,dx,dy,dz
      double precision x11,x12,x13,x14,x15
      double precision y11,y12,y13,y14,y15
      double precision z11,z12,z13,z14,z15
      double precision r1,fnh,dx15,dy15,dz15,fx,fy,fz
      double precision dx14,dy14,dz14
      double precision dx11,dy11,dz11
      double precision dx12,dy12,dz12
      double precision dx13,dy13,dz13
      double precision dbvx,dbvy,dbvz
      integer          start,istop,gnr,step
      integer          i,j,k
      integer          ix11,iy11,iz11
      integer          ix12,iy12,iz12
      integer          ix13,iy13,iz13
      integer          ix14,iy14,iz14
      integer          ix15,iy15,iz15
      real             rv(6),tv(6),tekin,tlx,tly,tlz
c!!
      double precision adfx,adfy,adfz,rcm,vc,r,e,delta,sigt,sigp,G
      real Mgal,Ms,rcut,gama
      integer nbod,numero
      double precision ucms,vcms,wcms
      parameter (G=2.)
      parameter(e=0.6)
      delta=e

c
c     Calculating the quadratic radii of the grids :
c     ==============================================
c
      rout1  = dble(rout(gnr)*rout(gnr))
      rcore1 = dble(rcore(gnr)*rcore(gnr))
      fnh    = dble(n/2)
c
c     Changing the center of density to double precision :
c     ====================================================
c
      xcms = dble(dgcms(1,gnr))
      ycms = dble(dgcms(2,gnr))
      zcms = dble(dgcms(3,gnr))
c!!
      ucms=dble(dgcms(4,gnr))
      vcms=dble(dgcms(5,gnr))
      wcms=dble(dgcms(6,gnr))
c!!


c
c     Changing the expanding factors of the grids to d.p. :
c     =====================================================
c
      sf1  = dble(enh1(gnr))
      sf2  = dble(enh2(gnr))
      sf3  = dble(enh3(gnr))
      sf4  = dble(enh4(gnr))
      sf5  = dble(enh5(gnr))
      sf1q = sf1**2
      sf2q = sf2**2
      sf3q = sf3**2
      sf4q = sf4**2
      sf5q = sf5**2

      if (gnr.gt.1) then
c!!   Calculo vrot gnum=1
         Mgal=(fh(5,1))
         Ms=fh(5,gnr)
         rcut=fh(2,1)
c!!   -----------------------
         gama=1.
c!!   -----------------------

         rcm=sqrt(xcms**2.+ycms**2.+zcms**2.)

         numero=0
         do i=1,gstno(1)
            r=sqrt(dble(star(1,i))**2.+dble(star(2,i))**2.+
     &           dble(star(3,i))**2.)
            if (r.le.rcm) numero=numero+1
         enddo
         vc=sqrt(G*dble(numero)*Mgal/gstno(1)/rcm)
         sigt=sqrt(3.*vc**2./((3.-delta)*2.))
         sigp=sqrt(3.*(1.-delta)*vc**2./((3.-delta)*2.))

         PRINT*, gnr,vc,sigt,sigp,e,Mgal,Ms,rcut,gama,delta
         PRINT*, xcms,ycms,zcms,ucms,vcms,wcms

      
c!!   calculo de la df, por ahora Binney
         call binney(vc,sigt,sigp,e,Mgal,Ms,rcut,gama,xcms,ycms,zcms,
     &        ucms,vcms,wcms,adfx,adfy,adfz)
c!!

      endif
c
c     Starting the loop over all galaxies :
c     =====================================
c
      do k = 1,gnum             ! starting galaxy loop
c
         if (k .eq. 1) then     ! finding the number of the first
            start = 1           ! and the last star of galaxy k
            istop  = gstno(1)
         else
            start = 0
            do i = 1,k-1
               start = start + gstno(i)
            enddo
            istop  = start + gstno(k)
            start  = start + 1
         endif
c
         hepot  = 0.d0        ! initialize potential energy
         tekin  = 0.0
         tlx    = 0.0
         tly    = 0.0
         tlz    = 0.0
c
c     Starting the particle loop of galaxy k :
c     ========================================
c
         do  i = start,istop  ! starting particle loop
c
c     Skip stars outside the simulation :
c     ===================================
c
            if (star(1,i) .gt. 1.e+10) goto 20 
c
c     Calculating the positions of the star for each grid :
c     =====================================================
c
            x    = dble(star(1,i))
            y    = dble(star(2,i))
            z    = dble(star(3,i))
            dx   = (x - xcms)
            dy   = (y - ycms)
            dz   = (z - zcms)
            x11  = sf1*dx + fnh
            x12  = x11
            x13  = sf3*dx + fnh
            x14  = sf4* x + fnh
            x15  = x14
            y11  = sf1*dy + fnh
            y12  = y11
            y13  = sf3*dy + fnh
            y14  = sf4* y + fnh
            y15  = y14
            z11  = sf1*dz + fnh
            z12  = z11
            z13  = sf3*dz + fnh
            z14  = sf4* z + fnh
            z15  = z14
c            
c     Quadratic radius of star i to the focus of grids :
c     ==================================================
c
            r1   = dx*dx + dy*dy + dz*dz
            
c
c     (1.) For all: phi 15
c     ====================
c
c     now the acceleration calculated out of the potential of
c     grid 5 of galaxy gnr is computed, this acceleration is added
c     to every star of all galaxies  
c
            ix15 = nint(x15)      ! finding mesh point (integer)
            iy15 = nint(y15)
            iz15 = nint(z15)
            dx15 = x15 - dble(ix15) ! finding distance to the 
            dy15 = y15 - dble(iy15) ! middle of the cell
            dz15 = z15 - dble(iz15)
c
c     Differentiating the potential numerically to second order :
c     ===========================================================
c            
            fx  = dble( rho(ix15+1,iy15,iz15,5) )
     $           -dble( rho(ix15-1,iy15,iz15,5) )
            fy  = dble( rho(ix15,iy15+1,iz15,5) )
     $           -dble( rho(ix15,iy15-1,iz15,5) )
            fz  = dble( rho(ix15,iy15,iz15+1,5) )
     $           -dble( rho(ix15,iy15,iz15-1,5) )
            fxx = 2.d0*( dble( rho(ix15+1,iy15,iz15,5) )
     $           +dble( rho(ix15-1,iy15,iz15,5) )
     $           -2.d0* dble( rho(ix15,iy15,iz15,5) ) )
            fyy = 2.d0*( dble( rho(ix15,iy15+1,iz15,5) )
     $           +dble( rho(ix15,iy15-1,iz15,5) )
     $           -2.d0* dble( rho(ix15,iy15,iz15,5) ) )
            fzz = 2.d0*( dble( rho(ix15,iy15,iz15+1,5) )
     $           +dble( rho(ix15,iy15,iz15-1,5) )
     $           -2.d0* dble( rho(ix15,iy15,iz15,5) ) )
            fxy = 0.5d0*( dble( rho(ix15+1,iy15+1,iz15,5) )
     $           -dble( rho(ix15-1,iy15+1,iz15,5) )
     $           -dble( rho(ix15+1,iy15-1,iz15,5) )
     $           +dble( rho(ix15-1,iy15-1,iz15,5) ) )
            fxz = 0.5d0*( dble( rho(ix15+1,iy15,iz15+1,5) )
     $           -dble( rho(ix15-1,iy15,iz15+1,5) )
     $           -dble( rho(ix15+1,iy15,iz15-1,5) )
     $           +dble( rho(ix15-1,iy15,iz15-1,5) ) )
            fyz = 0.5d0*( dble( rho(ix15,iy15+1,iz15+1,5) )
     $           -dble( rho(ix15,iy15-1,iz15+1,5) )
     $           -dble( rho(ix15,iy15+1,iz15-1,5) )
     $           +dble( rho(ix15,iy15-1,iz15-1,5) ) )
c
c     Calculating the accelerations :
c     ===============================
c
            ex  = (fx + fxx*dx15 + fxy*dy15 + fxz*dz15)*sf5q
            ey  = (fy + fxy*dx15 + fyy*dy15 + fyz*dz15)*sf5q
            ez  = (fz + fxz*dx15 + fyz*dy15 + fzz*dz15)*sf5q
c
c     (2.) for all stars with r1 > rout1 : phi 14
c     ===========================================
c
c     now the acceleration out of grid 4 of galaxy gnr is
c     calculated and added to all stars of all galaxies that are
c     more distant to the focus of grids of galaxy gnr than rout1
c
            if (r1 .gt. rout1) then
c
               ix14  = nint(x14)
               iy14  = nint(y14)
               iz14  = nint(z14)
               dx14  = x14 - dble(ix14)
               dy14  = y14 - dble(iy14)
               dz14  = z14 - dble(iz14)
c
c     Adding the potential energy :
c     =============================
c
               hepot = hepot + dble( rho(ix14,iy14,iz14,4) )*sf4
     $              + dble( rho(ix15,iy15,iz15,5) )*sf5
c
               fx  = dble( rho(ix14+1,iy14,iz14,4) )
     $              -dble( rho(ix14-1,iy14,iz14,4) )
               fy  = dble( rho(ix14,iy14+1,iz14,4) )
     $              -dble( rho(ix14,iy14-1,iz14,4) )
               fz  = dble( rho(ix14,iy14,iz14+1,4) )
     $              -dble( rho(ix14,iy14,iz14-1,4) )
               fxx = 2.d0*( dble( rho(ix14+1,iy14,iz14,4) )
     $              +dble( rho(ix14-1,iy14,iz14,4) )
     $              -2.d0* dble( rho(ix14,iy14,iz14,4) ) )
               fyy = 2.d0*( dble( rho(ix14,iy14+1,iz14,4) )
     $              +dble( rho(ix14,iy14-1,iz14,4) )
     $              -2.d0* dble( rho(ix14,iy14,iz14,4) ) )
               fzz = 2.d0*( dble( rho(ix14,iy14,iz14+1,4) )
     $              +dble( rho(ix14,iy14,iz14-1,4) )
     $              -2.d0* dble( rho(ix14,iy14,iz14,4) ) )
               fxy = 0.5d0*( dble( rho(ix14+1,iy14+1,iz14,4) )
     $              -dble( rho(ix14-1,iy14+1,iz14,4) )
     $              -dble( rho(ix14+1,iy14-1,iz14,4) )
     $              +dble( rho(ix14-1,iy14-1,iz14,4) ) )
               fxz = 0.5d0*( dble( rho(ix14+1,iy14,iz14+1,4) )
     $              -dble( rho(ix14-1,iy14,iz14+1,4) )
     $              -dble( rho(ix14+1,iy14,iz14-1,4) )
     $              +dble( rho(ix14-1,iy14,iz14-1,4) ) )
               fyz = 0.5d0*( dble( rho(ix14,iy14+1,iz14+1,4) )
     $              -dble( rho(ix14,iy14-1,iz14+1,4) )
     $              -dble( rho(ix14,iy14+1,iz14-1,4) )
     $              +dble( rho(ix14,iy14-1,iz14-1,4) ) )
c
               ex = ex + (fx + fxx*dx14 + fxy*dy14 + fxz*dz14)*sf4q
               ey = ey + (fy + fxy*dx14 + fyy*dy14 + fyz*dz14)*sf4q
               ez = ez + (fz + fxz*dx14 + fyz*dy14 + fzz*dz14)*sf4q
c
               goto 19  ! no further accelerations for these stars
c
            endif
c
c     (3.) for all stars with r1 < rout1 : phi 11
c     ===========================================
c
c     if the star is inside of rout1 of galaxy gnr the acceleration
c     calculated out of grid 1 with higher resolution is added
c
            ix11 = nint(X11)
            iy11 = nint(y11)
            iz11 = nint(z11)
            dx11 = x11 - dble(ix11)
            dy11 = y11 - dble(iy11)
            dz11 = z11 - dble(iz11)
c
            fx  = dble( rho(ix11+1,iy11,iz11,1) )
     $           -dble( rho(ix11-1,iy11,iz11,1) )
            fy  = dble( rho(ix11,iy11+1,iz11,1) )
     $           -dble( rho(ix11,iy11-1,iz11,1) )
            fz  = dble( rho(ix11,iy11,iz11+1,1) )
     $           -dble( rho(ix11,iy11,iz11-1,1) )
            fxx = 2.d0*( dble( rho(ix11+1,iy11,iz11,1) )
     $           +dble( rho(ix11-1,iy11,iz11,1) )
     $           -2.d0* dble( rho(ix11,iy11,iz11,1) ) )
            fyy = 2.d0*( dble( rho(ix11,iy11+1,iz11,1) )
     $           +dble( rho(ix11,iy11-1,iz11,1) )
     $           -2.d0* dble( rho(ix11,iy11,iz11,1) ) )
            fzz = 2.d0*( dble( rho(ix11,iy11,iz11+1,1) )
     $           +dble( rho(ix11,iy11,iz11-1,1) )
     $           -2.d0* dble( rho(ix11,iy11,iz11,1) ) )
            Fxy = 0.5d0*( dble( rho(ix11+1,iy11+1,iz11,1) )
     $           -dble( rho(ix11-1,iy11+1,iz11,1) )
     $           -dble( rho(ix11+1,iy11-1,iz11,1) )
     $           +dble( rho(ix11-1,iy11-1,iz11,1) ) )
            fxz = 0.5d0*( dble( rho(ix11+1,iy11,iz11+1,1) )
     $           -dble( rho(ix11-1,iy11,iz11+1,1) )
     $           -dble( rho(ix11+1,iy11,iz11-1,1) )
     $           +dble( rho(ix11-1,iy11,iz11-1,1) ) )
            fyz = 0.5d0*( dble( rho(ix11,iy11+1,iz11+1,1) )
     $           -dble( rho(ix11,iy11-1,iz11+1,1) )
     $           -dble( rho(ix11,iy11+1,iz11-1,1) )
     $           +dble( rho(ix11,iy11-1,iz11-1,1) ) )
c
            ex = ex + (fx + fxx*dx11 + fxy*dy11 + fxz*dz11)*sf1q
            ey = ey + (fy + fxy*dx11 + fyy*dy11 + fyz*dz11)*sf1q
            ez = ez + (fz + fxz*dx11 + fyz*dy11 + fzz*dz11)*sf1q
c
c     (4.) for all stars with r1 > rcore1 : phi 12
c     ============================================
c
c     if the star is outside of the highest resolution center of
c     galaxy gnr the accelerations calculated out of grid 2 are
c     added now
c
            if (r1 .gt. rcore1) then
c
               ix12 = nint(x12)
               iy12 = nint(y12)
               iz12 = nint(z12)
               dx12 = x12 - dble(ix12)
               dy12 = y12 - dble(iy12)
               dz12 = z12 - dble(iz12)
c
c     Adding the potential energy :
c     =============================
c
               hepot = hepot + ( dble( rho(ix11,iy11,iz11,1) )*sf1
     $              +dble( rho(ix12,iy12,iz12,2) )*sf2
     $              +dble( rho(ix15,iy15,iz15,5) )*sf5 )
c
               fx=dble( rho(ix12+1,iy12,iz12,2) )
     $              -dble( rho(ix12-1,iy12,iz12,2) )
               fy=dble( rho(ix12,iy12+1,iz12,2) )
     $              -dble( rho(ix12,iy12-1,iz12,2) )
               fz=dble( rho(ix12,iy12,iz12+1,2) )
     $              -dble( rho(ix12,iy12,iz12-1,2) )
               fxx=2.d0*(dble( rho(ix12+1,iy12,iz12,2) )
     $              +dble( rho(ix12-1,iy12,iz12,2) )
     $              -2.d0*dble( rho(ix12,iy12,iz12,2) ) )
               fyy=2.d0*(dble( rho(ix12,iy12+1,iz12,2) )
     $              +dble( rho(ix12,iy12-1,iz12,2) )
     $              -2.d0*dble( rho(ix12,iy12,iz12,2) ) )
               fzz=2.d0*(dble( rho(ix12,iy12,iz12+1,2) )
     $              +dble( rho(ix12,iy12,iz12-1,2) )
     $              -2.d0*dble( rho(ix12,iy12,iz12,2) ) )
               fxy=0.5d0*(dble( rho(ix12+1,iy12+1,iz12,2) )
     $              -dble( rho(ix12-1,iy12+1,iz12,2) )
     $              -dble( rho(ix12+1,iy12-1,iz12,2) )
     $              +dble( rho(ix12-1,iy12-1,iz12,2) ) )
               fxz=0.5d0*(dble( rho(ix12+1,iy12,iz12+1,2) )
     $              -dble( rho(ix12-1,iy12,iz12+1,2) )
     $              -dble( rho(ix12+1,iy12,iz12-1,2) )
     $              +dble( rho(ix12-1,iy12,iz12-1,2) ) )
               fyz=0.5d0*(dble( rho(ix12,iy12+1,iz12+1,2) )
     $              -dble( rho(ix12,iy12-1,iz12+1,2) )
     $              -dble( rho(ix12,iy12+1,iz12-1,2) )
     $              +dble( rho(ix12,iy12-1,iz12-1,2) ) )
c
               ex = ex + (fx+fxx*dx12+fxy*dy12+fxz*dz12)*sf2q
               ey = ey + (fy+fxy*dx12+fyy*dy12+fyz*dz12)*sf2q
               ez = ez + (fz+fxz*dx12+fyz*dy12+fzz*dz12)*sf2q
c
               goto 19 ! no further accelerations for these stars
c
            endif
c
c     (5.) for all stars with r1 < rcore1 : phi 13
c     ============================================
c
c     all stars inside of rcore1 get accelerations added, 
c     calculated out of the highest resolution grid 3 of galaxy gnr
c
            ix13 = nint(x13)
            iy13 = nint(y13)
            iz13 = nint(z13)
            dx13 = x13 - dble(ix13)
            dy13 = y13 - dble(iy13)
            dz13 = z13 - dble(iz13)
c
c     Adding potential energy :
c     =========================
c
            hepot = hepot + ( dble( rho(ix11,iy11,iz11,1) )*sf1
     $           +dble( rho(ix13,iy13,iz13,3) )*sf3
     $           +dble( rho(ix15,iy15,iz15,5) )*sf5 )
c
            fx=dble( rho(ix13+1,iy13,iz13,3) )
     $           -dble( rho(ix13-1,iy13,iz13,3) )
            fy=dble( rho(ix13,iy13+1,iz13,3) )
     $           -dble( rho(ix13,iy13-1,iz13,3) )
            fz=dble( rho(ix13,iy13,iz13+1,3) )
     $           -dble( rho(ix13,iy13,iz13-1,3) )
            fxx=2.d0*( dble( rho(ix13+1,iy13,iz13,3) )
     $           +dble( rho(ix13-1,iy13,iz13,3) )
     $           -2.d0*dble( rho(ix13,iy13,iz13,3) ) )
            fyy=2.d0*( dble( rho(ix13,iy13+1,iz13,3) )
     $           +dble( rho(ix13,iy13-1,iz13,3) )
     $           -2.d0*dble( rho(ix13,iy13,iz13,3) ) )
            fzz=2.d0*( dble( rho(ix13,iy13,iz13+1,3) )
     $           +dble( rho(ix13,iy13,iz13-1,3) )
     $           -2.d0*dble( rho(ix13,iy13,iz13,3) ) )
            fxy=0.5d0*( dble( rho(ix13+1,iy13+1,iz13,3) )
     $           -dble( rho(ix13-1,iy13+1,iz13,3) )
     $           -dble( rho(ix13+1,iy13-1,iz13,3) )
     $           +dble( rho(ix13-1,iy13-1,iz13,3) ) )
            fxz=0.5d0*( dble( rho(ix13+1,iy13,iz13+1,3) )
     $           -dble( rho(ix13-1,iy13,iz13+1,3) )
     $           -dble( rho(ix13+1,iy13,iz13-1,3) )
     $           +dble( rho(ix13-1,iy13,iz13-1,3) ) )
            fyz=0.5d0*( dble( rho(ix13,iy13+1,iz13+1,3) )
     $           -dble( rho(ix13,iy13-1,iz13+1,3) )
     $           -dble( rho(ix13,iy13+1,iz13-1,3) )
     $           +dble( rho(ix13,iy13-1,iz13-1,3) ) )
c
            ex = ex + (fx+fxx*dx13+fxy*dy13+fxz*dz13)*sf3q
            ey = ey + (fy+fxy*dx13+fyy*dy13+fyz*dz13)*sf3q
            ez = ez + (fz+fxz*dx13+fyz*dy13+fzz*dz13)*sf3q
c
c     End of calculation of accelerations
c
 19         continue
c
c     Calculating the old internal kin. energy and angular mom. :
c     ===========================================================
c
            if (gnr .eq. 1) then
               do j = 1,6
                  rv(j) = star(j,i) - gcms(j,k)
                  tv(j) = star(j,i) - tcms(j)
               enddo
               gekin(k) = gekin(k) + rv(4)**2 + rv(5)**2 + rv(6)**2
               tekin    = tekin    + tv(4)**2 + tv(5)**2 + tv(6)**2
               lgx(k)   = lgx(k)   + rv(2)*rv(6) - rv(3)*rv(5)
               lgy(k)   = lgy(k)   + rv(3)*rv(4) - rv(1)*rv(6)
               lgz(k)   = lgz(k)   + rv(1)*rv(5) - rv(2)*rv(4)
               tlx      = tlx      + tv(2)*tv(6) - tv(3)*tv(5)
               tly      = tly      + tv(3)*tv(4) - tv(1)*tv(6)
               tlz      = tlz      + tv(1)*tv(5) - tv(2)*tv(4)
            endif
c            
c     Convert velocities to double precission :
c     =========================================
c
            dbvx = dble(star(4,i))
            dbvy = dble(star(5,i))
            dbvz = dble(star(6,i))
c
c     Leap-Frog Integration :
c     =======================
c

c++
            if( i .le. 5 ) then 
c!!            write( 6,* ) i, r1, ex, ey, ez 
            endif 
c++
            dbvx = dbvx + dble(ex) * dble(dt(gnr))
            dbvy = dbvy + dble(ey) * dble(dt(gnr))
            dbvz = dbvz + dble(ez) * dble(dt(gnr))

c!!   alla vamos, anyado df a la fuerza heacha por halo, gnr=1
            if(gnr.eq.1) then
               dbvx = dbvx - adfx* dble(dt(gnr))
               dbvy = dbvy - adfy* dble(dt(gnr))
               dbvz = dbvz - adfz* dble(dt(gnr))
            endif
c!!   

c
c     Convert back to single precission:
c     ==================================
c
            star(4,i) = real(dbvx)
            star(5,i) = real(dbvy)
            star(6,i) = real(dbvz)
c            
c     Update center of velocity (galaxy) :
c     ====================================
c
            if (gnr .eq. gnum) then
               do j = 1,3
                  hcms(j,k) = hcms(j,k) + star(j+3,i)
               enddo
            endif
c
 20         continue
c            
         enddo                  ! end of particle loop
c
         if (gnr .eq. gnum) then
            do j = 1,3
               hcms(j,k) = hcms(j,k) / grest(k)
            enddo
         endif
c
c     Calculating potential energy per galaxy and of system :
c     =======================================================
c
         if (k .eq. gnr) gepot(gnr) = - mass(gnr) * hepot
         epot = epot - mass(k) * hepot
c
c     Add up kinetic energy and angular momentum (system) :
c     =====================================================
c
         if (gnr .eq. 1) then
            ekin = ekin + mass(k) * tekin
            lx   = lx   + mass(k) * tlx
            ly   = ly   + mass(k) * tly
            lz   = lz   + mass(k) * tlz
         endif               
c     
      enddo                     ! end of galaxy loop
c
c     Update center of velocity (system) :
c     ====================================
c
      if (gnr .eq. gnum) then
c
         do j = 4,6
            tcms(j) = 0.0
         enddo
         tmass = 0.0 
         do k = 1,gnum
            tmass = tmass + mass(k) * grest(k)
            do j = 4,6
               tcms(j) = tcms(j) + hcms(j-3,k) * mass(k) * grest(k)
            enddo
         enddo
         do j = 4,6
            tcms(j) = tcms(j) / tmass
         enddo
c
c     Calculating the new kin. energy and angular mom. :
c     ==================================================
c
         do k = 1,gnum          ! starting galaxy loop
c
            tekin = 0.0
            tlx   = 0.0
            tly   = 0.0
            tlz   = 0.0
c
            if (k .eq. 1) then  ! finding the number of the first
               start = 1        ! and the last star of galaxy k
               istop  = gstno(1)
            else
               start = 0
               do i = 1,k-1
                  start = start + gstno(i)
               enddo
               istop  = start + gstno(k)
               start  = start + 1
            endif
c
            do i = start,istop  ! start particle loop
               if (star(1,i) .gt. 1.e+10) goto 30
               do j = 1,3
                  rv(j)   = star(j,i) - gcms(j,k)
                  rv(j+3) = star(j+3,i) - hcms(j,k)
                  tv(j)   = star(j,i) - tcms(j)
                  tv(j+3) = star(j+3,i) - tcms(j+3)
               enddo
               gekin(k) = gekin(k) + rv(4)**2 + rv(5)**2 + rv(6)**2
               tekin    = tekin    + tv(4)**2 + tv(5)**2 + tv(6)**2
               lgx(k)   = lgx(k)   + rv(2)*rv(6) - rv(3)*rv(5)
               lgy(k)   = lgy(k)   + rv(3)*rv(4) - rv(1)*rv(6)
               lgz(k)   = lgz(k)   + rv(1)*rv(5) - rv(2)*rv(4)
               tlx      = tlx      + tv(2)*tv(6) - tv(3)*tv(5)
               tly      = tly      + tv(3)*tv(4) - tv(1)*tv(6)
               tlz      = tlz      + tv(1)*tv(5) - tv(2)*tv(4)
 30            continue
            enddo               ! end particle loop
c
c     Now the exact kin. energy and angular momentum can be derived :
c     ===============================================================
c
            gekin(k) = gekin(k) * 0.25 * mass(k)
            lgx(k)   = lgx(k) * 0.5 * mass(k)
            lgy(k)   = lgy(k) * 0.5 * mass(k)
            lgz(k)   = lgz(k) * 0.5 * mass(k)
            ekin     = ekin + mass(k) * tekin
            lx       = lx   + mass(k) * tlx
            ly       = ly   + mass(k) * tly
            lz       = lz   + mass(k) * tlz
c
         enddo                  ! end galaxy loop
c
         ekin = 0.25 * ekin
         lx   = 0.5 * lx
         ly   = 0.5 * ly
         lz   = 0.5 * lz
c
      endif
c
      return
c
      end



      subroutine binney(vc,sigt,sigp,e,Mgal,Ms,rcut,gamma,
     &     x,y,z,u,v,w,ax,ay,az)
      double precision x,y,z,u,v,w,ax,ay,az,vc,sigt,sigp
      parameter(pi=3.14159)
      parameter(G=2.)
      double precision vst,vsp,e,dintt,dintp,qf,ft,fp,q,Bt,Bp
      real Mgal,Ms,rcut,gama
      integer nstep
      parameter(nstep=100000)
      parameter (qf=1000.)
c---------------------------
c     .x es con theta,phi usando las vel.
c     si anyado un numero es con theta,phi usando las posic
c     nota: .x.1 es con halo fijo
      logl=3.75
c      logl=4.75 .x.2
c      logl=5.75     .x.3
c      logl=37.5     .x.4
c-----------------------------
      
      r=sqrt(x**2.+y**2.+z**2.)

c     densidad halo isotermo
      rhoh=Mgal/(2.*pi**(3./2.)*rcut)*exp(-r**2./rcut**2.)/
     &     (r**2.+gamma**2.)

c     constante
      cte=2.*sqrt(2.*pi)*rhoh*G**2.*logl*Ms*sqrt(1.-e**2.)

c     integracion

c           1_calculo vst,vsp
      vmodul=sqrt(u**2.+v**2.+w**2.)
      theta=acos(z/r)
      vst=vmodul*sin(theta)
      vsp=vmodul*cos(theta)
      vst2=sqrt(u**2.+v**2.)
      vsp2=w
      PRINT*, 'VELOCIDADES',vst,vst2,vsp,vsp2,theta*180./3.14159
      
c           2_integra
      
      
      dq=qf/real(nstep)
      
      dintt=dq/2.*( exp(-vst**2./(2.*sigt**2.)-
     &     vsp**2./(2.*sigt**2.)/(1.-e**2.))/sqrt(1.-e**2.) + 
     &     exp(-vst**2./(2.*sigt**2.)/(1+qf)-
     &     vsp**2./(2.*sigt**2.)/(1.-e**2.+qf))/
     &     ((1.+qf)**2.*sqrt(1.-e**2.+qf))  )
      
      dintp=dq/2.*( exp(-vst**2./(2.*sigt**2.)-
     &     vsp**2./(2.*sigt**2.)/(1.-e**2.))/sqrt(1.-e**2.) + 
     &     exp(-vst**2./(2.*sigt**2.)/(1+qf)-
     &     vsp**2./(2.*sigt**2.)/(1.-e**2.+qf))/
     &     ((1.+qf)**2.*sqrt(1.-e**2.+qf)**3.)  )
      
      q=dq
      do i=2,nstep-1
         q=q+dq
         fft=exp(-vst**2./(2.*sigt**2.)/(1+q)-
     &     vsp**2./(2.*sigt**2.)/(1.-e**2.+q))/
     &     ((1.+q)**2.*sqrt(1.-e**2.+q))
         ffp=exp(-vst**2./(2.*sigt**2.)/(1+q)-
     &     vsp**2./(2.*sigt**2.)/(1.-e**2.+q))/
     &     ((1.+q)**2.*sqrt(1.-e**2.+q)**3.)
         dintt=dintt+fft*dq
         dintp=dintp+ffp*dq
      enddo

      Bt=dintt
      Bp=dintp

      at=cte*Bt*vst/(sigt**2.*sigp)
      ap=cte*Bp*vsp/(sigt**2.*sigp)

      

      phi=atan2(v,u)
     
      ax=at*cos(phi)
      ay=at*sin(phi)
      az=ap
      
      PRINT*, 'Bt',Bt,'Bp',Bp,'at',at,'ap',ap,'axyz',ax,ay,az,
     &     'phi',phi,atan(v/u) 
      
      return
      end

      FUNCTION erfcc(x)
      REAL erfcc,x
      REAL t,z
      z=abs(x)
      t=1./(1.+0.5*z)
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*
     *(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*
     *(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.-erfcc
      return
      END

c! anyadiendo pot externo
c**********************************************************************
c
c     P U S H E R 1 . F
c
c**********************************************************************
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
      double precision xcms,ycms,zcms,rad
      double precision ex,ey,ez,fxx,fxy,fxz,fyy,fyz,fzz
      double precision x,y,z,dx,dy,dz,mfrac,mfracd
      double precision x11,x12,x13,x14,x15
      double precision y11,y12,y13,y14,y15
      double precision z11,z12,z13,z14,z15
      double precision r1,fnh,dx15,dy15,dz15,fx,fy,fz,ffx,ffy,ffz,dr
      double precision dx14,dy14,dz14
      double precision dx11,dy11,dz11
      double precision dx12,dy12,dz12
      double precision dx13,dy13,dz13
      double precision dbx,dby,dbz
      double precision dbvx,dbvy,dbvz
      real*8 tiem,tcurrent,thub,tf
      integer*4          start,istop,gnr,step
      integer*4          i,j,k,l,ibound(istop)
      integer*4          ix11,iy11,iz11
      integer*4          ix12,iy12,iz12
      integer*4          ix13,iy13,iz13
      integer*4          ix14,iy14,iz14
      integer*4          ix15,iy15,iz15
      real*4             rv(6),tv(6),tekin,tlx,tly,tlz
      integer*4 knum,kk,lf,kbound,knumd,nd,kunb
      real*4 detot
      real*8 pot,Ms0,Msf,Ms,desatk,desatp

c     particulas del disco!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nd=0

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
c
c     Starting the loop over all galaxies :
c     =====================================
c
c     contador de friccion a 0
      lf=0
c     masa inicial/final!!!!!!!!!!!!!!!
      Ms0=dble(fh(5,gnr))
      Msf=0.5* Ms0
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

c!       
c     numero de part ligadas
         knum =0
         knumd=0
         kunb=0
         do kk=start,istop            
            if (k.eq.gnr.and.de(1,kk)+de(2,kk).lt.0..and.
     &           tm(kk).lt.fh(4,gnr)) then             
               knum=knum+1                              
            endif
            if (k.eq.gnr.and.de(1,kk)+de(2,kk).gt.0..and.
     &           tm(kk).lt.fh(4,gnr)) then  
               tcurrent=  dble(ih(3))*dble(fh(4,gnr))*sqrt(2.)
               tm(kk)=tcurrent
               kunb=kunb+1
            endif
            if (k.eq.gnr) then  !energias a 0
               de(1,kk)=0.
               de(2,kk)=0.
            endif
            if (step.eq.1) tm(kk)=0.0
         enddo
         mfrac =dble(knum )/dble(istop-start)

         
c
         hepot  = 0.d0        ! initialize potential energy
         tekin  = 0.0
         tlx    = 0.0
         tly    = 0.0
         tlz    = 0.0

c     3 iteraciones para cualcula Phi_bound
         l=0
c
c     Starting the particle loop of galaxy k :
c     ========================================
c
         do  i = start,istop  ! starting particle loop
c
c     Skip stars outside the simulation :
c     ===================================
c
            if (star(1,i) .gt. 1.e+10) then
!               print*, 'estrella a tomar por saco ',i,star(1,i)
               goto 20 
            endif
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
               if (tm(I).lt.fh(4,gnr)) then
               de(2,i)=de(2,i) - dble( rho(ix14,iy14,iz14,4) )*sf4
     $              - dble( rho(ix15,iy15,iz15,5) )*sf5
               endif

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

               if (tm(I).lt.fh(4,gnr)) then
               de(2,i)=de(2,i)-( dble( rho(ix11,iy11,iz11,1) )*sf1
     $              +dble( rho(ix12,iy12,iz12,2) )*sf2
     $              +dble( rho(ix15,iy15,iz15,5) )*sf5 )
               endif
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
            if (tm(I).lt.fh(4,gnr)) then
            de(2,i)=de(2,i)-( dble( rho(ix11,iy11,iz11,1) )*sf1
     $           +dble( rho(ix13,iy13,iz13,3) )*sf3
     $           +dble( rho(ix15,iy15,iz15,5) )*sf5 )
            endif
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
            dbx  = dble(star(1,i))
            dby  = dble(star(2,i))
            dbz  = dble(star(3,i))
            dbvx = dble(star(4,i))
            dbvy = dble(star(5,i))
            dbvz = dble(star(6,i))
c
c     Leap-Frog Integration :
c     =======================
c

c! potencial externo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (k.eq.gnr) then

            tcurrent=  dble(ih(3))*dble(fh(4,gnr))*sqrt(2.)
            tf      =  dble(ih(4))*dble(fh(4,gnr))*sqrt(2.)
            thub    =  14. / 0.013
            tiem    =  thub - tf + tcurrent          
               
            call force 
     &           (tiem,tf,dbx,dby,dbz,dbvx,dbvy,dbvz,
     &           dble(fh(14,gnr)),dble(fh(15,gnr)),dble(fh(16,gnr)),
     &           dble(fh(17,gnr)),dble(fh(18,gnr)),dble(fh(19,gnr)),
     &           fx,fy,fz,pot)
!               fx=0.d0
!               fy=0.d0
!               fz=0.d0

            if (lf.eq.0) then
c     masa
               Ms=Ms0*mfrac
!               call masa(tiem,tf,thub,Ms0,Msf,Ms)
c     friccion en xcm,ycm,zcm
!               call friccion(tiem,tf,Ms,
!     &              dble(fh(14,gnr)),dble(fh(15,gnr)),dble(fh(16,gnr)),
!     &              dble(fh(17,gnr)),dble(fh(18,gnr)),dble(fh(19,gnr)),
!     &              ffx,ffy,ffz)
               ffx=0.d0
               ffy=0.d0
               ffz=0.d0
               print*,step,tiem,mfrac,kunb,
     &              sqrt(fh(14,gnr)**2 + fh(15,gnr)**2 + fh(16,gnr)**2)
               lf=1
            endif

            else !  k=!gnr
               fx=0.
               fy=0.
               fz=0.
               ffx=0.
               ffy=0.
               ffz=0.
            endif !k=gnr?

c     energia cinetica !G=2 !sat frame
            desatk=0.25*real( (dbvx-dble(fh(17,gnr)))**2.+
     &      (dbvy-dble(fh(18,gnr)))**2.+(dbvz-dble(fh(19,gnr)))**2.)
c     energia pot           !sat frame                        
            desatp=de(2,i)
            detot=desatp+desatk           
            de(1,i)=desatk
            de(2,i)=de(2,i)
c     


c     energia cinetica !G=2 !galaxy frame
            de(3,i)=0.25*real(dbvx**2.+dbvy**2.+dbvz**2.)
c     energia potencial     !galaxy frame
            de(4,i)=real(pot)


            if (detot.lt.0.) then
               dbvx = dbvx + (dble(ex)+fx +ffx) * dble(dt(gnr))
               dbvy = dbvy + (dble(ey)+fy +ffy) * dble(dt(gnr))
               dbvz = dbvz + (dble(ez)+fz +ffz) * dble(dt(gnr))
            else
               dbvx = dbvx + (dble(ex)+fx     ) * dble(dt(gnr))
               dbvy = dbvy + (dble(ey)+fy     ) * dble(dt(gnr))
               dbvz = dbvz + (dble(ez)+fz     ) * dble(dt(gnr))
            endif
        
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

      subroutine param(t,tf,Md,a,b,Mb,c,Mh,rch,qx,qy,qz)
      implicit real*8 (a-h,o-z)
      real*8 mf,Ms0,Md,Mb,Mh,rch(2),Mh0
      real*8 Munit
      parameter (Munit=5.6e10,runit=3.5,G=2.)
c     parametros de la galaxia

c     disco
      Md=7.5e10 /Munit*0.d0                 !Disc mass
      a =3.5 /runit   
      b =0.3 /runit

c     bulge
      Mb=1.3e10 /Munit*0.d0                 !bulge mass
      c =1.2    /runit

c     halo
      thub=14./ 0.013       
!      call nfwprop(t,thub,Mh0,rv0,ac,Mh,rch)
      Mh=1.e12   /Munit
      rch(2)=258./runit
      rch(1)=21.5/runit
      qx=1.d0
      qy=1.d0
      qz=1.d0
      return
      end

      subroutine friccion(t,tf,Ms,xcm,ycm,zcm,ucm,vcm,wcm,fx,fy,fz)
      implicit real*8 (a-h,o-z)
      real*8 Mh,rch(2),Md,Mb,Munit,Mh0,Ms0,Msf,Ms,mf
      real*8 G,fdf(3),fdfd(3)
      parameter (Munit=5.6e10,runit=3.5,G=2.)

      do i=1,3
         fdf(i)=0.
         fdfd(i)=0.
      enddo

      call param(t,tf,Md,a,b,Mb,c,Mh,rch,qx,qy,qz)

      call dfh(Ms,Mh,rch,qx,qy,qz,xcm,ycm,zcm,ucm,vcm,wcm,fdf)
      call dfdb(Ms,Md,a,b,Mb,c,xcm,ycm,zcm,ucm,vcm,wcm,fdfd)
      fx= (fdf(1) + fdfd(1) )*G
      fy= (fdf(2) + fdfd(2) )*G
      fz= (fdf(3) + fdfd(3) )*G
      return
      end


      subroutine force 
     &     (t,tf,x,y,z,u,v,w,xcm,ycm,zcm,ucm,vcm,wcm,fx,fy,fz,pot)
      implicit real*8 (a-h,o-z)
      real*8 Mh,rch(2),Md,Mb,Munit,Mh0,Ms0,Msf,Ms
      real*8 fd(3),fb(3),fh(3),G,fdf(3)
      parameter (Munit=5.6e10,runit=3.5,G=2.)

      call param(t,tf,Md,a,b,Mb,c,Mh,rch,qx,qy,qz)

c     potencial
      call forced( Md,a,b,x,y,z,fd,potd)
      call forceb( Mb,c,  x,y,z,fb,potb)
!      call forcehnfw(Mh,rch,qx,qy,qz,x,y,z,fh,poth)
      call force_an  (Mh,rch,qz,x,y,z,fh,poth)
!      if (abs(qz-1.).lt.0.05) then
!         print*,'Error: usamos force_an(q=1) con q=',q
!         print*,'Usar forcehnfw (cuidado con la relacion v(0)--ecc.)'
!         print*,'cambiar centreorbit.f'
!         stop
!      endif

c     suma de fuerzas

      fx=(fd(1)+fb(1)+fh(1))*G
      fy=(fd(2)+fb(2)+fh(2))*G
      fz=(fd(3)+fb(3)+fh(3))*G
c     potencial
      pot=potd+potb+poth
      return
      end

      subroutine masa(t,tf,thub,Ms0,Msf,Ms)
      implicit real*8 (a-h,o-z)
      real*8 Ms0,Msf,Ms,n
      t2=t*0.013
      thub2=thub*0.013
      tf2=tf*0.013
      Ms=Ms0*(1.- (exp(-(thub2-t2))**0.3-exp(-tf2)**0.3) *(1.-Msf/Ms0) )
      return
      end


      subroutine nfwprop(t,thub,Mh0,rv0,ac,Mh,rch)
      implicit real*8 (a-h,o-z)
      real*8 Mh0,Mh,rch(2)
      parameter(S=2.)
      h=0.7
      pi=4.*atan(1.)
c     densidades cosmol.
      om0=0.3
      ol0=0.7
c     conversion t <--> z
      a=0.9*t/thub+0.1
      z=1./a-1.
c     basado en Bullock et  al 2001 MNRAS 321,559
      e2=om0*(1.+z)**3.+ol0
      om=om0*(1.+z)**3./e2         
      del=(18.*pi**2. + 82.*(om-1.)-39.*(om-1.)**2.)/om
c     masa !!Wechsler et al 2002 
      Mh=Mh0*exp(-ac*S*(1./a-1.))
c     radio vir !!Bullock et  al 2001
      rch(2)=75./h/(1.+z)/3.5*
     &     (Mh*5.6e10/(1.e11/h) *200./om0/del)**(1./3.)
c     concentracion !!Wechsler et al 2002 
      cv=4.1 * a/ac
c     radio de escala
      rch(1)= rch(2)/cv
      return
      end

      subroutine forcehnfw(Mh,rc,qx,qy,qz,x,y,z,f,pot)
c     fuerzas NFW
      implicit real*8 (a-h,o-z)
      parameter (nsmax=5000)
      real*8 f(3),Mh,rc(2)
      parameter (np=20,np2=40)
      real*8 xx(np),ww(np),xx2(np2),ww2(np2)
      data xx/0.257677531,1.35105273,3.29120894,6.03311356,9.51255351,
     &13.647987,18.3424874,23.4860217,28.9580306,34.6302555,40.3697445,
     &46.0419694,51.5139783,56.6575126,61.352013,65.4874465,68.9668864,
     &71.7087911,73.6489473,74.7423225/
      data ww/0.660525268,1.52255362,2.35020181,3.12287781,3.82237949,
     &4.43229495,4.93832394,5.3286041,5.59398699,5.72825202,5.72825202,
     &5.59398699,5.3286041,4.93832394,4.43229495,3.82237949,3.12287781,
     &2.35020181,1.52255362,0.660525268/      
      data xx2/0.330429429,1.73883024,4.26375938,7.8905964,12.5975984,
     &18.3564737,25.1325931,32.8852067,41.5676904,51.1278272,61.5081216,
     &72.6461457,84.474914,96.9232872,109.916399,123.376108,137.221465,
     &151.369204,165.734237,180.230172,194.769828,209.265763,223.630796,
     &237.778535,251.623892,265.083601,278.076713,290.525086,302.353854,
     &313.491878,323.872173,333.43231,342.114793,349.867407,356.643526,
     &362.402402,367.109404,370.736241,373.26117,374.669571/
      data ww2/0.847739456,1.96842835,3.07894845,4.17109672,5.23818881,
     &6.27378662,7.2716565,8.22579528,9.13046393,9.98022131,10.7699567,
     &11.4949205,12.1507525,12.7335086,13.2396839,13.6662342,14.0105942,
     &14.2706929,14.4449659,14.5323652,14.5323652,14.4449659,14.2706929,
     &14.0105942,13.6662342,13.2396839,12.7335086,12.1507525,11.4949205,
     &10.7699567,9.98022131,9.13046393,8.22579528,7.2716565,6.27378662,
     &5.23818881,4.17109672,3.07894845,1.96842835,0.847739456/
      external rhoh
      pi=4.*atan(1.)
      

c     distancia
      Rpl=sqrt(x**2.+y**2.)
c     halo achatamiento
      dix=0.
      diy=0.
      diz=0.
      do i=1,np2
         am=sqrt( x**2./(qx**2.+xx2(i)) + y**2./(qy**2.+xx2(i)) + 
     &        z**2./(qz**2.+xx2(i)) )
         drhoh=rhoh(am,rc(1),rc(2),Mh)
         dix=dix + ww2(i)* drhoh/
     &   (qx**2.+xx2(i))**1.5/(qy**2.+xx2(i))**0.5/(qz**2.+xx2(i))**0.5
         diy=diy + ww2(i)* drhoh/
     &   (qx**2.+xx2(i))**0.5/(qy**2.+xx2(i))**1.5/(qz**2.+xx2(i))**0.5
         diz=diz + ww2(i)* drhoh/
     &   (qx**2.+xx2(i))**0.5/(qy**2.+xx2(i))**0.5/(qz**2.+xx2(i))**1.5
      enddo
      f(1)=-2.*pi*qx*qy*qz*x*dix 
      f(2)=-2.*pi*qx*qy*qz*y*diy 
      f(3)=-2.*pi*qx*qy*qz*z*diz 
      
c     potencial      
      rho0=Mh/(4.*pi*rc(1)**3.)/
     &     (log(1.+rc(2)/rc(1))-rc(2)/(rc(2)+rc(1)))
      pot=0.
      do i=1,np2
         am=sqrt( x**2./(qx**2.+xx2(i)) + y**2./(qy**2.+xx2(i)) + 
     &            z**2./(qz**2.+xx2(i)) )/rc(1)
         fun=am/(1.+am)
     &      /sqrt(qx**2.+xx2(i))/sqrt(qy**2.+xx2(i))/sqrt(qz**2.+xx2(i))
         pot=pot+ww2(i)*fun
      enddo
      pot=pot*2.*pi*qx*qy*qz*rc(1)**2.*rho0-4.*pi*rc(1)**2.*rho0
      return
      end

 
      subroutine force_an(Mh,rc,qh,x,y,z,f,pot)
      implicit real*8 (a-h,o-z)
      real*8 f(3),Mh,rc(2)
      rad=sqrt(x**2.+y**2.+z**2.)
      phi0=-Mh/(log(1.+rc(2)/rc(1))-rc(2)/(rc(2)+rc(1)))
      deriv=( rad/rc(1)/(1.+rad/rc(1)) - log(1.+rad/rc(1))  )/rad**2.
      f(1)=-phi0*deriv*x/rad
      f(2)=-phi0*deriv*y/rad
      f(3)=-phi0*deriv*z/rad
      pot=phi0*log(1.+rad/rc(1))/rad
      return
      end

      function rhoh(r,r0,rm,Mh)
      real*8 rhoh,r,r0,rm,Mh,pi
      pi=4.*atan(1.)
      rhoh=Mh/4./pi/r/(r+r0)**2./(log(1.+rm/r0)-rm/(rm+r0))
      return
      end

      
      subroutine forceb(Mb,ab,x,y,z,f,pot) 
      implicit real*8 (a-h,o-z)
      real*8 Mb,fx,fy,fz,f(3)
      rr=sqrt(x**2.+y**2.+z**2.)
      f(1)=-Mb/(rr+ab)**2. *x/rr
      f(2)=-Mb/(rr+ab)**2. *y/rr
      f(3)=-Mb/(rr+ab)**2. *z/rr
      pot=-Mb/(rr+ab)
      return
      end

      subroutine forced(Md,a,b,x,y,z,f,pot)
      implicit real*8 (a-h,o-z)
      real*8 Md,fx,fy,fz,f(3)
      Rpl=sqrt(x**2.+y**2.)

      Rm1=(Rpl**2.+a**2.*(1.+b/a*sqrt(1.+(z/b)**2.))**2.)**1.5
      f(1)=-Md/Rm1 * x 
      f(2)=-Md/Rm1 * y 
      f(3)=-Md/Rm1/(b/a)*(1.+b/a*sqrt(1.+(z/b)**2.))
     &     /sqrt(1.+(z/b)**2.) * z 
      pot=-Md/sqrt( Rpl**2.+ (a + sqrt(z**2.+b**2.) )**2. )
      return
      end

      subroutine dfh(Ms,Mh,rc,qx,qy,qz,x,y,z,u,v,w,f)
      implicit real*8 (a-h,o-z)
      real*8 f(3),Mh,Ms,rc(2),lnl
      parameter (lnl=2.1)
      parameter (np=20)
      real*8 xx(np),ww(np)
      data xx/0.257677531,1.35105273,3.29120894,6.03311356,9.51255351,
     &13.647987,18.3424874,23.4860217,28.9580306,34.6302555,40.3697445,
     &46.0419694,51.5139783,56.6575126,61.352013,65.4874465,68.9668864,
     &71.7087911,73.6489473,74.7423225/
      data ww/0.660525268,1.52255362,2.35020181,3.12287781,3.82237949,
     &4.43229495,4.93832394,5.3286041,5.59398699,5.72825202,5.72825202,
     &5.59398699,5.3286041,4.93832394,4.43229495,3.82237949,3.12287781,
     &2.35020181,1.52255362,0.660525268/
      external rhoh
      pi=4.*atan(1.)      

c     centro?
      rad=sqrt(x**2.+y**2.+z**2.)
      if (rad.le.rc(1)/100.)return
         
    
c     distancia
      am=sqrt( x**2./qx**2. + y**2./qy**2. + z**2./qz**2. )

      cte=2.*sqrt(2.*pi)*  rhoh(am,rc(1),rc(2),Mh)*lnl*Ms
      
c     velocidad de dispersion
      call vdisph(am,Mh,rc,s) !s_i=sigma_i^2

c     velocidades
      vpl=sqrt(u**2.+v**2.)/sqrt(2.)
      vz =w                /sqrt(2.)          
      

c     integracion en velocidades      
      dix=0.
      diy=0.
      diz=0.
      do i=1,np
         am=sqrt( u**2./(qx**2.+xx(i)) + v**2./(qy**2.+xx(i)) + 
     &            w**2./(qz**2.+xx(i)) )

         drhov=exp(-am**2./2./s**2.)
         dix=dix + ww(i)* drhov/
     &   (qx**2.+xx(i))**1.5/(qy**2.+xx(i))**0.5/(qz**2.+xx(i))**0.5
         diy=diy + ww(i)* drhov/
     &   (qx**2.+xx(i))**0.5/(qy**2.+xx(i))**1.5/(qz**2.+xx(i))**0.5
         diz=diz + ww(i)* drhov/
     &   (qx**2.+xx(i))**0.5/(qy**2.+xx(i))**0.5/(qz**2.+xx(i))**1.5
      enddo
      f(1)=-cte*dix/(s)**3.*u
      f(2)=-cte*diy/(s)**3.*v
      f(3)=-cte*diz/(s)**3.*w

      return
      end
 
      subroutine vdisph(r,Mh,rc,sr)
      implicit real*8 (a-h,o-z)
      real*8 Mh,rc(2),masah
      parameter (np=20)
      real*8 xx(np),ww(np)
      external rhoh,masah
      pi=4.*atan(1.)
      b=10.*r
      call gauleg(r,b,xx,ww,np)
      s=0.
      do k=1,np
         drhoh=  rhoh(xx(k),rc(1),rc(2),Mh)
         dm=masah(xx(k),Mh,rc)
         s=s + ww(k)*drhoh*dm/xx(k)**2. !1D sigma
      enddo
      drhoh=  rhoh(r,rc(1),rc(2),Mh)
      sr=s/drhoh

      return
      end

      function masah (r,Mh,rc)
      implicit real*8(a-h,o-z)
      parameter (np=15)
      real*8 Mh,xx(np),ww(np),masah,rc(2)
      external rhoh
      pi=4.*atan(1.)
c     calcula M(r) 
      a=0.
      call gauleg(a,r,xx,ww,np)
      dm=0.
      do k=1,np
         drhoh=  rhoh(xx(k),rc(1),rc(2),Mh)
         dm=dm+ ww(k)*4.*pi*xx(k)**2.*drhoh
      enddo
      masah=dm
      return
      end



      subroutine dfdb(Ms,Md,a,b,Mb,c,x,y,z,u,v,w,f)
      implicit real*8 (a-h,o-z)
      real*8 r(6),f(3),Md,Mb,Ms,lnld,lnlb
      parameter (lnld=0.5)
      external erf,rhod
      pi=4.*atan(1.)

c     distancia
      Rpl=sqrt(x**2.+y**2.)
      z  = r(3)
c     from Lewis&Freeman 1989, sigr(0)=100 km/s
      sigr=(100./262.)*exp(-Rpl/a/2.)
      sig=sigr !1D 
!      sig=sqrt(3.)*sigr !3D 
c     velocidad circular  eq.2-170 BT
      if (Rpl.gt.0.01) then
         vc=220./262.!sqrt(  Md/Rpl*( 1.-exp(-Rpl/a)*(1.+ Rpl/a) )  )
      else
         vc=0.
      endif
c     vel relativa
      phi=atan2(r(2),r(1))
      vrelx=u/sqrt(2.)-vc*sin(phi)
      vrely=v/sqrt(2.)+vc*cos(phi)
      vrelz=w/sqrt(2.)

      vrel =sqrt(vrelx**2.+vrely**2.+vrelz**2.)
c     calculo de la cte
      xvel=vrel/(sqrt(2.)*sig)   
      Xd=(erf(xvel)-2.*xvel/sqrt(pi)*exp(-xvel**2.))
      cte=4.*pi*Ms*rhod(Rpl,z,Md,a,b)*Xd*lnld/vrel**3.
c     fuerza
      f(1)=-cte*vrelx
      f(2)=-cte*vrely
      f(3)=-cte*vrelz
      return
      end

      function rhod(R,z,Md,a,b)
      implicit real*8 (a-h,o-z)
      real*8  Md
      pi=4.*atan(1.)
!      rhod=(Md*b**2./4./pi)*(a*R**2. + (a + 3.*sqrt(z**2.+b**2.))*
!     &     (a + sqrt(z**2.+b**2.) )**2.)/
!     &     (R**2. + (a + sqrt(z**2.+b**2.) )**2.)**2.5 / 
!     &     (z**2.+b**2.)**1.5
      rhod=a*Md/2./pi/(R**2.+a**2.)**(3./2.)*exp(-abs(z)/b)
      return
      end       
      

      function erf(x)
      double precision erf,x,pi,ss,zero
      integer np
      parameter(np=10)
      double precision xx(np),ww(np)
      external expon
      pi=4.*atan(1.)
      zero=0.
      ss=0.
      call gauleg(zero,x,xx,ww,np)
      do i=1,np
         ss=ss+ww(i)*expon(xx(i))
      enddo
      erf=2./sqrt(pi)*ss
      return
      end

      function expon(x)
      double precision x,expon
      expon=exp(-x**2.)
      return
      end

c     numerical recipies
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      real*8 x1,x2,x(n),w(n)
      real*8 EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      real*8 p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END    

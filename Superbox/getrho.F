c**********************************************************************
c
c     G E T R H O . F
c
c**********************************************************************
c
c     This subroutine determines the mass density of every mesh
c     point of all grids for galaxy gnr - last modified 15.05.2k (cmb).
c     This subtoutine decides the membership of stars to grids, and
c     whether to remove stars from the calculation; linear momentum
c     is subtracted to the rest of the particles when this happens. 
c
      subroutine getrho (gnr)
c
c     output of the mass density is written in the array
c     rho (i,j,k,m) ; i,j,k = 0,n-1 ; m = number of grids (= 5)
c
      include 'super.cb'
c
      integer gnr,start,istop,i,j,k,l
      integer ix1,ix2,ix3,ix4,ix5
      integer iy1,iy2,iy3,iy4,iy5
      integer iz1,iz2,iz3,iz4,iz5
      real    rcore1,rout1,fnh,xmass,x,y,z,dx,dy,dz,r
      real    rhogrid(5),volcore,volout,volsys, comf 
      real    linmo(3), lincount 

c
c
c     (1.) Initialize rho :
c     =====================
c
      do  l = 1,5
         rhomax(gnr,l) = -1. 
         rhogrid(l)    =  0.  
         do  k = 0,n-1
            do j = 0,n-1
               do i = 0,n-1
                  rho(i,j,k,l) = 0.
               enddo
            enddo
         enddo
      enddo
c      
c     (2.) Finding the number of the first and the last star of
c          the galaxy number gnr :
c     =========================================================
c
      if (gnr.eq.1) then
         start = 1
         istop = gstno(1)
      else
         start = 0
         do i = 1,gnr-1
            start = start + gstno(i)
         enddo
         istop  = start + gstno(gnr)
         start  = start + 1
      endif
c
c     (3.) Initialize the quadratic radii of the different grids of
c          the galaxy gnr :
c     =============================================================
c
      rout1  = rout(gnr)*rout(gnr)
      rcore1 = rcore(gnr)*rcore(gnr)
      xmass  = mass(gnr)
      fnh    = real(n/2)
c
c       Here's a common geometric factor, where the sampling is taken 
c       to be spherical : 

      comf = 4.*pi / 3. 

      volcore = rcore1 * rcore(gnr) * comf 
      volout  = rout1 *rout(gnr) * comf 
      volsys  = rsystem(gnr) * rsystem(gnr) * rsystem(gnr) * comf
c      
c     (4.) Looking for each star of galaxy gnr inside of which grid
c          and in which mesh point of these grids and adding the
c          mass of one star to that point of array rho :
c     ==============================================================
c
      do i = 1,3

         linmo(i) = 0.
      end do 

c     Set counter = 1 to avoid division by zero later. 

      lincount = 0.0 
            
c
      do 20 i = start,istop ! beginning star loop
c         
         if (star(1,i).gt.1.e+10) goto 20 ! star is outside of
                                          ! simulation range
c
c     (4.1) Finding postion of the star :
c     ===================================
c
         x   = star(1,i)              ! position of star i
         y   = star(2,i)
         z   = star(3,i)
c
         dx  = x - dgcms(1,gnr)       ! relative position of star i
         dy  = y - dgcms(2,gnr)       ! to the center of galaxy gnr
         dz  = z - dgcms(3,gnr)

c     (4.2) Finding the integer numbers of the mesh-points of the
c           star
c     ===========================================================
c
         ix1 = nint(enh1(gnr)*dx + fnh) ! middle grid without inner
         ix2 = ix1                      ! middle grid with inner 
         ix3 = nint(enh3(gnr)*dx + fnh) ! inner  grid
         ix4 = nint(enh4(gnr)* x + fnh) ! outer  grid with middle
         ix5 = ix4                      ! outer grid without middle
         iy1 = nint(enh1(gnr)*dy + fnh)
         iy2 = iy1
         iy3 = nint(enh3(gnr)*dy + fnh)
         iy4 = nint(enh4(gnr)* y + fnh)
         iy5 = iy4
         iz1 = nint(enh1(gnr)*dz + fnh)
         iz2 = iz1
         iz3 = nint(enh3(gnr)*dz + fnh)
         iz4 = nint(enh4(gnr)* z + fnh)
         iz5 = iz4

c
c     (4.3) Proofing if star is still inside the simulation :
c     =======================================================
c
         if( (ix4.lt.2).or.(ix4.ge.(n-2)) 
     &   .or. (iy4.lt.2).or.(iy4.ge.(n-2))
     &   .or. (iz4.lt.2).or.(iz4.ge.(n-2)) ) then

            star(1,i) = 1.e+30

            linmo(1) = linmo(1) + xmass*star(4,i)
            linmo(2) = linmo(2) + xmass*star(5,i)
            linmo(3) = linmo(3) + xmass*star(6,i) 
            lincount = lincount + 1. 

            goto 20
         endif
c++         if ((iy4.lt.2).or.(iy4.ge.(n-2))) then
c            star(1,i) = 1.e+30
c            goto 20
c         endif
c         if ((iz4.lt.2).or.(iz4.ge.(n-2))) then
c            star(1,i) = 1.e+30
c            goto 20
c         endif
c++
c
         r = dx*dx + dy*dy + dz*dz ! calculating the quadr. radius
                                   ! of star i
c
c     (4.4) Adding the particle mass to the corresponding mesh-points :
c     =============================================================
c
c     if the particle is outside of the range of the middle grid,
c     it is only calculated in the outer grid without the middle
c     grid window (grid 5) :
c
         rhogrid(5) = rhogrid(5) + xmass/volsys ! Average system density
 
         if (r.gt.rout1) then
            rho(ix5,iy5,iz5,5) = rho(ix5,iy5,iz5,5) + xmass
            rhogrid(3) = rhogrid(3) + xmass/(volsys-volout)
            goto 20
         endif
         
c         
c     else it is calculated in the outer grid which corresponds
c     to the middle grid (grid 4) :
c
         rho(ix4,iy4,iz4,4) = rho(ix4,iy4,iz4,4) + xmass
         rhogrid(4) = rhogrid(4) + xmass/volout ! Average grid-4 density
c        
c     if the particle is outside of the range of the inner grid,
c     it is calculated in the middle grid without the inner grid
c     window (grid 1) :
c
         
         if(r.gt.rcore1) then
            rho(ix1,iy1,iz1,1) = rho(ix1,iy1,iz1,1) + xmass
            rhogrid(2) = rhogrid(2) + xmass/(volout-volcore) 
            goto 20
         endif

c
c     else it is calculated in the middle grid which corresponds
c     to the inner grid (grid 2) :
c
         rho(ix2,iy2,iz2,2) = rho(ix2,iy2,iz2,2) + xmass
         rhogrid(1) = rhogrid(1) + xmass/volcore ! Average grid-1 density

c
c     and the particle is calculated in the innermost grid (grid 3) :
c
         rho(ix3,iy3,iz3,3) = rho(ix3,iy3,iz3,3) + xmass

c
 20   continue


c
c     Add linear momentum of lost particles to rest of the system
c      

c       ------------------- Removed because of malfunction. RE 5/Nov/13

c      if( lincount .gt. 0 ) then 

c      do i = start,istop 
c
c         star(4,i) = star(4,i) + linmo(1)/xmass/lincount 
c         star(5,i) = star(5,i) + linmo(2)/xmass/lincount
c         star(6,i) = star(6,i) + linmo(3)/xmass/lincount 
c      end do 
c      endif 
c       ------------------- -------------------------------------------


c
c       Store the integrated densities for each galaxy - 
c
      do i = 1,5 
         rhomax(gnr,i) = rhogrid(i) 
      enddo 
c
c++      write(unit,'(a,$)' ) ' In GETRHO : grid densities = '
c++      write(unit,90210 ) ( i, rhomax(gnr,i), i = 1,5 ) 
#ifdef FORT90 
90210 format( 5(I2,1x,es8.2) ) 
#else 
90210 format( 5(I2,1x,e8.3) ) 
#endif 
      return
c
      end
c
c******************************************************************
c     end of file
c******************************************************************

      







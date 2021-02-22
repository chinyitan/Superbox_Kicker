      program orb_dynfric
C Based on comporb.f.
*
C Computes orbits in an isothermal halo under dynamical friction
C (Binney & Tremaine).
*
C Integration is via Runga-Kutta method.
*
      implicit none
*
      integer                i,j,n,maxline,include
      parameter              (maxline=3000)
C store every nth step
      real                   st
      real                   sx,sy,sr,svx,svy,svel
      real                   sforce,seta,sdens,scoul
*
      real*8                 tn,xn,yn,vxn,vyn,h,tmax,tmin,r,vel
      real*8                 a,b,c,d
      real*8                 f,eta,force,dens
      real*8                 kx1,kx2,kx3,kx4
      real*8                 ky1,ky2,ky3,ky4
      real*8                 kvx1,kvx2,kvx3,kvx4
      real*8                 kvy1,kvy2,kvy3,kvy4
      real*8                 Vx,Vy,Xp,Yp
      real*8                 gamma,Vc,Msat,G,Mgal
      real*8                 time,crosstime
      real*8                 veltransf,lengthtransf,masstransf
      real*8                 timetransf
      real*8                 gammastar,Vcstar,Msatstar,Gstar,Mgalstar
      real*8                 timestar,crosstimestar
      real*8                 coulomb,bmax,bmin,hole
      real*8                 centrtime,alpha,factor
      real*8                 pi
      character*80           junk,filename
      character*1            f_b_integr

      parameter(pi=3.1415927D0)
*
*
* ===========Model PARAMETERS: =========
* Gravitational constant:
      Gstar = 2.D0
* core radius of infinitely extended isothermal halo:
      gammastar = 2.D0
* Mass scale:
      Mgalstar = 1.D0
*
*======INITIAL CONDITIONS: ======
* Read parameters and initial conditions in km/sec, pc, solar masses:
      read(5,'(a)')junk
      write(6,*)' forward/backward integration? (f/b)'
      read(5,'(a)')f_b_integr
      write(6,'(a)')f_b_integr
      read(5,'(a)')junk
      write(6,*)' step size in Gyr?'
      read(5,*)h
      write(6,*)h
      write(6,*)' store every nth step: n?'
      read(5,'(I6)')n
      write(6,'(I6)')n
      write(6,*)' start time ? [Gyr]'
      read(5,*)tn
      write(6,*)tn
      tmin = tn
      write(6,*)' max. computing time?'
      read(5,*)tmax
      write(6,*)tmax
      read(5,'(a)')junk
*
      write(6,*)' G [pc**3/(solar masses*yr**2)] ?'
      read(5,*) G
      write(6,*) G
      write(6,*)' gamma [kpc]?'
      read(5,*) gamma
      write(6,*) gamma
      read(5,'(a)')junk
      read(5,'(a)')junk
*
      gammastar = gammastar * gamma/3.D0
*
      write(6,*)' initial Vx (km/sec)?'
      read(5,*)Vx
      write(6,*)Vx
      write(6,*)' initial Vy (km/sec)?'
      read(5,*)Vy
      write(6,*)Vy
      if (f_b_integr.EQ.'b') then
         Vx = -1.D0*Vx
         Vy = -1.D0*Vy
      end if
      write(6,*)' initial X (kpc)?'
      read(5,*)Xp
      write(6,*)Xp
      write(6,*)' initial Y (kpc)?'
      read(5,*)Yp
      write(6,*)Yp
      write(6,*)' circular halo velocity (km/sec)?'
      read(5,*)Vc
      write(6,*)Vc
      write(6,*)' Mass of satelite (solar masses)?'
      read(5,*)Msat
      write(6,*)Msat
      write(6,*)' bmin (rad. of sat.) [kpc]?'
      read(5,*)bmin
      write(6,*)bmin
      read(5,'(a)')junk
      write(6,*)' include dynamical friction? (y=1,n=0)'
      read(5,*) include
      write(6,*) include
      read(5,'(a)')junk
*
      write(6,*)' initial hole?'
      read(5,*)factor
      write(6,*)factor
      write(6,*)' alpha?'
      read(5,*)alpha
      write(6,*)alpha
      read(5,'(a)')junk
*
      write(6,*)' filename?'
      read(5,'(a)')filename
      write(6,'(a)')filename
*
* Mgal in solar masses for an isothermal halo within a radius of 100kpc:
      Mgal = 5.7D11*(Vc/220.D0)**2.D0 * (100.D0/50.D0)
* transform velocities to to pc/yr and distances to pc
      Vx = Vx*1.02274D-6
      Vy = Vy*1.02274D-6
      Vc = Vc*1.02274D-6
      Xp = Xp*1000.D0
      Yp = Yp*1000.D0
      gamma = gamma*1000.D0
      bmin = bmin * 1000.D0
*================================
* Transform to model units:
      veltransf = DSQRT((Gstar*Mgalstar*gamma)/(G*Mgal*gammastar))
* in km/sec:
      Vcstar = veltransf * Vc
      vxn = veltransf * Vx
      vyn = veltransf * Vy
      crosstimestar = 2.D0*DSQRT(2.D0) * gammastar/Vcstar
      crosstime =  2.D0*DSQRT(2.D0) * gamma/Vc
      lengthtransf = ((Gstar*Mgalstar*crosstimestar*crosstimestar) /
     &                (G * Mgal * crosstime*crosstime))**(1.D0/3.D0)
      xn = lengthtransf * Xp
      yn = lengthtransf * Yp
      bmin = lengthtransf * bmin
      masstransf = ((gammastar**3.D0)*G*crosstime*crosstime) /
     &              ((gamma**3.D0)*Gstar*crosstimestar*crosstimestar)
      Msat = masstransf * Msat
      timetransf = (gammastar*Vc)/(Vcstar*gamma)

      write(6,*)
      write(6,*)' In model units:'
      write(6,*)' xn=',xn
      write(6,*)' yn=',yn
      write(6,*)' vxn=',vxn
      write(6,*)' vyn=',vyn
      write(6,*)' Vc=',Vcstar
      write(6,*)' Mgal=',Mgalstar
      write(6,*)' Msat=',Msat
      write(6,*)' bmin=',bmin
      tmax = tmax * timetransf*1.D9
      tmin = tmin * timetransf*1.D9
      tn = tn * timetransf*1.D9
      h = h * timetransf*1.D9
      write(6,*)' tmax= ',tmax
      write(6,*)' tmin= ',tmin
      write(6,*)' tmax-tmin= ',tmax-tmin
      write(6,*)' h= ',h
      write(6,*)
*
      open(unit=11,file=filename)
*
* For Coulomb term: lambda = bmax/bmin, take bmin to be radius of satellite.
* Uncertain is bmax:
*   taking it to be instantaneous radial distance of satellite from origin, 
*      as in Tremaine (1976), leads to serious oscillations in vel. once satellite
*      reaches centre.
*   take it to be the initial radial distance of satellite from origin, 
*      bmax is thus not time dependent but a constant:
c       bmax = DSQRT(xn*xn+yn*yn)
*
c      bmin = 0.1D0
*
c      n = 20*NINT(tmax/maxline)
      i=0
      j = 0
      centrtime = 0.D0
      hole = factor
*
*
*
      do while (DABS(tn).LE.DABS(tmax))
         j=j+1
*
* See above comments on bmax:
         bmax = 0.1D0*DSQRT(xn*xn+yn*yn)
         coulomb = DLOG(bmax/bmin)
c         coulomb = 1.7D0
*
         if (coulomb.LE.0.D0) coulomb = 0.D0
* model emptying of central region of halo due to satellite:
* here central region will be:
*    1) emtied only if satellite actually passes through central region
*    2) the void will grow with time
         if (bmax.LE.hole) then
            centrtime = DABS(centrtime + h)
            hole = factor * centrtime**alpha
         end if
         if (bmax.LE.hole) coulomb=0.D0
*
         tn = tn + h
*
         kx1 = vxn
         kx2 = vxn + 0.5D0*h*kx1
         kx3 = vxn + 0.5D0*h*kx2
         kx4 = vxn + h*kx3
*
         ky1 = vyn
         ky2 = vyn + 0.5D0*h*ky1
         ky3 = vyn + 0.5D0*h*ky2
         ky4 = vyn + h*ky3
*
         a = xn
         b = yn
         c = vxn
         d = vyn
         kvx1 = f(c,a,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
         kvy1 = f(d,b,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
*
         a = xn + 0.5D0*h*kx1
         b = yn + 0.5D0*h*ky1
         c = vxn + 0.5D0*h*kvx1
         d = vyn + 0.5D0*h*kvy1
         kvx2 = f(c,a,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
         kvy2 = f(d,b,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
*
         a = xn + 0.5D0*h*kx2
         b = yn + 0.5D0*h*ky2
         c = vxn + 0.5D0*h*kvx2
         d = vyn + 0.5D0*h*kvy2
         kvx3 = f(c,a,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
         kvy3 = f(d,b,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
*
         a = xn + h*kx3
         b = yn + h*ky3
         c = vxn + h*kvx3
         d = vyn + h*kvy3
         kvx4 = f(c,a,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
         kvy4 = f(d,b,a,b,c,d,
     &            Gstar,Vcstar,gammastar,Msat,coulomb,tn,include,
     &            f_b_integr)
*
         xn = xn + (h/6.D0) * (kx1 + 2.D0*kx2 + 2.D0*kx3 + kx4)
         yn = yn + (h/6.D0) * (ky1 + 2.D0*ky2 + 2.D0*ky3 + ky4)
         vxn = vxn + (h/6.D0) * (kvx1 + 2.D0*kvx2 + 2.D0*kvx3 + kvx4)
         vyn = vyn + (h/6.D0) * (kvy1 + 2.D0*kvy2 + 2.D0*kvy3 + kvy4)
*
         r = DSQRT(xn*xn+yn*yn)
         vel = DSQRT(vxn*vxn+vyn*vyn)
c         if (vel.LE.1.D-4) then
c            write(6,*)
c            write(6,*)'--------vel<1E-2--------'
c            write(6,*)
c            goto 1000
c         end if
         if (mod(j,n).EQ.0) then
            dens = Vcstar*Vcstar /
     &             (4.D0*pi*Gstar*(gammastar*gammastar+r*r))
c*
c* Write: time frac. comp.,time,r:',
c            write(6,100)j/tmax,tn,r,vel,hole,coulomb
c100         format(6(1x,F9.4))
* Write:',
            st = SNGL(tn)
            if (f_b_integr.EQ.'b') st = -1.*st
            sr = SNGL(r)
            sx = SNGL(xn)
            sy = SNGL(yn)
            svel = SNGL(vel)
            seta = SNGL(eta(xn,yn,vxn,vyn,
     &                  Gstar,Vcstar,gammastar,Msat,coulomb))
            scoul = SNGL(coulomb)
            sforce = SNGL(force(xn,yn,Gstar,Vcstar,gammastar))
            sdens = SNGL(dens)
            svx =SNGL(vxn)
            svy =SNGL(vyn)
* Transform to physical units (Gyr,km/sec,kpc):
* In years:
c            st = 1.E6*st/timetransf
* In Gyr:
            st = st/(timetransf*1.E9)
* In pc/yr:
            svx = svx/(veltransf*1.02274E-6)
            svy = svy /(veltransf*1.02274E-6)
            svel = svel/(veltransf*1.02274E-6)
* In kpc:
            sx = sx/(lengthtransf*1000.)
            sy = sy/(lengthtransf*1000.)
            sr = sr/(lengthtransf*1000.)
*
* Write: fraction of time computed, time, r, etc.:',
            write(6,100)(tn-tmin)/(tmax-tmin),st,sr,svel,hole,coulomb
100         format(6(1x,F9.4))
            write(11,1100) st,sr,sx,sy,svel,svx,svy,
     &                     seta,sforce,sdens,scoul
         end if
      end do
*
1000  continue
      write(6,*)
      write(6,*)'           THE END'
1100  format(1x,F12.3,6(1x,F12.5),3(1x,E12.4),1x,F6.3)
      close(11)
*
      END
C------------------------------------------------------
      real function f(alphaf,betaf,x,y,vx,vy,
     &                G,Vc,gamma,Msat,coulomb,tn,include,f_b_integr)

*
      implicit none
      integer        include
      real*8         alphaf,betaf
      real*8         x,y,vx,vy
      real*8         tn,pk1,pk2,pk3
      real*8         Vc,gamma,veldisp,eta
      real*8         force,msat
      real*8         pi,coulomb,G
      character*1    f_b_integr
      parameter(pi=3.1415927D0)
*
*
C If integrate backward (f_b_integr='b') then dyn. friction becomes an
C acceleration.
      if (include.EQ.1) then
         if (f_b_integr.EQ.'f') then
            f = -1.D0*eta(x,y,vx,vy,G,Vc,gamma,Msat,coulomb) * alphaf
         else if (f_b_integr.EQ.'b') then
            f = eta(x,y,vx,vy,G,Vc,gamma,Msat,coulomb) * alphaf
         end if
         f = f - force(x,y,G,Vc,gamma) * betaf
      else
         f = -1.D0 * force(x,y,G,Vc,gamma) * betaf
      end if

*
c-----for testing:------------
c      if (tn.GE.118.9D0) then
c         pk1 = eta(x,y,vx,vy,G,Vc,gamma,Msat,coulomb)
c         pk2 = force(x,y,G,Vc,gamma)
c         pk3 = (DSQRT(vx*vx+vy*vy))**3.D0
c         write(6,*)pk1,pk2,pk3
c      end if
c-----------------------------

*
      return
      end
*
C---------------------------------------------------------
      real function eta(x,y,vx,vy,G,Vc,gamma,Msat,coulomb)
* dynamical friction
*
      implicit none
      real*8         x,y,vx,vy,r,vel
      real*8         Vc,gamma,veldisp
      real*8         dens,msat
      real*8         pi,coulomb,G,xerf
      parameter(pi=3.1415927D0)
*
      r = DSQRT(x*x+y*y)
      vel = DSQRT(vx*vx+vy*vy)
      if (Vc.EQ.0.D0) then
         veldisp = 1.D0
      else
         veldisp = Vc/DSQRT(2.D0)
      end if
      xerf = vel/(DSQRT(2.D0)*veldisp)
*
* density of halo at distance r:
      dens = Vc*Vc/(4.D0*pi*G*(gamma*gamma+r*r))
*
* dynamical friction:
      if (vel.LE.1.D-5) then
         eta = 0.D0
      else
         eta = 4.D0*pi*coulomb * G*G* dens * Msat
         eta = eta * (ERF(xerf)
     &         - 2.D0*xerf / (DSQRT(pi)*DEXP(xerf*xerf)))
         eta = eta / (vel*vel*vel)
      end if
*
      return
      end
*
C-----------------------------------------------------------
      real function force(x,y,G,Vc,gamma)
* gravitational attraction in isothermal infinite halo:
*
      implicit none
      real*8         x,y,r
      real*8         Vc,gamma
      real*8         pi,G
      parameter(pi=3.1415927D0)
*
      r = DSQRT(x*x+y*y)
*
* gravitational accelaration (force/mass) in isothermal infinite halo:
      force = Vc*Vc/(gamma*gamma+r*r)
*
      return
      end












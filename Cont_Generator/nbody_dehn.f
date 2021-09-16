      implicit real*8 (a-h,o-z)
      real*8 mu
      parameter (mu=5.6e10,ru=3.5,vu=262.)
      parameter (ne=250,nr=4000,np=120)
      parameter (nb=5e4)
      real*8 xx(np),ww(np)
      real*8 rf(nr),qf(nr),fst(nr),fdm(nr),mle(nr)
      real*8 rad(3,nb),vel(3,nb),Ms,ml(nb),mlaux(nb),mpart
      integer listml(nb),idum
      external phinfw,phicore,phihern,rhonfw,rhocore,rhohern
      external phidehn,rhodehn
      external phicl
      character*60 namenb,nameasc,namelist
      character namemod*18
!      character namemod*10
      character namebar*23                !namemod*26
      character namedir*7
      common gh,rsh,rho0h
      pi=4.d0*atan(1.d0)
      
c     this code generates equilibrium ensembles of particles (Dehnen profile) moving in a potential
c     potential = self-grav Dehnen (stars) + external Dehnen (DM)
      
c     N-body units: stellar component
      rho0=1.d0
      rs=1.d0
      
      
c     DM Dehnen (1993)
      dmdm=            1.4e9 /mu  !total mass 
      rsdm=                1.6/ru !scale radius (kpc)
      gh  =  0.d0               !slope 
c     stars Dehnen (1993)
      dmp=      1.35e3/mu
      drp=         0.78  /1000.   /ru  !kpc
      g   =  0.d0
      
      namemod='m1e3r05_m1e9r10g00'
      rho0dm=dmdm*(3.d0-gh)/(4*pi*rsdm**3.)
      rho0p=dmp*(3.d0-g)/(4*pi*drp**3.)
      
      print*,'Mhalo (msol)=',dmdm*mu,log10(dmdm*mu)
      print*,'rhalo (kpc)=',rsdm*ru      

      print*,'Mstars (msol)=',dmp*mu,log10(dmp*mu)
      print*,'rstars (pc)=',drp*ru*1000
            
c     DM in units of stars
      rho0h=rho0dm/rho0p
      rsh  =rsdm/drp
      print*,'dm rho0,rs in N-body units ',rho0h,rsh
      
      
c     stars alpha-beta-gamma;  Plummer (2,5,0), Hern (1,4,1)
      as=2.
      bs=5.
      gs=0.1


      rc=rs* 1.d0
      ra=rc* 1.d20      
c     
      idum=124869

      namedir='./mods/'

      namenb =namedir//namemod//'.CONT'
      nameasc=namedir//namemod//'.dat'

      print*,namemod

c     DF
      rmax=100.       *rs  
      rmin=1.e-5     *rs

c      open(1,file='df.dat')       
c      open(2,file='prof.dat')

      call df(phicl,rhodehn,rmin,rmax,ne,g,rs,rho0,rc,ra,as,bs,gs,
     &     nr,rf,qf,fst,fdm,mle,io)
      call nbody(idum,phicl,rhodehn,rmin,rmax,
     &     g,rs,rho0,rc,ra,as,bs,gs,
     &     nr,rf,qf,fdm,mle,nb,rad,vel,ml)   
!      call prof(phidehn,rhodehn,rmin,rmax,g,rs,rho0,rc,ra,as,bs,gs,
!     &     nr,rf,qf,fdm,nb,rad,vel)

      close(1)
      close(2)


c     out ********
      
c     masa satelite
      call gauleg(rmin,rmax,xx,ww,np)
      dmtot=0.d0
      do i=1,np
         rhot=rhodehn (xx(i),g,rs,rho0)
         dmtot=dmtot+ww(i)*rhot*4.*pi*xx(i)**2.
      enddo
      Ms=dmtot! masa en N-body
      mpart=dmtot/dble(nb) ! masa una particula

      dmhalo=4.*pi*rho0/(3.d0-g)*rs**3.
      print*,'Mhalo=',Ms,dmhalo


c     out particulas ASCII / ML list
      open(3,file=nameasc)
      write(3,*) nb,mpart,dmp,drp,g
      do i=1,nb  
!         write(3,100) rad(1,i),rad(2,i),rad(3,i),
!     &        vel(1,i),vel(2,i),vel(3,i)
         write(3,100) rad(1,i),rad(2,i),rad(3,i),
     &        vel(1,i),vel(2,i),vel(3,i)
         ml(i)=1.d0
      enddo
      close(3)

c     SUPERBOX
      call superbox(namenb,Ms,g,rho0,rs,rc,
     &     nb,rmax,rad,vel,ml,dmp,drp)
 100  format(7e12.4) 
      end

      subroutine nbody(idum,potfun,rhofun,rmin,rmax,
     &     g,rs,rho0,rc,ra,as,bs,gs,
     &     nr,rf,qf,pf,mle,nb,rad,vel,ml)
      implicit real*8 (a-h,o-z)
      real*8 rf(nr),qf(nr),pf(nr),mle(nr)
      real*8 rad(3,nb),vel(3,nb),rmod(nb),ml(nb)
      integer list(nr),listmax(nr),idum
      real*8 pnmax(nr), er(3),ep(3),et(3)
      real*8 rhomax(nr)
      parameter (nmax=5e7)
      external rhon,phin,rhoq,phinfw,rhonfw,rtbis,potfun,rhofun,herndf
      external ran3
      pi=4.d0*atan(1.d0)
      twopi=2.d0*pi
c     
      print*,'ra/rc=',ra/rc
      print*,'qmin,qmax=',qf(1),qf(nr)
      print*,'rmin,rmax',rmin,rmax
      print*,'f_min,f_max=',pf(1),pf(nr)


c     GENERANDO POSICIONES
      dr=(rmax-rmin)/real(nr)
      rr=rmin-dr
      do i=1,nr
         rr=rr+dr
         rqa=1.+rr**2./ra**2.
         rhomax(i)=rr**2.*rhofun(rr,g,rs,rho0)/rqa

         listmax(i)=i
         
      enddo
      call sort2(nr,rhomax,listmax)
      fmax=rhomax(nr)*1.1
      print*,'rho*r^2 max=',fmax,ra/rc,listmax(nr)

      k10=int(nb/10)
      xra=0.d0
      do i=1,nb
         j=0
 60      continue
c     r
         xr=ran3(idum)*(rmax-rmin)+rmin
         rqa=1.+xr**2./ra**2.

         pn=rhofun(xr,g,rs,rho0)/rqa* xr**2./fmax


         f=ran3(idum)
         if (pn.gt.1.) then
            print*,'error p',pn,i,xr
            stop
         endif

         if (f.le.pn) then
            rmod(i)=xr
            xra=xra+xr/real(nb)
            goto 50
         else
            j=j+1
            if (j.ge.nmax) then
               print*,'error en generacion de p(x) ',j,nmax
               j=0
               goto 60
!               end
            endif
            goto 60
         endif
 50      continue
c     k10 print i=1,nb incrementos de 10
         if (i.eq.k10) then
            print*,k10
            k10=k10+int(nb/10)
         endif
      enddo
      print*,'<xr>=',xr

c     GENERANDO VELOCIDADES      
      do i=1,nr
         rqa=1.+(rf(i)/ra)**2.
         rhomax(i)=pf(i)*sqrt(2.*qf(i)) *rf(i)**2./rqa
         listmax(i)=i
      enddo
      call sort2(nr,rhomax,list)
      fmax=rhomax(nr)

c     boost
      xi=600.
 66   pnboost=1./fmax*xi        !(rs/rc)**1.5
      print*,'pnboost=',pnboost,xi

      Phimax=potfun(rmin,g,rs,rho0)       
      vemax=sqrt(2.d0*Phimax)

      k10=int(nb/10)
      ibad=0
      do i=1,nb       
         j=0
c     r
         xr=rmod(i)

         Phi=potfun(xr,g,rs,rho0)
         ve=sqrt(2.d0*Phi)
         rqa=(1.+xr**2./ra**2.)

 20      continue

c     Q
         xp=qf(1)+ran3(idum)*(Phi-qf(1))      
c     beta
         beta=(xr/ra)**2./rqa
c     theta, angle vr=v*cth, vt=v*sth
         dxrand=ran3(idum)
         cth=2.0*(dxrand-0.5)
         sth=SQRT(1.0-cth*cth)
         dxrand=ran3(idum)
         signs=2.0*(dxrand-0.5)
         cth=signs*cth/ABS(signs)
         phiran=twopi*dxrand

c     vel
        viso=sqrt( 2.*(Phi-xp) )
        vr=viso*cth
        vt=viso*sqrt(1.d0-beta)*sth
        vth =vt*cos(phiran)
        vphi=vt*sin(phiran)
        v=sqrt(vr**2.+vth**2.+vphi**2.)
        
        if (v.ge.ve) print*,'sdfsdf'
c     f(Q) interpolacion
         qq=xp

         do k=nr,1,-1
            if (qf(k).le.qq) then
               pfi=pf(k+1)+(pf(k)-pf(k+1))/(qf(k)-qf(k+1))*(qq-qf(k+1))
c               print*,'int ',k,qf(k+1),qq,qf(k),pf(k+1),pfi,pf(k),Phi
c                print*,xr,pfi/herndf(qq,rs,ra,rho0) !hernquist
c               stop
               goto 30
            endif
         enddo
         print*,'error interp. DF ',qq,qf(1),qf(nr),pf(1),pf(nr),y
         stop
 30      continue

c     f(Q) sqrt(2(Q-Phi)) r^2  --> likelihood

         pn=pfi*sqrt(2.*(Phi-qq)) *pnboost*xr**2./rqa


         f=ran3(idum)

         if (pn.gt.1.) then
            print*,'error pv',pn,i
            xi=xi/2.
            goto 66
         endif    

         if (f.le.pn) then
c     radius
            dxrand=ran3(idum)
            cth=2.0*(dxrand-0.5)
            sth=SQRT(1.0-cth*cth)
            dxrand=ran3(idum)
            signs=2.0*(dxrand-0.5)
            cth=signs*cth/ABS(signs)
            dxrand=ran3(idum)
            phi=twopi*dxrand

            rad(1,i)=xr*sth*cos(phi)
            rad(2,i)=xr*sth*sin(phi)
            rad(3,i)=xr*cth         
c     velocity
            vel(1,i)= vr*sth*cos(phi)+vth*cth*cos(phi)-
     &           vphi*sin(phi)
            vel(2,i)= vr*sth*sin(phi)+vth*cth*sin(phi)+
     &           vphi*cos(phi)
            vel(3,i)= vr*cth-vth*sth
!            print*,xr,beta,(xr/ra)**2./rqa,
!     &           (vr**2.+vt**2.)/(vr**2.+vth**2.+vphi**2.),vr

            goto 40
         else
            j=j+1            
            if (j.ge.nmax) then
               print*,'error en generacion de pv(x) ',i,j,nmax,xr/rc,pn
               rad(1,i)=1.d10
               rad(2,i)=1.d10
               rad(3,i)=1.d10
               vel(1,i)=1.d10
               vel(2,i)=1.d10
               vel(3,i)=1.d10
               ibad=ibad+1
               goto 40
!               j=0
!               goto 20
               stop
            endif
            goto 20
         endif
 40      continue

c     k10 print i=1,nb incrementos de 10
!         print*,i,pn,ibad
         if (i.eq.k10) then
            print*,k10
            k10=k10+int(nb/10)
         endif
      enddo
c     
      print*,'part generadas ',nb,ibad
      return
      end




      subroutine prof(potfun,rhofun,rmin,rmax,g,rs,rho0,rc,ra,as,bs,gs,
     &     nr,rf,qf,pf,nb,rad,vel)
      implicit real*8 (a-h,o-z)
      parameter (np=50)
      real*8 rf(nr),qf(nr),pf(nr)
      real*8 ml(nb),rad(3,nb),vel(3,nb),rr(nb),vv(nb),vr(nb),vt(nb)      
      real*8 ya(np),ys(np),xa(np),xn(np),x1(np),x2(np)
      real*8 vra(np),vta(np),vrs(np),vts(np)
      external rhon,phin,rhoq,phinfw,rhonfw,rtbis,potfun,rhofun
      external srhern
      pi=4.d0*atan(1.d0)
c
      print*,'Prof ',nb,np,nr
c     distrib
      do i=1,nb
         ml(i)=1.d0

         rr(i)=sqrt(rad(1,i)**2.+rad(2,i)**2.+rad(3,i)**2.)
         vv(i)=sqrt(vel(1,i)**2.+vel(2,i)**2.+vel(3,i)**2.)
         vr(i) =
     &   (rad(1,i)*vel(1,i)+rad(2,i)*vel(2,i)+rad(3,i)*vel(3,i))/rr(i)
         vt(i) = sqrt(vv(i)**2. - vr(i)**2.)
      enddo


c     
      a=rmin
      b=rmax    
      call ldistrib(a,b,ml,nb,rr,vr,xn,x1,x2,xa,vra,vrs,np)

c
      sr0=2.*pi*rc**2.*rho0* (xa(1)/rc)*log(rc/xa(1))
      do i=1,np
         rqa=1.+xa(i)**2./ra**2.
         sd=(xn(i)+1.e-10)/(4.*pi/3.*(x2(i)**3.-x1(i)**3.) )
!         print*,xa(i)/rc,int(xn(i)),
!     &        log10(sd)-log10(rhoq(xa(i),rc,ra,as,bs,gs))
!         print*,xa(i)/rs,vrs(i)**2.,srhern(xa(i),rc,rho0)
         write(2,*)log10(xa(i)/rc),log10(sd),
     &        log10(rhoq(xa(i),rc,ra,as,bs,gs)/rqa),
     &        vrs(i)**2.,srhern(xa(i),rc,rho0)
      enddo
      return
      end

      subroutine 
     & df(potfun,rhofun,rmin,rmax,n,g,rs,rho0,rc,ra,as,bs,gs,
     &     nr,rf,qf,fst,fdm,ml,io)
      implicit real*8 (a-h,o-z)
      real*8 mu
      parameter (mu=5.6e10,ru=3.5,vu=262.)
      parameter (eps=1.e-4,eps2=1.e-8)
      real*8 xx(n),ww(n)
      real*8 rf(nr),qf(nr),fst(nr),fdm(nr),ml(nr),mlaux(nr),thetaq(nr)
      real*8 dnedm(nr),dnest(nr)
      integer listml(nr)
      external rhon,phin,rhoq,phinfw,rhonfw,rtbis,potfun,rhofun
      external rhohern,phihern,phidehn,rhodehn
      external fqhern
      pi=4.d0*atan(1.d0)

c            
      io=0
      iod=0

c     Intervalos R/rc
      dr=(rmax-rmin)/real(nr)
      qmin=0.d0!potfun(rmax,rs,rho0)
      qmax=potfun(rmin/10.,g,rs,rho0)

      print*,'dr/rc=',dr/rc
      print*,'qmin,qmax ',qmin,qmax
      print*,'rmin,rmax ',rmin,rmax
      print*,'g=',g
      print*,'NE=',n
      print*,'Nr=',nr

c     f(Q), g(Q)
      rneg=0.d0
      rdneg=0.d0
      krdneg=0
      krneg =0
      k10=int(nr/10)
      rr=rmax
      do i=1,nr       
                  
c     f(Q)
         qf(i)=potfun(rr,g,rs,rho0)     
         rf(i)=rr                  
         q1=qf(i)
         
         call gauleg(qmin,q1,xx,ww,n)
         didm=0.d0
         dist=0.d0
         do j=1,n
            q=xx(j)

c     r=r(q) Intervalo r
            rmin2=eps2       *rs
            rmax2=1000000000.*rs
            r=rtbis(q,potfun,g,rs,rho0,rmin2,rmax2,eps2)
            

            p=potfun(r,g,rs,rho0)         
            d=rhofun(r,g,rs,rho0)
            dq= rhoq(r,rc,ra,as,bs,gs)
            

            r2=r*(1.+eps)
            p2=potfun(r2,g,rs,rho0)
            d2=rhofun(r2,g,rs,rho0)
            dq2= rhoq(r2,rc,ra,as,bs,gs)

            r1=r*(1.-eps)
            p1=potfun(r1,g,rs,rho0)
            d1=rhofun(r1,g,rs,rho0)
            dq1= rhoq(r1,rc,ra,as,bs,gs) 
         
c     deriv
            dx=(r2-r1)/2.
            dphir  =(p2-p1)/(2.*dx)
            drhor  =(d2-d1)/(2.*dx)
            drhoqr =(dq2-dq1)/(2.*dx)

            drho2r2 =(d2+d1-2.*d)/(dx)**2.
            drhoq2r2=(dq2+dq1-2.*dq)/(dx)**2.
            dphi2r2 =(p2+p1-2.*p)/(dx)**2.

            drhop =drhor /dphir
            drhoqp=drhoqr/dphir

            drho2phi2 =1./dphir**3.*(drho2r2 *dphir-drhor *dphi2r2)
            drhoq2phi2=1./dphir**3.*(drhoq2r2*dphir-drhoqr*dphi2r2)
            if (j.eq.n) thetaq(i)=drhoq2r2 - drhoqr * dphi2r2/dphir

c     isothermal potential theta teorica
            thetaqt=(1.+r**as)**((gs-bs-2.)/as)/r**(2.+gs)*
     &      (gs**2.+r**as*(gs*(as+2.*bs)-as*bs)+bs*r**(2.*as))
c            
c
            if (drhoq2phi2.lt.0.d0.and.r.le.50.*rc) then
               iod=1
               if (krdneg.eq.0) then
                  krdneg=1
                  rdneg=rr
               endif
             endif
             if (j.eq.n) derivq=drhoq2phi2

c     rho.dphi at q=0
            if (i.eq.1.and.j.eq.1) then
               drhop0 =drhop
               drhoqp0=drhoqp
            endif
!            print*,q1,q,r/rm,rm,p/q,drho2phi2

c     Integra f(Q)=1/sqrt(8.)/pi**2 * (d2rhoq/dp2 /sqrt(q-p) + 1/sqrt(q)drHodp=0)
            fqdm=1./sqrt(8.d0)/pi**2.*
     &           (drho2phi2 /sqrt(q1-q) + drhop0 /sqrt(q)) 
            fqst=1./sqrt(8.d0)/pi**2.*
     &           (drhoq2phi2/sqrt(q1-q) + drhoqp0/sqrt(q)) 

            didm=didm+ww(j)*fqdm
            dist=dist+ww(j)*fqst
            
            if (fqdm.ne.fqdm) then
               print*,'DF NaN ',j,drho2phi2,q1-q,qmin,q1,q
               stop
            endif

         enddo
c     DF
         fst(i)=dist+1.e-20
         fdm(i)=didm+1.e-20

c     Integra g(Q)=16pi^2 int_0^r(p) sqrt(2*(p-q))r^2
         call gauleg(rmin,rr,xx,ww,n)
         gf=0.d0
         do j=1,n
            p=potfun(xx(j),g,rs,rho0)     
            gf=gf+ww(j)*16.*pi**2.*sqrt(2.*(p-q))*xx(j)**2.
         enddo
c     N(E)dE
         dnedm(i)=gf*fdm(i)
         dnest(i)=gf*fst(i)
                  
         
c     f(E)<0
         if (dist.lt.0.d0.and.rr.le.50.*rc.and.i.lt.nr) then
            io=1
            if (krneg.eq.0) then
               krneg=1
               rneg=rr
            endif
         endif

c     k10 print i=1,nr incrementos de 10
         if (i.eq.k10) then
            print*,k10,nr,log10(rf(i)/rc)
            k10=k10+int(nr/10)
         endif
c
         rr=rr-dr
      enddo
      print*,'iod,io=',iod,io,rdneg/rc,rneg/rc

c     normalizacion m/l
      dmdm=0.d0
      dmst=0.d0
      do i=1,nr
         dmdm=dmdm+dnedm(i)
         dmst=dmst+dnest(i)
      enddo
      do i=1,nr
         ml(i)=(fst(i)/dmst ) / (fdm(i)/dmdm) 
         mlaux(i)=ml(i)
         listml(i)=i
      enddo
      call sort2(nr,mlaux,listml)
      do i=1,nr
         ml(i)=ml(i)/mlaux(nr)
         if (rf(i).lt.2.*rc)print*,qf(i)/qmax,rf(i)/rc,ml(i),dmst,dmdm
      enddo
      
c     f(E)>0
      if (io.eq.1) then
         do i=1,nr
            print*,rf(i)/rc,fst(i),fdm(i)
         enddo
         print*,'Negative DF'
         stop         
      endif

c     Test
!      call test(potfun,rhofun,g,rs,rho0,rc,ra,as,bs,gs,nr,rf,qf,fst,fdm)      
c     
 199  format (1i6,10e16.6)
      end

 
     
      function rhoq(r,rm,ra,a,b,g)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
      rt=10.*rm !king tidal radius
      xr=sqrt(  (1.+(r/rm)**2.) / (1.+(rt/rm)**2.)  )
      rqa=(1.+r**2./ra**2.)
      rho0=1.d0

      rhoq=rqa*rho0 /(r/rm)**g/(1.+(r/rm)**a)**((b-g)/a)

!      rhoq=rqa*rho0 /xr**2.*(acos(xr)/xr-sqrt(1.-xr**2.))
!      if (r.ge.rt) rhoq=1.e-20

      return
      end

 
      subroutine 
     &     test(potfun,rhofun,g,rs,rho0,rc,ra,as,bs,gs,nr,rf,qf,fst,fdm)
      implicit real*8 (a-h,o-z)
      real*8 rf(nr),qf(nr),fst(nr),fdm(nr)
      parameter (nq=60,eps2=1.e-8)
      real*8 xx(nq),ww(nq)
      external rtbis,potfun,rhofun
      pi=4.d0*atan(1.d0)
      qmin=0.d0
      

      do j=1,nr    

         di=0.
         Phi=qf(j)
         call gauleg(qmin,Phi,xx,ww,nq)

         do i=1,nq
 
c     interpola
            do k=1,nr
               if (xx(i).lt.qf(1)) then
                  fsti=fst(1)
                  q11=qf(1)
                  q22=qf(1)
                  goto 10
               endif
               if (qf(k).ge.xx(i)) then
                  fsti=fst(k-1)+
     &                 (fst(k)-fst(k-1))/(qf(k)-qf(k-1))*
     &                 (xx(i)-qf(k-1))
                  q11=qf(k-1)
                  q22=qf(k)
                  goto 10
               endif
               
            enddo
 10         continue
!            print*,xx(i),q11,q22,fsti
            
            fun=4.*pi*fsti*sqrt(2.*(Phi-xx(i)))
            di=di+ww(i)*fun
         enddo
         dq= rhoq(rf(j),rc,ra,as,bs,gs)
         if (rf(j)/rs.lt.20.) print*,'test ',rf(j)/rs,dq,di,di/dq
      enddo
      
!      r=rtbis(q,potfun,g,rs,rho0,rmin,rmax,eps2)
      return
      end


      function phin(r,rm,vm)
      implicit real*8 (a-h,o-z)
c     jaffe
      phin=-log(r/(r+rm))
      return
      end
      
      function phinfw(r,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.*atan(1.d0)
c     nfw
      phinfw=4.*pi*rho0*rs**2. * log(1.+r/rs)/r
      return
      end

      function phicore(r,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.*atan(1.d0)
c     core
      phicore=4.*pi*rho0*rs**2. * (log(1.+r/rs)/r - 1./2./(rs+r) )
      return
      end

      function phihern(r,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.*atan(1.d0)
c     core
      phihern=2.*pi*rho0*rs**2. / (1.+r/rs)
      return
      end
      
      function rhon(r,rm,vm)
      implicit real*8 (a-h,o-z)
c     jaffe
      rhon=rm**4./r**2./(r+rm)**2.
      return
      end
      
      function rhonfw(r,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
c     nfw*(1+r^2/ra^2)
      rhonfw=rho0/(r/rs)/(1.+r/rs)**2.
      return
      end

      function rhocore(r,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
c     nfw*(1+r^2/ra^2)
      rhocore=rho0/(1.+r/rs)**3.
      return
      end

      function rhohern(r,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
c     nfw*(1+r^2/ra^2)
      rhohern=rho0/(r/rs)/(1.+r/rs)**3.
      return
      end

      function srhern(r,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
c     nfw*(1+r^2/ra^2)
      srhern=pi*rho0*rs**2./6.*(12.*r*(r+rs)**3./rs**4.*log(1.+rs/r)-
     &     r/(r+rs)*(25.+52.*r/rs+42.*(r/rs)**2.+12.*(r/rs)**3.) )
      return
      end

      function herndf(qsn,rs,ra,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
c     
      vg=sqrt(2.*pi*rs**2.*rho0)
      dm=2.*pi*rs**3.*rho0
      q=sqrt(qsn*rs/dm)
c
      herndf=dm/8./sqrt(2.)/pi**3./rs**3./vg**3./(1.-q**2.)**(5./2.)*
     c     (3.*asin(q)+q*sqrt(1.-q**2.)*(1.-2.*q**2.)*
     c     (8.*q**4.-8*q**2.-3.)) + 
     c     dm/(sqrt(2.)*pi**3.*rs**3*vg**3.)*(rs/ra)**2.*q*(1.-2.*q**2.)
      return
      end

      function rhodehn(r,g,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
      rhodehn=rho0/(r/rs)**g/(1.+r/rs)**(4.-g)
      return
      end

      function phidehn(r,g,rs,rho0)
      implicit real*8 (a-h,o-z)
      pi=4.d0*atan(1.d0)
      dm=4.*pi*rs**3.*rho0/(3.-g)
      if (int(g).eq.2) then
         phidehn=dm/rs*(-log(r/(rs+r))  )
      else
         phidehn=dm/rs*(1.d0- (r/(rs+r))**(2.-g) )/(2.-g)
      endif
      return
      end

      function phicl(r,g,rs,rho0)
      implicit real*8 (a-h,o-z)
      common gh,rsh,rho0h      
      pi=4.d0*atan(1.d0)
c     stars
      dm=4.*pi*rs**3.*rho0/(3.-g)
      if (int(g).eq.2) then
         phis=dm/rs*(-log(r/(rs+r))  )
      else
         phis=dm/rs*(1.d0- (r/(rs+r))**(2.-g) )/(2.-g)
      endif
c     dm
      dmh=4.*pi*rsh**3.*rho0h/(3.-gh)
      if (int(g).eq.2) then
         phih=dmh/rsh*(-log(r/(rsh+r))  )
      else
         phih=dmh/rsh*(1.d0- (r/(rsh+r))**(2.-gh) )/(2.-gh)
      endif
c     
      phicl=phis+phih      
      return
      end
      
      subroutine superbox(name,Ms,g,rho0,rs,rc,
     &     nb,rmax,rad,vel,ml,dmp,drp)
      implicit real*8 (a-h,o-z)
      real*8 mu
      parameter (mu=5.6e10,ru=3.5,vu=262.,tu=0.013)
      real*8 rad(3,nb),vel(3,nb),Ms
      integer ih(60)
      real fh(60),xx(nb),yy(nb),zz(nb),uu(nb),vv(nb),ww(nb)
c
      real*8 ml(nb),rr(nb),vv2(nb),vr(nb),vt(nb)
      parameter (np=20)
      real*8 ya(np),ys(np),xa(np),xn(np),x1(np),x2(np),vra(np),vrs(np)
c
      character*(*) name
      pi=4.*atan(1.)

      ngalaxies=1

c     unidades
      dr=drp/rs
      dm=dmp/Ms
      dv=sqrt(dm/dr)

      dt=real( drp**(3./2.)/sqrt(dmp)/100.   )
      print*,'dt (Myr)=',dt*tu*1000.
 
c     SUPERBOX****************
      print*,'n body____'

      open (2,file=name,access='direct',
     $     recl=240)

c     HEADER
      tstep=dt  / sqrt(2.)
      tnum= 5.  /0.013/(tstep*sqrt(2.)) ! Gyr
      print*,'time steps=',tnum,dt,tstep
      
      ih(1) = 1                 ! number of galaxies
      ih(2) = nb
      ih(3)=0                   !tiempo a 0
      ih(4)=int(tnum)           ! last integration step
      ih(5)=1                   ! physical units (0) or model (1) 
      ih(6)=1                   !c.o.d
      ih(7)=50                   !save data each ih(7) steps
      ih(8) = 0                 ! no sticky particle code
      ih(9) = 0                 ! no black hole code : = 0 no BH. 
      ih(10)= 0 
      ih(11)= 1                 ! flag for comoving (1) or fixed (0) coord.
      ih(12)= 0                 ! constant (0) or adjustable (1) dt 
      ih(13)= 1                 ! unit length = pc (0) or kpc (1) 

      do i = 14,60
         ih(i) = 0
      enddo

      write(2,rec=1) ih

      
      print*,'time step (ph units)=',tstep,rs,dmtot
c     reescalando la grid
      fh(1)=drp*2.
      fh(2)=drp*20.
      fh(3)=rmax*4.
      fh(4)=tstep
c     cambiando masa
      fh(5) =real(dmp)          !total mass
      fh(52)=real(dmp)
      fh(6) = real(nb)     ! number of particles (this galaxy)
      fh(7) = real(nb)     ! number of particles left      
      do j = 1,6
         fh( 7+j) = 0.d0      ! center of mass and velocity
         fh(13+j) = 0.d0      ! focus of grids (center of mass)
      end do
      fh(20) = 0.0              ! current time 
      fh(21) = 0.0              ! don't save any star data
      fh(22) = 0.0
      fh(23) = 0.0
      fh(24) = 0.0              ! don't calculate moment of interia
      fh(25) = 0.0
      fh(26) = 0.0
      fh(27) = 0.0              ! don't calc. velocity dispersion
      fh(28) = 0.0
      fh(29) = 0.0
      fh(30) = 0.0              ! don't calculate angular momentum
      fh(31) = 0.0
      fh(32) = 0.0

      fh(33) = 0.0              ! number of stickies (this galaxy)
      fh(34) = 0.0              ! multiplier for sticky grid
      fh(35) = 0.0              ! multiplier for deadtime
      fh(36) = 0.0              ! multiplier for rad. inelasticity
      fh(37) = 0.0              ! multiplier for tang. inelasticity
      fh(38) = 0.0              ! maximum collision distance

      fh(39) = 0.0              ! empty

      fh(40) = 0.0              ! reserved for black hole data
      fh(41) = 0.0
      fh(42) = 0.0
      fh(43) = 0.0
      fh(44) = 0.0
      fh(45) = 0.0
      fh(46) = 0.0
      fh(47) = 0.0
      fh(48) = 0.0
      fh(49) = 0.0
      fh(50)  = 0.0 
      fh(51) = 1.
      fh(53) = 1.
      fh(54) = 1. 
      fh(55) = 1.
      fh(56) = 1. 
      fh(57) = 1.
      fh(58) = 1. 
      fh(59) = 1. 
      fh(60) = 1. 
      do i=1,ngalaxies
         write(2,rec=1+i) fh
      enddo
      close(2)

c     PARTICLES
      open (2,file=name,access='direct',recl=24)

c     NO cambiando posicion de estrellas
      rsum=0.d0
      vsum=0.d0
      do k=1,nb
         i=k
         xx(i)=rad(1,i)*dr
         yy(i)=rad(2,i)*dr
         zz(i)=rad(3,i)*dr

         uu(i)=vel(1,i)*dv *sqrt(2.)
         vv(i)=vel(2,i)*dv *sqrt(2.)
         ww(i)=vel(3,i)*dv *sqrt(2.)

         write(2,rec=10+ngalaxies*10+i)
     &        xx(i),yy(i),zz(i),uu(i),vv(i),ww(i)   
c     
         rsum=rsum+xx(i)**2.+yy(i)**2.+zz(i)**2.
         vsum=vsum+uu(i)**2.+vv(i)**2.+ww(i)**2.
!         print*,sqrt(xx(i)**2.+yy(i)**2.+zz(i)**2.)
      enddo
      close(2)
      print*,'<r>,<v> ',sqrt(rsum/real(nb))*ru,
     &     sqrt(vsum/real(nb))*vu
c     ---------------------------------------------------
c     distrib
      do i=1,nb
         rr(i)=sqrt(xx(i)**2.+yy(i)**2.)*ru 
         vr(i)=ww(i)*vu/sqrt(2.)                 
      enddo
c     
      a=0.01*rs
      b=50.*rs
      call ldistrib(a,b,ml,nb,rr,vr,xn,x1,x2,xa,vra,vrs,np)

c
      sr0=2.*pi*rs**2.*rho0* (xa(1)/rc)*log(rc/xa(1))
      do i=1,np
         rqa=1.+xa(i)**2./ra**2.
         sd=(xn(i)+1.e-10)/(pi*(x2(i)**2.-x1(i)**2.) )
         sdt=20000./(1.+xa(i)**2./(rc)**2.)**2.
         print*,xa(i)/rs,vrs(i),vra(i),log10(sd)
      enddo


 100  format(15e16.6)
      end
      subroutine ldistrib(a,b,ml,n,x,y,xn,x1,x2,xa,ya,ys,np)
c     calcula la distrib en el interv LOG x\in (log a,log b) del vector y(n)
c     xn: numero de part en x
c     ya: valor medio de y en x
c     ys: dispersion en torno a ya
      implicit real*8 (a-h,o-z)
      parameter(cero=1.e-10,eps=1.e-3,ismax=50)
      real*8 y(n),x(n),ya(np),ys(np),xa(np),xn(np),ml(n),x1(np),x2(np)
      real*8 tm(n)
      real*8 mls(n),xs(n),vs(n)
      dx=(log10(b+cero)-log10(a+cero))/real(np)
      xx=log10(a+cero)

      do i=1,np         
         xx=xx+dx
c     valor medio
         dk=0.
         ya(i)=0.
         num=0
         do j=1,n
            if (log10(x(j)).le.xx.and.log10(x(j)).gt.xx-dx) then
               num=num+1
               dk=dk+ml(j)
               ya(i)=ya(i)+y(j)*ml(j)
c     guardamos estas particulas para el 3sigma clipping
               mls(num)=ml(j)
               xs(num)=x(j)
               vs(num)=y(j)
            endif            
         enddo
         xa(i)=10.**(xx-dx/2.)
         x1(i)=10.**(xx-dx)
         x2(i)=10.**xx
         ya(i)=ya(i)/real(dk)
         xn(i)=dk

c     dispersion
         ys(i)=0.
         do j=1,n
            if (tm(j).lt.0.001) then
            if (log10(x(j)).le.xx.and.log10(x(j)).gt.xx-dx) then
               ys(i)=ys(i)+(y(j)-ya(i))**2.*ml(j)
            endif
            endif
         enddo
         if (num.gt.1) then
            ys(i)=sqrt( ys(i)/dk )
         else
            ya(i)=-1000.
            ys(i)=-1000.
         endif
      enddo

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
      FUNCTION rtbis(q,func,g,rm,vm,x1,x2,xacc)
      INTEGER JMAX
      REAL*8 rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=60)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      real*8 q,rm,vm,g

      fmid=q-func(x2,g,rm,vm)
      f=q-func(x1,g,rm,vm)
      if(f*fmid.ge.0.) then
         print*, 'root must be bracketed in rtbis'
         print*,x1,x2,f,fmid
         stop
      endif
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=q-func(xmid,g,rm,vm)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      print*, 'too many bisections in rtbis'
      stop
      END

      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
      subroutine SORT2( n,ra,list )
      integer list, ir, l, jb
      real*8  ra
      dimension ra(n), list(n) 

c This is a routine meant to sort the particles in increasing
c order of the radius. Cf. Press et al. 1987, p. 231. 

      l = n/2 + 1
      ir = n

 10   if( l .gt. 1 ) then 
         l = l-1
         rra = ra(l)
         jb = list(l)
      else
         rra = ra(ir)
         jb = list(ir)
         ra(ir) = ra(1)
         list(ir) = list(1)
         ir = ir - 1
         if( ir .eq. 1 ) then 
            ra(1) = rra
            list(1) = jb
            return
         endif
      endif

      i = l
      j = l + l 

 20   if( j .le. ir ) then 
         if( j .lt. ir ) then 
            if( ra(j) .lt. ra(j+1) ) j = j+1
         endif
         if( rra .lt. ra(j) ) then 
            ra(i) = ra(j) 
            list(i) = list(j)
            i = j
            j = j + j
         else
            j = ir + 1
         endif
         go to 20
      endif

      ra(i ) = rra
      list(i) = jb 
      go to 10

      stop
      end

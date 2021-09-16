C***********************************************************************
C
C
        PROGRAM buildgal
C
C
C***********************************************************************
C
C
C   Program to generate initial conditions for simulations with disk-
C   bulge-halo models.  This version allows gas to be included in the
C   disk and enables the user to add a satellte.
C
C   The code assumes the following dimensionless system of units:
C
C      1. G = 1
C      2. h_disk = 1 (exponential scale-length of disk)
C      3. M_disk = 1 (disk mass)
C
C                                                                      
C  VARIABLES :                                                         
C ------------                                                         
C
C     dadrtab  : radial acceleration gradient table.
C     gasmass  : mass of gas in the disk.
C     gastemp  : temperature of the gas.
C     lulog    : logical unit for log output.
C     lustat   : logical unit for statistics output.
C     maxtab   : maximum number of radial bins in tables.
C     nbodies  : total number of particles.
C     ndgas    : number of gas particles in the disk.
C     qsolar   : Toomre Q parameter in the solar neighborhood.             
C     rmingas  : minimum radius for gas particles.
C     rotcirc  : circular equilibrium velocity for each particle.
C     rotmean  : mean circular velocity for each particle.
C     selfggas : option to include gas self-grav.
C     sigphi   : phi velocity dispersion for each particle.
C     sigr     : radial velocity dispersion for each particle.
C     sigt     : critical Toomre dispersion (xQ) for each particle.
C     sigz     : z velocity dispersion for each particle.
C     z0       : disk scale height.
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        INCLUDE 'buildgal.h'

C   Declaration of local variables:
C   -------------------------------

        INTEGER lulog,lustat
        PARAMETER (lulog=8,lustat=18)
        real*4 mlx,mly,mlz,mpx,mpy,mpz
       
        CHARACTER*40 filename,statfile,yesno
        PARAMETER (filename='SPIRAL',statfile='STATS')

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        CALL initsys
c           -------
        CALL initdisk
c           --------
        CALL initbulg
c           --------
        CALL inithalo

C           --------
        CALL diskvel
C           -------
        CALL diskstat(lustat,statfile)
C           --------
        IF(usehalo.AND.selfghal) CALL halovel
c                                -------
        IF(usebulge.AND.selfgbul) CALL bulgevel
c                                 --------

        call nonspherical


        CALL cmtv
c          ----
        CALL addsat
c          ------
        CALL stackmod
c          --------
        CALL outbods(filename)
c           -------

        end





c  ----------------------------------------------------------------
        subroutine nonspherical
c   --------------------------------------------------------------        
        include 'buildgal.h'
        character*40 yesno
        real mpx,mpy,mpz,mlx,mly,mlz
        
c     Builgal is too slow in non-symetric bulge calculation, so we change 
c     the scale length in z (cbulge) and plane x-y (abulge), having an
c     oblate in positions and velocities. The original datas about r,v 
c     come from spherical bulge calculation.

        write(6,*) ' '
        write(6,*) 'do you use orginal datas from buildgal?(y/n)'
        read (5,*) yesno
        if (yesno.eq.'n') goto 5
        if (yesno.eq.'y') then
           open (15,file='datos.dat')
           do i=ndisk+nbulge+1,nbodies
              write(15,200) x(i),y(i),z(i),vx(i),vy(i),vz(i)
           enddo
           close (15) 
           write(6,*) ' '
           write(6,*) 'writing data from original buildgal in datos.dat'
           goto 6
        endif
          

 5      write(6,*) 'do you want a non-spherical halo, man (y/n)?'
        read (5,*) yesno
        if (yesno.eq.'n') goto 10
        if (yesno.eq.'y') then
           write(6,*) 'give me (b/a), bitte'
           read(5,*) epsilon

           e=(1.-epsilon**2.)**0.5

           open(10,file='datos.dat')
           open(14,file='delta.dat')
      
           if (e.ne.0)  fact=asin(e)/e
           if (e.eq.0)  fact=1.
           gamhalo2=gamhalo*fact
           g2=gamhalo2 
c     tomamos b del pot como b=a2, pues si b<<1 tenemos
c     que delta=0 para casi todas las r
           b=gamhalo2*(1.-e**2.)**0.5
         
          
           write(6,*) 'ellepticity:',e,'','  gammahalo',gamhalo,'',
     &      '  gammahalo2',gamhalo2,'  b',b,'  fact',fact

           do i=ndisk+nbulge+1,nbodies
              x(i)=x(i)*fact
              y(i)=y(i)*fact
              z(i)=z(i)*fact*(1.-e**2.)**0.5
              
              
              phipot=ABS(1./(x(i)**2.+y(i)**2.+(abs(z(i))+g2)**2.)**0.5)
              delta=1.-(-g2 + 1./phipot)/(1./phipot**2. - g2**2.)**0.5
              if (delta.gt.1.-epsilon) delta=1.-epsilon

              vx(i)=(3./(3.-delta))**0.5*vx(i)
              vy(i)=(3./(3.-delta))**0.5*vy(i)
              vz(i)=((3.-3.*delta)/(3.-delta))**0.5*vz(i)
              r=(x(i)**2.+y(i)**2.+z(i)**2.)**0.5
              write(14,300) r,delta,phipot
              write(10,200) x(i),y(i),z(i),vx(i),vy(i),vz(i)
              
          enddo

          close(10)
          close(14)
          write(6,*) ' '
          write(6,*) 'halo data w/o copying r,v in datos.dat (ascii)'
       endif 

c     We copy positions and velocities of particles in differents cuadrants
c     to have null-flow of particles (we want no rotation (L=0) and no total 
c     displacement (p=0)). Vease qie la masa es ahora 8*M, luego las vel
c     aumentan en sqrt(8)

       fact=sqrt(8.)

 6     open (11,file='buildjorge.dat')
      
       mlx=0.
       mly=0.
       mlz=0.
       mpx=0.
       mpy=0.
       mpz=0.
       
       do i=ndisk+nbulge+1,nbodies
          write(11,200) x(i),y(i),z(i),
     &            fact*vx(i),fact*vy(i),fact*vz(i)
           mpx=mpx+vx(i)
           mpy=mpy+vy(i)
           mpz=mpz+vz(i)
           mlx=mlx+y(i)*vz(i)-vy(i)*z(i)
           mly=mly+z(i)*vx(i)-x(i)*vz(i)
           mlz=mlz+x(i)*vy(i)-y(i)*vx(i)
            
       
           write(11,200) x(i),y(i),-z(i),
     &            fact*vx(i),fact*vy(i),-fact*vz(i)
           mpx=mpx+vx(i)
           mpy=mpy+vy(i)
           mpz=mpz-vz(i)
           mlx=mlx-y(i)*vz(i)+vy(i)*z(i)
           mly=mly-z(i)*vx(i)+x(i)*vz(i)
           mlz=mlz+x(i)*vy(i)-y(i)*vx(i)
      

           write(11,200) -x(i),-y(i),z(i),
     &            -fact*vx(i),-fact*vy(i),fact*vz(i)
           mpx=mpx-vx(i)
           mpy=mpy-vy(i)
           mpz=mpz+vz(i)
           mlx=mlx-y(i)*vz(i)+vy(i)*z(i)
           mly=mly-z(i)*vx(i)+x(i)*vz(i)
           mlz=mlz+x(i)*vy(i)-y(i)*vx(i)   
   

   
           write(11,200) -x(i),y(i),z(i),
     &            -fact*vx(i),fact*vy(i),fact*vz(i)
           mpx=mpx-vx(i)
           mpy=mpy+vy(i)
           mpz=mpz+vz(i)
           mlx=mlx+y(i)*vz(i)-vy(i)*z(i)
           mly=mly-z(i)*vx(i)+x(i)*vz(i)
           mlz=mlz-x(i)*vy(i)+y(i)*vx(i)
       

       
           write(11,200) x(i),-y(i),z(i),
     &            fact*vx(i),-fact*vy(i),fact*vz(i)
           mpx=mpx+vx(i)
           mpy=mpy-vy(i)
           mpz=mpz+vz(i)
           mlx=mlx-y(i)*vz(i)+vy(i)*z(i)
           mly=mly+z(i)*vx(i)-x(i)*vz(i)
           mlz=mlz-x(i)*vy(i)+y(i)*vx(i)
   

       
           write(11,200) -x(i),-y(i),-z(i),
     &            -fact*vx(i),-fact*vy(i),-fact*vz(i)
           mpx=mpx-vx(i)
           mpy=mpy-vy(i)
           mpz=mpz-vz(i)
           mlx=mlx+y(i)*vz(i)-vy(i)*z(i)
           mly=mly+z(i)*vx(i)-x(i)*vz(i)
           mlz=mlz+x(i)*vy(i)-y(i)*vx(i)
      
       
       
           write(11,200) -x(i),y(i),-z(i),
     &            -fact*vx(i),fact*vy(i),-fact*vz(i)
           mpx=mpx-vx(i)
           mpy=mpy+vy(i)
           mpz=mpz-vz(i)
           mlx=mlx-y(i)*vz(i)+vy(i)*z(i)
           mly=mly+z(i)*vx(i)-x(i)*vz(i)
           mlz=mlz-x(i)*vy(i)+y(i)*vx(i)
       

       
          write(11,200) x(i),-y(i),-z(i),
     &            fact*vx(i),-fact*vy(i),-fact*vz(i)
          mpx=mpx+vx(i)
          mpy=mpy-vy(i)
          mpz=mpz-vz(i)
          mlx=mlx+y(i)*vz(i)-vy(i)*z(i)
          mly=mly-z(i)*vx(i)+x(i)*vz(i)
          mlz=mlz-x(i)*vy(i)+y(i)*vx(i) 

          
       
       enddo

       write(6,*) ' '
       write(6,*) 'px',mpx/nhalo,'  py',mpy/nhalo,'  pz',mpz/nhalo
       write(6,*) 'lx',mlx/nhalo,'  ly',mly/nhalo,'  lz',mlz/nhalo 
       write(6,*) 'total nbodies in buildjorge.dat',nhalo*8
       write(6,*) ' '
       write(6,*) 'copied bulge data in buildjorge.dat'
       close(11)
      

         
              
 200   format(1x,6F10.5)
 300   format(1x,3f10.5)
 10    continue
       
       return

        END


C***********************************************************************
C
C
        SUBROUTINE addsat
C
C
C***********************************************************************
C
C
C     Subroutine to add an optional satellite to the data file.  They
C     are modeled by the density profile analyzed by Hernquist (1990;
C     Ap. J. 356, 359).  Namely,
C
C        rho[S](r) = M_sat * a_sat /(r * (a_sat + r)**3),
C
C     where rho[S](r) is the volume mass density of the satellite at
C     spherical coordinate (r), a_sat is the scale-length of the 
C     satellite and M_sat is the satellite mass.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C=======================================================================

        CALL readsat
C            -------

        IF(usesat) THEN

           CALL insmass
C               -------
           CALL setsat
C               ------
           CALL satvel
C               ------
           CALL cmsat
C               -----
           CALL cmsatmod
C               --------

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE berror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to terminate the program as the result of a fatal
C     error, close the output files, and dump timing information.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message

C=======================================================================

        WRITE(6,*) message

        STOP
        END
C***********************************************************************
C
C
        SUBROUTINE bessel(X,BI0,BI1,BK0,BK1)
C
C
C***********************************************************************
C
C
C   SUBROUTINE TO COMPUTE THE VALUE OF THE MODIFIED BESSEL FUNCTIONS 
C   I0,I1,K0 AND K1 AT X.
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        REAL*8 t,sqx,ct1,ct2,t2,hx,hx2,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,
     &         b6,b7,b8,b9,c1,c2,c3,c4,c5,c6,d1,d2,d3,d4,d5,d6,d7,d8,d9,
     &         e1,e2,e3,e4,e5,e6,e7,f1,f2,f3,f4,f5,f6,f7,g1,g2,g3,g4,g5,
     &         g6,h1,h2,h3,h4,h5,h6,h7,X,BI0,BI1,BK0,BK1,ti,hxi

C=======================================================================
           
        T=X/3.75
        SQX=SQRT(X)
        CT1=SQX*EXP(-X)
        CT2=SQX*EXP(X)
        T2=T*T
        HX=X/2.
        HX2=HX*HX

        A1=3.5156229
        A2=3.0899424
        A3=1.2067492
        A4=0.2659732
        A5=0.0360768
        A6=0.0045813

        B1=0.39894228
        B2=0.01328592
        B3=0.00225319
        B4=-0.00157565
        B5=0.00916281
        B6=-0.02057706
        B7=0.02635537
        B8=-0.01647633
        B9=0.00392377

        C1=0.87890594
        C2=0.51498869
        C3=0.15084934
        C4=0.02658733
        C5=0.00301532
        C6=0.00032411

        D1=0.39894228
        D2=-0.03988024
        D3=-0.00362018
        D4=0.00163801
        D5=-0.01031555
        D6=0.02282967
        D7=-0.02895312
        D8=0.01787654
        D9=-0.00420059

        E1=-0.57721566
        E2=0.42278420
        E3=0.23069756
        E4=0.03488590
        E5=0.00262698
        E6=0.00010750
        E7=0.00000740

        F1=1.25331414
        F2=-0.07832358
        F3=0.02189568
        F4=-0.01062446
        F5=0.00587872
        F6=-0.00251540
        F7=0.00053208

        G1=0.15443144
        G2=-0.67278579
        G3=-0.18156897
        G4=-0.01919402
        G5=-0.00110404
        G6=-0.00004686

        H1=1.25331414
        H2=0.23498619
        H3=-0.03655620
        H4=0.01504268
        H5=-0.00780353
        H6=0.00325614
        H7=-0.00068245

        IF(X.LE.3.75)THEN
           BI0=1.+T2*(A1+T2*(A2+T2*(A3+T2*(A4+T2*(A5+T2*A6)))))
           BI1=0.5+T2*(C1+T2*(C2+T2*(C3+T2*(C4+T2*(C5+T2*C6)))))
           BI1=BI1*X
        ELSE
           TI=1./T
           BI0=B1+TI*(B2+TI*(B3+TI*(B4+TI*(B5+TI*(B6+TI*
     &                      (B7+TI*(B8+TI*B9)))))))
           BI0=BI0/CT1
           BI1=D1+TI*(D2+TI*(D3+TI*(D4+TI*(D5+TI*(D6+TI*
     &                      (D7+TI*(D8+TI*D9)))))))
           BI1=BI1/CT1
        ENDIF

        IF(X.LE.2.)THEN
           BK0=-LOG(HX)*BI0+E1+HX2*
     &                 (E2+HX2*(E3+HX2*(E4+HX2*(E5+HX2*(E6+
     &                  HX2*E7)))))
           BK1=X*LOG(HX)*BI1+1.+HX2*(G1+HX2*(G2+HX2*(G3+HX2*
     &                               (G4+HX2*(G5
     &                                +HX2*G6)))))
           BK1=BK1/X
        ELSE
           HXI=2./X
           BK0=F1+HXI*(F2+HXI*(F3+HXI*(F4+HXI*(F5+HXI*(F6+HXI*F7)))))
           BK0=BK0/CT2
           BK1=H1+HXI*(H2+HXI*(H3+HXI*(H4+HXI*(H5+HXI*(H6+HXI*H7)))))
           BK1=BK1/CT2
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE bulgepot
C
C
C***********************************************************************
C
C
C  Subroutine to compute potentials of bulge particles.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ninterp,j,smindex,i,ihalo
        REAL*8 deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,drdeldrg,
     &         drsm,drdotdr,phsmooth,one,phsm,drhalo

        DIMENSION phsmooth(0:30001)

C=======================================================================

        ninterp=30000

        deldrg=2./ninterp
        one=1.0

        DO 5 i=0,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.1.0) THEN
              phsmooth(i)=-2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
           ELSE
              phsmooth(i)=-one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-
     &                    xw3/30.)
           ENDIF
           IF(xw.GE.2.0) THEN
              phsmooth(i)=one
           ENDIF
 5      CONTINUE

        DO 30 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
           pot(i)=0.0
 30     CONTINUE

C   Compute bulge-disk interaction.
C   ------------------------------

        DO 50 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
c----    CDIR$ IVDEP
           DO 40 J=1,ndisk
              AA=X(I)-X(J)
              BB=Y(I)-Y(J)
              CC=Z(I)-Z(J)
              drdotdr=aa**2+bb**2+cc**2
              sdrdotdr=SQRT(drdotdr)
              rinveff=1./(sdrdotdr+1.e-10)
              drdeldrg=sdrdotdr*ninterp/(epsdisk+epsdisk)
              smindex=drdeldrg
              IF(ninterp.LT.smindex) smindex=ninterp
              IF(1.0.LT.drdeldrg-smindex) THEN
                 drsm=1.0
              ELSE
                 drsm=drdeldrg-smindex
              ENDIF
              phsm=(1.-drsm)*phsmooth(smindex)+drsm*phsmooth(1+smindex)
              rinveff=phsm*rinveff
              pot(i)=pot(i)-pmass(J)*rinveff
 40        CONTINUE

           IF(usebulge.AND.selfgbul.AND.axibulge) THEN
c--     CDIR$ IVDEP
              DO 45 J=ndisk+1,I-1
                 AA=X(I)-X(J)
                 BB=Y(I)-Y(J)
                 CC=Z(I)-Z(J)
                 drdotdr=aa**2+bb**2+cc**2
                 sdrdotdr=SQRT(drdotdr)
                 rinveff=1./(sdrdotdr+1.e-10)
                 drdeldrg=sdrdotdr*ninterp/(epsbulge+epsbulge)
                 smindex=drdeldrg
                 IF(ninterp.LT.smindex) smindex=ninterp
                 IF(1.0.LT.drdeldrg-smindex) THEN
                    drsm=1.0
                 ELSE
                    drsm=drdeldrg-smindex
                 ENDIF
                 phsm=(1.-drsm)*phsmooth(smindex)+drsm*
     &                phsmooth(1+smindex)
                 rinveff=phsm*rinveff
                 pot(i)=pot(i)-pmass(J)*rinveff
 45           CONTINUE

c--    CDIR$ IVDEP
              DO 47 J=I+1,ndisk+nbulge
                 AA=X(I)-X(J)
                 BB=Y(I)-Y(J)
                 CC=Z(I)-Z(J)
                 drdotdr=aa**2+bb**2+cc**2
                 sdrdotdr=SQRT(drdotdr)
                 rinveff=1./(sdrdotdr+1.e-10)
                 drdeldrg=sdrdotdr*ninterp/(epsbulge+epsbulge)
                 smindex=drdeldrg
                 IF(ninterp.LT.smindex) smindex=ninterp
                 IF(1.0.LT.drdeldrg-smindex) THEN
                    drsm=1.0
                 ELSE
                    drsm=drdeldrg-smindex
                 ENDIF
                 phsm=(1.-drsm)*phsmooth(smindex)+drsm*
     &                phsmooth(1+smindex)
                 rinveff=phsm*rinveff
                 pot(i)=pot(i)-pmass(J)*rinveff
 47           CONTINUE

           ENDIF

 50     CONTINUE

C   Compute bulge-halo interaction.
C   -------------------------------

        IF(halotype.EQ.'IS') THEN

           drhalo=rhalo(5)-rhalo(4)

           DO 60 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
              ihalo=radsph(i)/drhalo
              ihalo=ihalo+1

              IF(radsph(i).GT.rhalo(ntabhalo)) THEN
                 pot(i)=pot(i)+uhalo(ntabhalo)*rhalo(ntabhalo)/radsph(i)
              ELSE
                 IF(radsph(i).LE.rhalo(2)) THEN
                    pot(i)=pot(i)+uhalo(2)
                 ELSE
                    pot(i)=pot(i)+((radsph(i)-rhalo(ihalo))*
     &                     uhalo(ihalo+1)/drhalo-(radsph(i)-
     &                     rhalo(ihalo+1))*uhalo(ihalo)/drhalo)
                 ENDIF
              ENDIF

 60        CONTINUE

        ELSE

           DO 70 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
              pot(i)= pot(i)-halomass/(radsph(i)+ahalo)
 70        CONTINUE

        ENDIF

C   Compute bulge self-interaction.
C   -------------------------------

        IF(usebulge.AND.selfgbul.AND.(.NOT.axibulge)) THEN

           DO 80 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
              pot(i)= pot(i)-bulgmass/(radsph(i)+abulge)
 80        CONTINUE

        ENDIF

        WRITE(6,*) 'bulgepot>> Bulge potentials computed <<'

        RETURN
        END



C***********************************************************************
C
C
        SUBROUTINE bulgevel
C
C
C***********************************************************************
C
C
C   Subroutine to initialize velocities of bulge particles, using the
C   spherical Jeans equations.  That is, the radial velocity dispersion
C   at radius r is:
C
C                            infinity
C                               /
C                         1     |  G M(r)
C      <v_r ^2 (r)> =  -------  |  ------ rho_h(r) dr
C                      rho_h(r) /   r^2
C                               r
C
C
C   where rho_h(r) is the bulge density at r.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER maxnbin
        PARAMETER(maxnbin=10000)

        INTEGER nbin,i,j,ir,ntest,irlower,irupper,irindex,nturn,ntot,
     &          nxsimp
        REAL*8 xintbin(0:maxnbin+1000),xmbin(0:maxnbin+1000),one,pi,
     &         radmax,zero,ria,ri,rj,r,dr,rhoi,rhoj,rja,twopi,xv,gspeed,
     &         p,cth,sth,signs,phi,gridmass(0:maxnbin+1000),epsilon,
     &         rlower,rupper,tmass,fracmass,vescape,sigfix,erfcc,
     &         vescsig,vphib,radcylb,rgauss,zb,sigzb,sigpb,sigrcylb,
     &         vzb,vrcylb,vbtot,axbulge,cxbulge,bxmass
        REAL*8 dxrand

C=======================================================================

        zero=0.0
        one=1.0
        pi=4.0*ATAN(one)

        nbin=1000

        IF(nbin.GE.maxnbin) CALL berror(' nbin error in bulgevel ')
C                              ------

        IF(.NOT.axibulge) THEN
           radmax=1.01*rmaxbulg
        ELSE
           radmax=1.01*SQRT(rmaxbulg**2+zmaxbulg**2)
        ENDIF

        dr=radmax/nbin

        DO 10 i=0,nbin
           xintbin(i)=0.0
           xmbin(i)=0.0
           gridmass(i)=0.0
 10     CONTINUE

        DO 30 i=1,nbodies
           r=SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
           epsilon=epsdisk
           IF(i.GT.ndisk.AND.i.LT.nbodies-nhalo+1) epsilon=epsbulge
           IF(i.GT.nbodies-nhalo) epsilon=epshalo

           IF(axibulge.AND.i.GT.ndisk.AND.i.LT.nbodies-nhalo+1) GO TO 30

           rlower=r-epsilon/2.0
           rupper=r+epsilon/2.0
           ir=r/dr
           ir=ir+1
           irlower=rlower/dr
           irlower=irlower+1

           IF(rlower.LT.0.0) THEN
              irlower=ABS(rlower)/dr
              irlower= -irlower
           ENDIF

           irupper=rupper/dr
           irupper=irupper+1
           tmass=0.0

           DO 20 j=irlower,irupper
              irindex=j
              IF(j.LT.ir) irindex=ir
              IF(j.GT.maxnbin+1000) irindex=maxnbin+1000

              IF(irlower.EQ.irupper) THEN
                 tmass=tmass+pmass(i)
                 gridmass(irindex)=gridmass(irindex)+pmass(i)
              ELSE
                 IF(j.NE.irlower.AND.j.NE.irupper) THEN
                    fracmass=(dr/epsilon)*pmass(i)
                    tmass=tmass+fracmass
                    gridmass(irindex)=gridmass(irindex)+fracmass
                 ELSE
                    IF(j.EQ.irupper) THEN
                       fracmass=((rupper-(irupper-1)*dr)/epsilon)*
     &                          pmass(i)
                       tmass=tmass+fracmass
                       gridmass(irindex)=gridmass(irindex)+fracmass
                    ELSE
                       fracmass=((irlower*dr-rlower)/epsilon)*
     &                          pmass(i)
                       tmass=tmass+fracmass
                       gridmass(irindex)=gridmass(irindex)+fracmass
                    ENDIF
                 ENDIF
              ENDIF

 20        CONTINUE

           IF(ABS(tmass-pmass(i))/pmass(i).GT.1.e-4) THEN
              CALL berror(' mass assignment error in bulgevel ')
           ENDIF

 30     CONTINUE

        DO 100 i=1,nbin
           xmbin(i)=xmbin(i-1)+gridmass(i)
 100    CONTINUE

        xmbin(nbin+1)=xmbin(nbin)

        DO 130 i=1,nbin

           ri=i*dr
           ria=ri/abulge
           rhoi=bulgmass/(2.*pi*abulge**3*ria*(1.+ria)**3)

           DO 125 j=i,nbin
              rj=j*dr+0.5*dr
              rja=rj/abulge
              rhoj=bulgmass/(2.*pi*abulge**3*rja*(1.+rja)**3)
              xintbin(i)=xintbin(i)+rhoj*0.5*(xmbin(j)+
     &                   xmbin(j+1))*dr/(rj*rj)
 125       CONTINUE

           xintbin(i)=xintbin(i)/rhoi

 130    CONTINUE
   
        xintbin(nbin+1)=0.0

        DO 140 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
           r=SQRT(x(i)**2+y(i)**2+z(i)**2)
           ir=r/dr
           ir=MIN(ir,nbin)
           sigr(i)=((xintbin(ir+1)-xintbin(ir))*r+xintbin(ir)*
     &             (ir+1)*dr-xintbin(ir+1)*ir*dr)/dr
           sigr(i)=SQRT(sigr(i))
 140    CONTINUE

C   Initialize velocities isotropically so that the distribution of
C   speeds is proportional to v^2 EXP[-v^2/(2*sigma_r^2)].  Limit
C   speed to the local escape speed.
C   ---------------------------------------------------------------

        CALL bulgepot
C            --------

        ntest=nbodies-nhalo-nbulge+1

        one=1.0
        twopi=2.0*4.0*ATAN(one)

        IF(.NOT.axibulge) THEN

 160       CONTINUE

           vescape=SQRT(2.0*ABS(pot(ntest)))
           xv=dxrand()*vescape

           vescsig=vescape/(SQRT(2.0)*sigr(ntest))
           sigfix=1.0-erfcc(vescsig)-8.0*vescsig*(0.75+0.5*vescsig**2)*
     &            EXP(-vescsig**2)/(3.0*SQRT(pi))
           sigfix=sigr(ntest)/SQRT(sigfix)

           gspeed=0.5*xv*xv*EXP(1.0-xv*xv/(2.0*sigfix**2))/sigfix**2

           p=dxrand()

           IF(p.LE.gspeed) THEN
              cth=2.0*(dxrand()-0.5)
              sth=SQRT(1.0-cth*cth)
              signs=2.0*(dxrand()-0.5)
              cth=signs*cth/ABS(signs)
              phi=twopi*dxrand()
              vx(ntest)=xv*sth*COS(phi)
              vy(ntest)=xv*sth*SIN(phi)
              vz(ntest)=xv*cth
              ntest=ntest+1
           ELSE
              GO TO 160
           ENDIF

           IF(ntest.LE.nbodies-nhalo) GO TO 160

        ELSE

           axbulge=abulge
           cxbulge=cbulge
           bxmass=bulgmass
           nxsimp=nsimpson

           DO 180 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
              radcylb=SQRT(x(i)**2+y(i)**2)
              zb=z(i)

              CALL obsigma(radcylb,zb,sigrcylb,sigzb,sigpb,axbulge,
C                  -------
     &                     cxbulge,bxmass,nxsimp,rmaxbulg,zmaxbulg)

              sigrcylb=SQRT(sigrcylb**2+sigr(i)**2)
              sigpb=SQRT(sigpb**2+sigr(i)**2)
              sigzb=SQRT(sigzb**2+sigr(i)**2)

 175          vrcylb=rgauss(zero,sigrcylb)
              vphib=rgauss(zero,sigpb)
              vzb=rgauss(zero,sigzb)
              vx(i)=(vrcylb*x(i)-vphib*y(i))/radcylb
              vy(i)=(vrcylb*y(i)+vphib*x(i))/radcylb
              vz(i)=vzb
              vbtot=vx(i)**2+vy(i)**2+vz(i)**2
              IF(vbtot.GE.2.0*ABS(pot(i))) GO TO 175
 180       CONTINUE

        ENDIF

        IF(bulgerot) THEN
           ntot=0
           nturn=0.5*nbulge*(2.0*brotfrac-1.0)
           IF(brotfrac.EQ.0.5) nturn=0
           IF(brotfrac.EQ.1.0) nturn=nbulge

           DO 300 i=nbodies-nhalo-nbulge+1,nbodies-nhalo
              radcylb=SQRT(x(i)**2+y(i)**2)
              vphib=(x(i)*vy(i)-y(i)*vx(i))/radcylb
              IF(vphib.LT.0.0.AND.ntot.LT.nturn) THEN
                 ntot=ntot+1
                 vx(i) = -vx(i)
                 vy(i) = -vy(i)
              ENDIF
 300       CONTINUE

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        FUNCTION cfunc(xrotcirc,xrad,xkappa)
C
C
C***********************************************************************
C
C
C   Function defining the ratio of the squares of azimuthal and radial
C   velocity dispersion.                            
C                                                                
C   OLD:                                                           
C         A(x) = A(0) exp(-x) + 0.5                              
C                                                                
C   where A(x)=<theta^2>/<phi^2>  , A(0)=0.5 and x=R/h.
C                                                                
C   NEW:                                                           
C         cfunc(rotcirc,rad,kappa) = 0.25 * kappa^2/omega^2         
C                                                                 
C   where  omega = rotcirc/rad                                  
C                                                                
C   This function ensures a circular velocity ellipsoid in the center 
C   and an asymptotic behaviour like that predicted by the epicyclic 
C   approximation.                               
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        REAL*8 xrotcirc,xrad,xkappa,xomega,cfunc

C=======================================================================

        xomega=xrotcirc/xrad

        IF(xomega.EQ.0.0.OR.xkappa.EQ.0.0) THEN
           cfunc=1.0
        ELSE
           cfunc=0.25*xkappa*xkappa/(xomega*xomega)
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE circv
C
C
C***********************************************************************
C
C
C   Subroutine to compute centripital rotation velocities for each
C   particle.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 smallnum

C=======================================================================

        smallnum=1.e-07

        DO 10 i=1,ndisk
           IF(radcyl(i).GT.smallnum) THEN
              IF(aradcyl(i).LT.0.0) THEN
                 rotcirc(i)=SQRT(ABS(aradcyl(i))*radcyl(i))
              ELSE
                 WRITE(6,*) 'circv>>particle:',i,' had arad=',aradcyl(i)
                 rotcirc(i)=0.
              ENDIF
           ELSE
              WRITE(6,*) 'circv>> particle :',i,' had rad = ',radcyl(i)
              rotcirc(i)=0.
           END IF
 10     CONTINUE

        WRITE(6,*) 'circv>> Centripital velocities computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE cmbulge
C
C
C***********************************************************************
C
C
C   Subroutine to transform bulge coordinates to center of mass.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j
        REAL*8 xsum,ysum,zsum,psum,xcm,ycm,zcm

C=======================================================================

        psum=0.0
        xsum=0.0
        ysum=0.0
        zsum=0.0
        
        DO 10 i=nbodies-nbulge+1,nbodies
           psum=psum+pmass(i)
           xsum=xsum+pmass(i)*x(i)
           ysum=ysum+pmass(i)*y(i)
           zsum=zsum+pmass(i)*z(i)
 10     CONTINUE

        xcm=xsum/psum
        ycm=ysum/psum
        zcm=zsum/psum

        DO 20 j=nbodies-nbulge+1,nbodies
           x(j)=x(j)-xcm
           y(j)=y(j)-ycm
           z(j)=z(j)-zcm
           radcyl(j)=SQRT(x(j)*x(j)+y(j)*y(j))
           radsph(j)=SQRT(radcyl(j)*radcyl(j)+z(j)*z(j))
 20     CONTINUE

        WRITE(6,*) 'cmbulge>> CMBULGE completed <<'
        WRITE(6,*) 'cmbulge>> Corrections : x y z : ',xcm,ycm,zcm

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE cmdisk
C
C
C***********************************************************************
C
C
C   Subroutine to transform disk coordinates to center of mass.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j
        REAL*8 xsum,ysum,zsum,psum,xcm,ycm,zcm

C=======================================================================

        psum=0.0
        xsum=0.0
        ysum=0.0
        zsum=0.0
        
        DO 10 i=1,ndisk
           psum=psum+pmass(i)
           xsum=xsum+pmass(i)*x(i)
           ysum=ysum+pmass(i)*y(i)
           zsum=zsum+pmass(i)*z(i)
 10     CONTINUE

        xcm=xsum/psum
        ycm=ysum/psum
        zcm=zsum/psum

        DO 20 j=1,ndisk
           x(j)=x(j)-xcm
           y(j)=y(j)-ycm
           z(j)=z(j)-zcm
           radcyl(j)=SQRT(x(j)*x(j)+y(j)*y(j))
      	   radsph(j)=SQRT(radcyl(j)*radcyl(j)+z(j)*z(j))
C	   print*,'radcyl',radcyl(j)

   20	CONTINUE

        WRITE(6,*) 'cmdisk>> CMDISK completed <<'
        WRITE(6,*) 'cmdisk>> Corrections : x y z : ',xcm,ycm,zcm

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE cmhalo
C
C
C***********************************************************************
C
C
C   Subroutine to transform halo coordinates to center of mass.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j
        REAL*8 xsum,ysum,zsum,psum,xcm,ycm,zcm

C=======================================================================

        psum=0.0
        xsum=0.0
        ysum=0.0
        zsum=0.0
        
        DO 10 i=nbodies-nhalo+1,nbodies
           psum=psum+pmass(i)
           xsum=xsum+pmass(i)*x(i)
           ysum=ysum+pmass(i)*y(i)
           zsum=zsum+pmass(i)*z(i)
 10     CONTINUE

        xcm=xsum/psum
        ycm=ysum/psum
        zcm=zsum/psum

        DO 20 j=nbodies-nhalo+1,nbodies
           x(j)=x(j)-xcm
           y(j)=y(j)-ycm
           z(j)=z(j)-zcm
           radcyl(j)=SQRT(x(j)*x(j)+y(j)*y(j))
           radsph(j)=SQRT(radcyl(j)*radcyl(j)+z(j)*z(j))
 20     CONTINUE

        WRITE(6,*) 'cmhalo>> CMHALO completed <<'
        WRITE(6,*) 'cmhalo>> Corrections : x y z : ',xcm,ycm,zcm

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE cmsat
C
C
C***********************************************************************
C
C
C   Subroutine to transform satellite coordinates to center of mass
C   and to place satellite at its starting point.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j
        REAL*8 xsum,ysum,zsum,vxsum,vysum,vzsum,psum,xcm,ycm,zcm,vxcm,
     &         vycm,vzcm

C=======================================================================

        psum=0.0
        xsum=0.0
        ysum=0.0
        zsum=0.0
        vxsum=0.0
        vysum=0.0
        vzsum=0.0
        
        DO 10 i=nbodies-nsat+1,nbodies
           psum=psum+pmass(i)
           xsum=xsum+pmass(i)*x(i)
           ysum=ysum+pmass(i)*y(i)
           zsum=zsum+pmass(i)*z(i)
           vxsum=vxsum+pmass(i)*vx(i)
           vysum=vysum+pmass(i)*vy(i)
           vzsum=vzsum+pmass(i)*vz(i)
 10     CONTINUE

        xcm=xsum/psum
        ycm=ysum/psum
        zcm=zsum/psum
        vxcm=vxsum/psum
        vycm=vysum/psum
        vzcm=vzsum/psum

        DO 20 j=nbodies-nsat+1,nbodies
           x(j)=x(j)-xcm+xsat
           y(j)=y(j)-ycm+ysat
           z(j)=z(j)-zcm+zsat
           vx(j)=vx(j)-vxcm+vxsat
           vy(j)=vy(j)-vycm+vysat
           vz(j)=vz(j)-vzcm+vzsat
 20     CONTINUE

        WRITE(6,*) 'cmsat>> CMSAT completed <<'
        WRITE(6,*) 'cmsat>> Corrections : x y z vx vy vz : ',
     &             xcm,ycm,zcm,vxcm,vycm,vzcm

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE cmsatmod
C
C
C***********************************************************************
C
C
C   Subroutine to transform entire model to center-of-mass coordinates,
C   including the halo.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j
        REAL*8 xsum,ysum,zsum,vxsum,vysum,vzsum,psum,vxcm,vycm,vzcm,xcm,
     &         ycm,zcm

C=======================================================================

        psum=0.0

        IF(usehalo.AND.(.NOT.selfghal)) psum=psum+halomass
        IF(usebulge.AND.(.NOT.selfgbul)) psum=psum+bulgmass

        xsum=0.0
        ysum=0.0
        zsum=0.0
        vxsum=0.0
        vysum=0.0
        vzsum=0.0
        
        DO 10 i=1,nbodies
           psum=psum+pmass(i)
           xsum=xsum+pmass(i)*x(i)
           ysum=ysum+pmass(i)*y(i)
           zsum=zsum+pmass(i)*z(i)
           vxsum=vxsum+pmass(i)*vx(i)
           vysum=vysum+pmass(i)*vy(i)
           vzsum=vzsum+pmass(i)*vz(i)
 10     CONTINUE

        xcm=xsum/psum
        ycm=ysum/psum
        zcm=zsum/psum
        vxcm=vxsum/psum
        vycm=vysum/psum
        vzcm=vzsum/psum

        DO 20 j=1,nbodies
           x(j)=x(j)-xcm
           y(j)=y(j)-ycm
           z(j)=z(j)-zcm
           vx(j)=vx(j)-vxcm
           vy(j)=vy(j)-vycm
           vz(j)=vz(j)-vzcm
 20     CONTINUE

        WRITE(6,*) 'cmsatmod>> CMSATMOD completed <<'
        WRITE(6,*) 'cmsatmod>> Corrections : x y z vx vy vz : ',
     &             xcm,ycm,zcm,vxcm,vycm,vzcm

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE cmtv
C
C
C***********************************************************************
C
C
C   Subroutine to center-of-mass transform velocities.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j
        REAL*8 vxsum,vysum,vzsum,psum,vxcm,vycm,vzcm

C=======================================================================

        vxsum=0
        vysum=0
        vzsum=0
        psum=0
        
        DO 10 i=1,nbodies
           vxsum=vxsum+pmass(i)*vx(i)
           vysum=vysum+pmass(i)*vy(i)
           vzsum=vzsum+pmass(i)*vz(i)
           psum=psum+pmass(i)
 10     CONTINUE
	print*,'pmass in cmtv',pmass(10)
        vxcm=vxsum/psum
        vycm=vysum/psum
        vzcm=vzsum/psum

        DO 20 j=1,nbodies
           vx(j)=vx(j)-vxcm
           vy(j)=vy(j)-vycm
           vz(j)=vz(j)-vzcm
 20     CONTINUE

        WRITE(6,*) 'cmtv>> CMTV completed <<'
        WRITE(6,*) 'cmtv>> Corrections : vx vy vz : ',vxcm,vycm,vzcm

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE dadr
C
C
C***********************************************************************
C
C
C   Subroutine to compute the azimuthally averaged gradient of the 
C   radial acceleration and store it as a table.  The table is 
C   contructed by computing the accelerations on a uniform polar mesh 
C   (maxtab bins in radius and nphi bins in angle). The acceleration at
C   each mesh point is computed using the n**2 method with softening 
C   length epsdisk. The accelerations are the azimuthally averaged
C   and differenced to compute a mean acceleration gradient across 
C   each radial bin.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER maxbuff,maxmesh,nphi

        PARAMETER (maxbuff=25000,maxmesh=25000,nphi=100)

        INTEGER k
        REAL*8 xmesh(maxmesh),ymesh(maxmesh),rmesh(maxmesh)
        REAL*8 tabbuff(maxbuff),delr

C=======================================================================

        IF(maxbuff.LT.maxtab*nphi.OR.maxmesh.LT.maxtab*nphi)
     &     CALL berror(' Dimension error in dadr ')
C               ------

        CALL setmesh(nphi,xmesh,ymesh,rmesh)
C            -------
        CALL meshacc(nphi,xmesh,ymesh,rmesh,tabbuff)
C            -------
    
C  tabbuff contains the azimuthally averaged radial accelerations.
C  ---------------------------------------------------------------

        delr= MAX(rmax,rmaxgas)/FLOAT(maxtab)

        DO 30 k=1,maxtab
           IF(k.eq.1)THEN
              dadrtab(k)=tabbuff(1)/delr
           ELSE
              dadrtab(k)=(tabbuff(k)-tabbuff(k-1))/delr 
           END IF
  30    CONTINUE 
        
        WRITE(6,*) 'dadr>> Radial acceleration gradient computed <<'  

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE diskstat(lustat,statfile)
C
C
C***********************************************************************
C
C
C   Subroutine to produce a table of disk characteristics.  The 
C   analytic values for the rotational velocity and epicyclic frequency
C   (softening free) are computed for the particular disk and halo 
C   characteristic and stored in vexact and exactk.       
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER maxbin,lutemp

        PARAMETER (maxbin=1000)
        PARAMETER (lutemp=17)
        
        CHARACTER*40 statfile
        INTEGER lustat,nbin(maxbin),kk,ilower,iupper,i,j,n
        REAL*8 rbin(maxbin),surfbin(maxbin),vcbin(maxbin),
     &         sigrbin(maxbin),sigzbin(maxbin),vexact(maxbin),
     &         exactk(maxbin),sigpbin(maxbin),binmean(maxbin),
     &         sigtbin(maxbin),qtoomre(maxbin),c1bin(maxbin),
     &         c2bin(maxbin),c3bin(maxbin),kapsum,kapbin(maxbin),
     &         drhalo,atemp,delr,r,va2,surfsum,vcsum,sigrsum,sigzsum,
     &         sigpsum,sigtsum,summean,csum1,csum2,csum3,omega,rs,bi0,
     &         bi1,bk0,bk1,vcirc2d,vcirc2h,halokapp,diskkapp,vcirc2b,
     &         bulgkapp

C=======================================================================

        drhalo=rhalo(5)-rhalo(4)

        OPEN(unit=lustat,file=statfile,status='unknown')
        OPEN(unit=lutemp,file='DDATA',status='unknown')

        DO 200 kk=1,2

           IF(kk.EQ.1) THEN
              ilower=1
              iupper=ndgas
              atemp=acorrgas
           ELSE
              ilower=ndgas+1
              iupper=ndisk
              atemp=acorr
           ENDIF

           delr=0.5
           r=0.

C va2 is left from before; this needs to be fixed.
C ------------------------------------------------

           va2=1.
C          va2=vhalo*vhalo

           DO 10 i=1,40
              n=0
              kapsum=0.
              surfsum=0.
              vcsum=0.
              sigrsum=0.
              sigzsum=0.
              sigpsum=0.
              sigtsum=0.
              summean=0.
              csum1=0.
              csum2=0.
              csum3=0.

              DO 20 j=ilower,iupper
                 IF(radcyl(j).gt.r.and.radcyl(j).le.r+delr)THEN
                    n=n+1
                    kapsum=kapsum+kappa(j)
                    surfsum=surfsum+surfd(j)
                    vcsum=vcsum+rotcirc(j)
                    summean=summean+rotmean(j)
                    sigrsum=sigrsum+sigr(j)
                    sigzsum=sigzsum+sigz(j)
                    sigpsum=sigpsum+sigphi(j)
                    sigtsum=sigtsum+sigt(j)
                    omega=rotcirc(j)/radcyl(j)
                    IF(omega.ne.0.and.kappa(j).ne.0.)THEN
                       csum1=csum1+0.25*kappa(j)*kappa(j)/(omega*omega)
                       IF(sigr(j).NE.0.0) THEN
                          csum2=csum2+
     &                       rotcirc(j)*rotcirc(j)/(sigr(j)*sigr(j))+
     &                       1.-(radcyl(j)/h)-(radcyl(j)*radcyl(j)/(
     &                       h*sqrt(radcyl(j)*radcyl(j)+2.*atemp*
     &                       atemp)))
                       ENDIF
                       csum3=csum3+
     &                    1.-(radcyl(j)/h)-(radcyl(j)*radcyl(j)/(
     &                    h*sqrt(radcyl(j)*radcyl(j)+2.*atemp*atemp))) 
                    END IF
                 END IF
  20          CONTINUE

              nbin(i)=n
              rbin(i)=r+delr/2.
              rs=rbin(i)/(2.*h)

              CALL bessel(rs,bi0,bi1,bk0,bk1)
C                  ------
              vcirc2d=(rbin(i)*rbin(i)*G*diskmass/(2.*h**3))*
     &                (bi0*bk0-bi1*bk1)

              vcirc2h=0.0
              halokapp=0.0
              vcirc2b=0.0
              bulgkapp=0.0

              IF(usehalo) THEN
                 IF(halotype.EQ.'IS') THEN
                    vcirc2h=va2*(1.-(gamhalo/rbin(i))*atan(rbin(i)/
     &                      gamhalo))
                    halokapp=(va2/(rbin(i)*rbin(i)))*
     &                       (2.-(gamhalo/rbin(i))*atan(rbin(i)/
     &                       gamhalo)-1./(1.+(rbin(i)/gamhalo)*
     &                       (rbin(i)/gamhalo)))
                 ELSE
                    vcirc2h=(halomass*rbin(i))/(rbin(i)+ahalo)**2
                    halokapp=(halomass*(rbin(i)+3.*ahalo)/(rbin(i)*
     &                       (rbin(i)+ahalo)**3))
                 ENDIF
              ENDIF

              IF(usebulge) THEN
                 vcirc2b=(bulgmass*rbin(i))/(rbin(i)+abulge)**2
                 bulgkapp=(bulgmass*(rbin(i)+3.*abulge)/(rbin(i)*
     &                    (rbin(i)+abulge)**3))
              ENDIF

              vexact(i)=sqrt(vcirc2d+vcirc2h+vcirc2b)
              diskkapp=(G*diskmass/(h**3))*(
     &                 2.*bi0*bk0-bi1*bk1+(rbin(i)/(2.*h))*(
     &                 bi1*bk0-bi0*bk1))
              exactk(i)=sqrt(halokapp+diskkapp+bulgkapp)

              IF(n.ne.0)THEN
                 kapbin(i)=kapsum/float(n)
                 surfbin(i)=surfsum/float(n)
                 vcbin(i)=vcsum/float(n)
                 sigrbin(i)=sigrsum/float(n)
                 sigzbin(i)=sigzsum/float(n)
                 sigpbin(i)=sigpsum/float(n)
                 sigtbin(i)=sigtsum/float(n)
                 IF(qsolar.ne.0.)THEN
                    qtoomre(i)=qsolar*sigrbin(i)/sigtbin(i)
                 ELSE
                    qtoomre(i)=0.
                 END IF
                 binmean(i)=summean/float(n)
                 c1bin(i)=csum1/float(n)
                 c2bin(i)=csum2/float(n)
                 c3bin(i)=csum3/float(n)
              ELSE
                 kapbin(i)=0.
                 surfbin(i)=0.
                 vcbin(i)=0.
                 sigrbin(i)=0.
                 sigzbin(i)=0.
                 sigpbin(i)=0.
                 sigtbin(i)=0.
                 binmean(i)=0.
                 c1bin(i)=0.
                 c2bin(i)=0.
                 c3bin(i)=0.
              END IF

              r=r+delr

  10       CONTINUE 

           WRITE(lustat,32)
  32       FORMAT(1x,'   R [kpc]   N      SIGMA     KAPPA',
     &            1x,'    Kexact',
     &            1x,'   Vcirc    Vexact     Vmean    sigma r ',
     &            1x,' sigma T',
     &            1x,'     Q ',
     &            1x,'   sigma phi',
     &            1x,' sigma z')
           WRITE(lustat,33)
  33       FORMAT(1x,' ')

           DO 30 i=1,40
              WRITE(lustat,31)rbin(i),nbin(i),surfbin(i),kapbin(i),
     &                        exactk(i),vcbin(i),vexact(i),
     &                        binmean(i),sigrbin(i),
     &                        sigtbin(i),qtoomre(i),sigpbin(i),
     &                        sigzbin(i)
  31          FORMAT(1x,f10.5,i6,11f10.5)
  30       CONTINUE

           DO 55 i=1,40
  55          WRITE(lutemp,56)rbin(i),c3bin(i),c1bin(i),c2bin(i)
  56       FORMAT(1x,4f10.5)

 200    CONTINUE

        CLOSE(unit=lutemp)
        CLOSE(unit=lustat)

        RETURN
        END 
C***********************************************************************
C
C
        SUBROUTINE diskvel
C
C
C***********************************************************************
C
C
C   Subroutine to initialize velocities of disk particles.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C=======================================================================

        CALL forced
C            ------
        IF(usebulge) CALL forceb
C                         ------
        IF(usehalo) CALL forceh
C                        ------
        CALL radacc
C            ------
        CALL dadr
C            ----
        CALL getkappa
C            --------
        CALL sigalar
C            -------
        CALL sigmar
C            ------
        CALL sigmaz
C            ------
        CALL circv
C            -----
        CALL sigcheck
C            --------
        CALL sigmap
C            ------
        CALL setsigma
C            --------
        CALL meanrot
C            -------
        CALL setrot
C            ------

        RETURN
        END
C***********************************************************************
C
C
                          FUNCTION erfcc(x)
C
C
C***********************************************************************
C
C
C      Function to compute the complementary error function erfc(x),
C      with fractional error everywhere less than 1.2x10^-7.  This
C      algorithm uses a Chebyshev fitting. (cf: Numerical Recipes
C      p.164).
C
C
C=======================================================================
     
C   Declaration of local variables.
C   -------------------------------
     
        REAL*8 erfcc,x,z,t
     
C=======================================================================
     
        z = ABS(x)
        t = 1./(1.+0.5*z)
        erfcc = t*EXP( -z*z -1.26551223 + t*(1.00002368 + t*
     &               ( 0.37409196 + t*(0.09678418 + t*(-0.18628806 + t*
     &               (0.27886807 + t*(-1.13520398 + t*(1.48851587 + t*
     &               (-0.82215223 + t* 0.17087277   )))))))))
     
        IF (x.LT.0.0) erfcc = 2.0 - erfcc
     
        RETURN
        END
C***********************************************************************
C
C
                          FUNCTION expint(x)
C
C
C***********************************************************************
C
C
C      Function to compute the exponential integral function E1(x),
C      with fractional error everywhere less than 2x10^-7.  This
C      algorithm is adapted from Abramowitz and Stegun, sec 5.1.
C
C
C=======================================================================
     
C   Declaration of local variables.
C   -------------------------------
     
        REAL*8 expint,x
        REAL*8 aux1,arg,aux2,aux3
     
C=======================================================================
     
        arg=x

        IF(x.LE.0.0) CALL berror(' arg error in expint ')
C                         ------

        IF(x.LE.1.0) THEN
           aux1=-LOG(arg)-0.57721566+0.99999193*arg-0.24991055*arg**2+
     &          0.05519968*arg**3-0.00976004*arg**4+0.00107857*arg**5
           expint=aux1
        ELSE
           aux1=arg**4+8.5733287401*arg**3+18.0590169730*arg**2+
     &          8.6347608925*arg+0.2677737343
           aux2=arg**4+9.5733223454*arg**3+25.6329561486*arg**2+
     &          21.0996530827*arg+3.9584969228
           aux3=aux1/(aux2*EXP(arg)*arg)
           expint=aux3
        ENDIF
     
        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE forceb
C
C
C***********************************************************************
C
C
C  Subroutine to compute forces on disk particles from the bulge.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ninterp,i,j,smindex
        REAL*8 acsmooth,deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,
     &         r3inveff,drdeldrg,drsm,accsm,drdotdr,arad

        DIMENSION acsmooth(0:30001)

C=======================================================================

        ninterp=30000

        deldrg=2./ninterp

        DO 5 i=0,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.1.0) THEN
              acsmooth(i)=xw3*(4./3.-6.*xw2/5.+0.5*xw3)
           ELSE
              acsmooth(i)=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-
     &                    xw4*xw2/6.
           ENDIF
           IF(xw.GE.2.0) THEN
              acsmooth(i)=1.0
           ENDIF
 5      CONTINUE

        IF(usebulge.AND.selfgbul.AND.axibulge) THEN

           DO 10 I=1,ndisk
c--    CDIR$ IVDEP
              DO 20 J=ndisk+1,ndisk+nbulge
                 AA=X(I)-X(J)
                 BB=Y(I)-Y(J)
                 CC=Z(I)-Z(J)
                 drdotdr=aa**2+bb**2+cc**2
                 sdrdotdr=SQRT(drdotdr)
                 rinveff=1./(sdrdotdr+1.e-10)
                 r3inveff=rinveff/(drdotdr+1.e-10)
                 drdeldrg=sdrdotdr*ninterp/(epsdisk+epsbulge)
                 smindex=drdeldrg
                 IF(ninterp.LT.smindex) smindex=ninterp
                 IF(1.0.LT.drdeldrg-smindex) THEN
                    drsm=1.0
                 ELSE
                    drsm=drdeldrg-smindex
                 ENDIF
                 accsm=(1.-drsm)*acsmooth(smindex)+drsm*
     &                 acsmooth(1+smindex)
                 r3inveff=accsm*r3inveff
                 ax(I)=ax(I)-AA*pmass(J)*r3inveff
                 ay(I)=ay(I)-BB*pmass(J)*r3inveff
                 az(I)=az(I)-CC*pmass(J)*r3inveff
 20           CONTINUE
 10        CONTINUE

        ELSE

           DO 30 i=1,ndisk
              arad= -bulgmass/(radsph(i)+abulge)**2
              ax(i)=ax(i)+arad*x(i)/(radsph(i)+1.e-10)
              ay(i)=ay(i)+arad*y(i)/(radsph(i)+1.e-10)
              az(i)=az(i)+arad*z(i)/(radsph(i)+1.e-10)
 30        CONTINUE

        ENDIF

        WRITE(6,*) 'forceb>> Bulge accelerations computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE forcebs
C
C
C***********************************************************************
C
C
C  Subroutine to compute forces on satellite from the bulge.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ninterp,j,smindex,i
        REAL*8 acsmooth,deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,
     &         r3inveff,drdeldrg,drsm,accsm,drdotdr,phsmooth,one,phsm,
     &         arad

        DIMENSION acsmooth(0:30001),phsmooth(0:30001)

C=======================================================================

        ninterp=30000

        deldrg=2./ninterp
        one=1.0

        DO 5 i=0,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.1.0) THEN
              phsmooth(i)=-2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
              acsmooth(i)=xw3*(4./3.-6.*xw2/5.+0.5*xw3)
           ELSE
              phsmooth(i)=-one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-
     &                    xw3/30.)
              acsmooth(i)=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-
     &                    xw4*xw2/6.
           ENDIF
           IF(xw.GE.2.0) THEN
              phsmooth(i)=one
              acsmooth(i)=1.0
           ENDIF
 5      CONTINUE

        IF(usebulge.AND.selfgbul.AND.axibulge) THEN

c--    CDIR$ IVDEP
           DO 20 J=ndisk+1,ndisk+nbulge
              AA=xsat-X(J)
              BB=ysat-Y(J)
              CC=zsat-Z(J)
              drdotdr=aa**2+bb**2+cc**2
              sdrdotdr=SQRT(drdotdr)
              rinveff=1./(sdrdotdr+1.e-10)
              r3inveff=rinveff/(drdotdr+1.e-10)
              drdeldrg=sdrdotdr*ninterp/(epsbulge+epsbulge)
              smindex=drdeldrg
              IF(ninterp.LT.smindex) smindex=ninterp
              IF(1.0.LT.drdeldrg-smindex) THEN
                 drsm=1.0
              ELSE
                 drsm=drdeldrg-smindex
              ENDIF
              phsm=(1.-drsm)*phsmooth(smindex)+drsm*phsmooth(1+smindex)
              accsm=(1.-drsm)*acsmooth(smindex)+drsm*acsmooth(1+smindex)
              rinveff=phsm*rinveff
              r3inveff=accsm*r3inveff
              potsat=potsat-pmass(J)*rinveff
              axsat=axsat-AA*pmass(J)*r3inveff
              aysat=aysat-BB*pmass(J)*r3inveff
              azsat=azsat-CC*pmass(J)*r3inveff
 20        CONTINUE

        ELSE

           radsat=SQRT(xsat**2+ysat**2+zsat**2)

           potsat= potsat-bulgmass/(radsat+abulge)

           arad= -bulgmass/(radsat+abulge)**2
           axsat=axsat+arad*xsat/(radsat+1.e-10)
           aysat=aysat+arad*ysat/(radsat+1.e-10)
           azsat=azsat+arad*zsat/(radsat+1.e-10)

        ENDIF

        WRITE(6,*) 'forcebs>> Bulge-satellite accelerations computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE forced
C
C
C***********************************************************************
C
C
C  Subroutine to compute forces on disk particles from the disk's self-
C  gravity.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ninterp,i,l,j,smindex
        REAL*8 acsmooth,deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,
     &         r3inveff,drdeldrg,drsm,accsm,drdotdr

        DIMENSION acsmooth(0:30001)

C=======================================================================

        ninterp=30000

        deldrg=2./ninterp

        DO 5 i=0,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.1.0) THEN
              acsmooth(i)=xw3*(4./3.-6.*xw2/5.+0.5*xw3)
           ELSE
              acsmooth(i)=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-
     &                    xw4*xw2/6.
           ENDIF
           IF(xw.GE.2.0) THEN
              acsmooth(i)=1.0
           ENDIF
 5      CONTINUE

        DO 30 L=1,ndisk
           ax(L)=0.
           ay(L)=0.
           az(L)=0.
 30     CONTINUE

        DO 10 I=1,ndisk-1
c--    CDIR$ IVDEP
           DO 20 J=I+1,ndisk
              AA=X(I)-X(J)
              BB=Y(I)-Y(J)
              CC=Z(I)-Z(J)
              drdotdr=aa**2+bb**2+cc**2
              sdrdotdr=SQRT(drdotdr)
              rinveff=1./(sdrdotdr+1.e-10)
              r3inveff=rinveff/(drdotdr+1.e-10)
              drdeldrg=sdrdotdr*ninterp/(epsdisk+epsdisk)
              smindex=drdeldrg
              IF(ninterp.LT.smindex) smindex=ninterp
              IF(1.0.LT.drdeldrg-smindex) THEN
                 drsm=1.0
              ELSE
                 drsm=drdeldrg-smindex
              ENDIF
              accsm=(1.-drsm)*acsmooth(smindex)+drsm*acsmooth(1+smindex)
              r3inveff=accsm*r3inveff
              ax(J)=ax(J)+AA*pmass(I)*r3inveff
              ay(J)=ay(J)+BB*pmass(I)*r3inveff
              az(J)=az(J)+CC*pmass(I)*r3inveff
              ax(I)=ax(I)-AA*pmass(J)*r3inveff
              ay(I)=ay(I)-BB*pmass(J)*r3inveff
              az(I)=az(I)-CC*pmass(J)*r3inveff
 20        CONTINUE
 10     CONTINUE

        WRITE(6,*) 'forced>> Disk accelerations computed <<'
        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE forceds
C
C
C***********************************************************************
C
C
C  Subroutine to compute forces on satellite center from disk particles.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ninterp,j,smindex,i
        REAL*8 acsmooth,deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,
     &         r3inveff,drdeldrg,drsm,accsm,drdotdr,phsmooth,one,phsm

        DIMENSION acsmooth(0:30001),phsmooth(0:30001)

C=======================================================================

        ninterp=30000

        deldrg=2./ninterp
        one=1.0

        DO 5 i=0,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.1.0) THEN
              phsmooth(i)=-2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
              acsmooth(i)=xw3*(4./3.-6.*xw2/5.+0.5*xw3)
           ELSE
              phsmooth(i)=-one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-
     &                    xw3/30.)
              acsmooth(i)=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-
     &                    xw4*xw2/6.
           ENDIF
           IF(xw.GE.2.0) THEN
              phsmooth(i)=one
              acsmooth(i)=1.0
           ENDIF
 5      CONTINUE

        potsat=0.0
        axsat=0.0
        aysat=0.0
        azsat=0.0

c--   CDIR$ IVDEP
        DO 20 J=1,ndisk
           AA=xsat-X(J)
           BB=ysat-Y(J)
           CC=zsat-Z(J)
           drdotdr=aa**2+bb**2+cc**2
           sdrdotdr=SQRT(drdotdr)
           rinveff=1./(sdrdotdr+1.e-10)
           r3inveff=rinveff/(drdotdr+1.e-10)
           drdeldrg=sdrdotdr*ninterp/(epsdisk+epsdisk)
           smindex=drdeldrg
           IF(ninterp.LT.smindex) smindex=ninterp
           IF(1.0.LT.drdeldrg-smindex) THEN
              drsm=1.0
           ELSE
              drsm=drdeldrg-smindex
           ENDIF
           phsm=(1.-drsm)*phsmooth(smindex)+drsm*phsmooth(1+smindex)
           accsm=(1.-drsm)*acsmooth(smindex)+drsm*acsmooth(1+smindex)
           rinveff=phsm*rinveff
           r3inveff=accsm*r3inveff
           potsat=potsat-pmass(J)*rinveff
           axsat=axsat-AA*pmass(J)*r3inveff
           aysat=aysat-BB*pmass(J)*r3inveff
           azsat=azsat-CC*pmass(J)*r3inveff
 20     CONTINUE

        WRITE(6,*) 'forceds>> Disk-satellite accelerations computed <<'
        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE forceh
C
C
C***********************************************************************
C
C
C  Subroutine to compute forces on disk particles from the halo.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,ihalo
        REAL*8 xmhr,drhalo,arad

C=======================================================================

        IF(halotype.EQ.'IS') THEN

           drhalo=rhalo(5)-rhalo(4)

           DO 10 i=1,ndisk
              ihalo=radsph(i)/drhalo
              ihalo=ihalo+1

              IF(radsph(i).GT.rhalo(ntabhalo)) THEN
                 xmhr=xmhalo(ntabhalo)
              ELSE
                 IF(radsph(i).LE.rhalo(2)) THEN
                    xmhr=xmhalo(2)*radsph(i)**3/rhalo(2)**3
                 ELSE
                    xmhr=(radsph(i)-rhalo(ihalo))*xmhalo(ihalo+1)/
     &                   drhalo-(radsph(i)-rhalo(ihalo+1))*
     &                   xmhalo(ihalo)/drhalo
                 ENDIF
              ENDIF

              arad=-xmhr/(radsph(i)**2+1.e-10)

              ax(i)=ax(i)+arad*x(i)/(radsph(i)+1.e-10)
              ay(i)=ay(i)+arad*y(i)/(radsph(i)+1.e-10)
              az(i)=az(i)+arad*z(i)/(radsph(i)+1.e-10)

 10        CONTINUE

        ELSE

           DO 20 i=1,ndisk
              arad= -halomass/(radsph(i)+ahalo)**2
              ax(i)=ax(i)+arad*x(i)/(radsph(i)+1.e-10)
              ay(i)=ay(i)+arad*y(i)/(radsph(i)+1.e-10)
              az(i)=az(i)+arad*z(i)/(radsph(i)+1.e-10)
  20       CONTINUE

        ENDIF

        WRITE(6,*) 'forceh>> Halo accelerations computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE forcehs
C
C
C***********************************************************************
C
C
C  Subroutine to compute forces on satellite from the halo.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ihalo
        REAL*8 xmhr,drhalo,arad

C=======================================================================

        radsat=SQRT(xsat**2+ysat**2+zsat**2)

        IF(halotype.EQ.'IS') THEN

           drhalo=rhalo(5)-rhalo(4)

           ihalo=radsat/drhalo
           ihalo=ihalo+1

           IF(radsat.GT.rhalo(ntabhalo)) THEN
              xmhr=xmhalo(ntabhalo)
              potsat=potsat+uhalo(ntabhalo)*rhalo(ntabhalo)/radsat
           ELSE
              IF(radsat.LE.rhalo(2)) THEN
                 xmhr=xmhalo(2)*radsat**3/rhalo(2)**3
                 potsat=potsat+uhalo(2)
              ELSE
                 xmhr=(radsat-rhalo(ihalo))*xmhalo(ihalo+1)/
     &                drhalo-(radsat-rhalo(ihalo+1))*
     &                xmhalo(ihalo)/drhalo
                 potsat=potsat+((radsat-rhalo(ihalo))*uhalo(ihalo+1)/
     &                  drhalo-(radsat-rhalo(ihalo+1))*
     &                  uhalo(ihalo)/drhalo)
              ENDIF
           ENDIF

           arad=-xmhr/(radsat**2+1.e-10)

           axsat=axsat+arad*xsat/(radsat+1.e-10)
           aysat=aysat+arad*ysat/(radsat+1.e-10)
           azsat=azsat+arad*zsat/(radsat+1.e-10)

        ELSE

           potsat= potsat-halomass/(radsat+ahalo)

           arad= -halomass/(radsat+ahalo)**2
           axsat=axsat+arad*xsat/(radsat+1.e-10)
           aysat=aysat+arad*ysat/(radsat+1.e-10)
           azsat=azsat+arad*zsat/(radsat+1.e-10)

        ENDIF

        WRITE(6,*) 'forcehs>> Halo-satellite accelerations computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE getkappa
C
C
C***********************************************************************
C
C
C   Subroutine to compute the epicyclic frequency distribution.  The 
C   epicyclic frequency is given by :
C
C                               3     d PHI         d^2  PHI
C           KAPPA(R) = sqrt [ -----  --------  +   ---------- ]
C                               R     d R           d  R^2
C
C   where PHI is the potential at R.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,iplace
        REAL*8 smallnum,brack,dadr

C=======================================================================

        smallnum=1.e-07

        DO 10 i=1,ndisk
           iplace=(radcyl(i)/MAX(rmax,rmaxgas))*FLOAT(maxtab)+1.
           IF(iplace.EQ.maxtab+1) THEN
              WRITE(6,*) 'getkappa>> iplace=',iplace,' for i=',i,' <<'
              iplace=maxtab
           ENDIF

           dadr=dadrtab(iplace)
C	   print*,'dadrtab',dadrtab(i)
C	   print*,'aradcyl',aradcyl(i)
           IF(radcyl(i).GT.smallnum) THEN
              brack=-3.*(aradcyl(i)/radcyl(i))-dadr
C	      print*,'brack',brack
	      IF(brack.lt.0.)THEN
                 WRITE(6,*) 'getkappa>> particle :',i,' kappa^2=',brack
                 WRITE(6,*) 'getkappa>>                     ',
     &                      '     rad = ',radcyl(i)
                 kappa(i)=0.

              ELSE
                 kappa(i)=SQRT(brack)
              ENDIF
           ELSE
              WRITE(6,*) 'getkappa>> particle:',i,' had rad=',radcyl(i)
              kappa(i)=0.
           ENDIF

  10    CONTINUE       
       
        WRITE(6,*) 'getkappa>> Kappa computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE halopot
C
C
C***********************************************************************
C
C
C  Subroutine to compute potentials of halo particles.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ninterp,j,smindex,i,ihalo
        REAL*8 deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,drdeldrg,
     &         drsm,drdotdr,phsmooth,one,phsm,drhalo

        DIMENSION phsmooth(0:30001)

C=======================================================================

        ninterp=30000

        deldrg=2./ninterp
        one=1.0

        DO 5 i=0,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.1.0) THEN
              phsmooth(i)=-2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
           ELSE
              phsmooth(i)=-one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-
     &                    xw3/30.)
           ENDIF
           IF(xw.GE.2.0) THEN
              phsmooth(i)=one
           ENDIF
 5      CONTINUE

        DO 30 i=nbodies-nhalo+1,nbodies
           pot(i)=0.0
 30     CONTINUE

C   Compute halo-disk interaction.
C   ------------------------------

        DO 50 i=nbodies-nhalo+1,nbodies
c--   CDIR$ IVDEP
           DO 40 J=1,ndisk
              AA=X(I)-X(J)
              BB=Y(I)-Y(J)
              CC=Z(I)-Z(J)
              drdotdr=aa**2+bb**2+cc**2
              sdrdotdr=SQRT(drdotdr)
              rinveff=1./(sdrdotdr+1.e-10)
              drdeldrg=sdrdotdr*ninterp/(epsdisk+epsdisk)
              smindex=drdeldrg
              IF(ninterp.LT.smindex) smindex=ninterp
              IF(1.0.LT.drdeldrg-smindex) THEN
                 drsm=1.0
              ELSE
                 drsm=drdeldrg-smindex
              ENDIF
              phsm=(1.-drsm)*phsmooth(smindex)+drsm*phsmooth(1+smindex)
              rinveff=phsm*rinveff
              pot(i)=pot(i)-pmass(J)*rinveff
 40        CONTINUE

           IF(usebulge.AND.selfgbul.AND.axibulge) THEN
c--   CDIR$ IVDEP
              DO 45 J=ndisk+1,ndisk+nbulge
                 AA=X(I)-X(J)
                 BB=Y(I)-Y(J)
                 CC=Z(I)-Z(J)
                 drdotdr=aa**2+bb**2+cc**2
                 sdrdotdr=SQRT(drdotdr)
                 rinveff=1./(sdrdotdr+1.e-10)
                 drdeldrg=sdrdotdr*ninterp/(epsbulge+epsbulge)
                 smindex=drdeldrg
                 IF(ninterp.LT.smindex) smindex=ninterp
                 IF(1.0.LT.drdeldrg-smindex) THEN
                    drsm=1.0
                 ELSE
                    drsm=drdeldrg-smindex
                 ENDIF
                 phsm=(1.-drsm)*phsmooth(smindex)+drsm*
     &                phsmooth(1+smindex)
                 rinveff=phsm*rinveff
                 pot(i)=pot(i)-pmass(J)*rinveff
 45           CONTINUE

           ENDIF

 50     CONTINUE

C   Compute halo self-interaction.
C   ------------------------------

        IF(halotype.EQ.'IS') THEN

           drhalo=rhalo(5)-rhalo(4)

           DO 60 i=nbodies-nhalo+1,nbodies
              ihalo=radsph(i)/drhalo
              ihalo=ihalo+1

              IF(radsph(i).GT.rhalo(ntabhalo)) THEN
                 pot(i)=pot(i)+uhalo(ntabhalo)*rhalo(ntabhalo)/radsph(i)
              ELSE
                 IF(radsph(i).LE.rhalo(2)) THEN
                    pot(i)=pot(i)+uhalo(2)
                 ELSE
                    pot(i)=pot(i)+((radsph(i)-rhalo(ihalo))*
     &                     uhalo(ihalo+1)/drhalo-(radsph(i)-
     &                     rhalo(ihalo+1))*uhalo(ihalo)/drhalo)
                 ENDIF
              ENDIF

 60        CONTINUE

        ELSE

           DO 70 i=nbodies-nhalo+1,nbodies
              pot(i)= pot(i)-halomass/(radsph(i)+ahalo)
 70        CONTINUE

        ENDIF

C   Compute halo-bulge interaction.
C   ------------------------------

        IF(usebulge.AND.selfgbul.AND.(.NOT.axibulge)) THEN

           DO 80 i=nbodies-nhalo+1,nbodies
              pot(i)= pot(i)-bulgmass/(radsph(i)+abulge)
 80        CONTINUE

        ENDIF

        WRITE(6,*) 'halopot>> Halo potentials computed <<'

        RETURN
        END



C***********************************************************************
C
C
        SUBROUTINE halovel
C
C
C***********************************************************************
C
C
C   Subroutine to initialize velocities of halo particles, using the
C   spherical Jeans equations.  That is, the radial velocity dispersion
C   at radius r is:
C
C                            infinity
C                               /
C                         1     |  G M(r)
C      <v_r ^2 (r)> =  -------  |  ------ rho_h(r) dr
C                      rho_h(r) /   r^2
C                               r
C
C
C   where rho_h(r) is the halo density at r.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER maxnbin
        PARAMETER(maxnbin=10000)

        INTEGER nbin,i,j,ir,ntest,irlower,irupper,irindex
        REAL*8 xintbin(0:maxnbin+1000),xmbin(0:maxnbin+1000),one,pi,
     &         radmax,zero,alphhalo,qhalo,ria,ri,rj,erfcc,r,dr,rhoi,
     &         rhoj,rja,twopi,p,xv,gspeed,cth,sth,signs,phi,tmass,
     &         fracmass,gridmass(0:maxnbin+1000),epsilon,rlower,rupper,
     &         vescape,sigfix,vmaxsig,rkludge,vmax,potrmax,vkludge
        REAL*8 dxrand

C=======================================================================

        zero=0.0
        one=1.0
        pi=4.0*ATAN(one)

        qhalo=gamhalo/rthalo
        alphhalo=1./(1.-SQRT(pi)*qhalo*EXP(qhalo*qhalo)*erfcc(qhalo))

        nbin=1000

        IF(nbin.GE.maxnbin) CALL berror(' nbin error in halovel ')
C                              ------

        radmax=1.01*rmaxhalo

        dr=radmax/nbin

        DO 10 i=0,nbin
           xintbin(i)=0.0
           xmbin(i)=0.0
           gridmass(i)=0.0
 10     CONTINUE

        rkludge=0.75

        DO 30 i=1,nbodies
           r=SQRT(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
           epsilon=epsdisk
           IF(i.GT.ndisk.AND.i.LT.nbodies-nhalo+1) epsilon=epsbulge
           IF(i.GT.nbodies-nhalo) epsilon=epshalo

           IF(i.LE.ndisk) r=r*rkludge

           rlower=r-epsilon/2.0
           rupper=r+epsilon/2.0
           ir=r/dr
           ir=ir+1
           irlower=rlower/dr
           irlower=irlower+1

           IF(rlower.LT.0.0) THEN
              irlower=ABS(rlower)/dr
              irlower= -irlower
           ENDIF

           irupper=rupper/dr
           irupper=irupper+1
           tmass=0.0

           DO 20 j=irlower,irupper
              irindex=j
              IF(j.LT.ir) irindex=ir
              IF(j.GT.maxnbin+1000) irindex=maxnbin+1000

              IF(irlower.EQ.irupper) THEN
                 tmass=tmass+pmass(i)
                 gridmass(irindex)=gridmass(irindex)+pmass(i)
              ELSE
                 IF(j.NE.irlower.AND.j.NE.irupper) THEN
                    fracmass=(dr/epsilon)*pmass(i)
                    tmass=tmass+fracmass
                    gridmass(irindex)=gridmass(irindex)+fracmass
                 ELSE
                    IF(j.EQ.irupper) THEN
                       fracmass=((rupper-(irupper-1)*dr)/epsilon)*
     &                          pmass(i)
                       tmass=tmass+fracmass
                       gridmass(irindex)=gridmass(irindex)+fracmass
                    ELSE
                       fracmass=((irlower*dr-rlower)/epsilon)*
     &                          pmass(i)
                       tmass=tmass+fracmass
                       gridmass(irindex)=gridmass(irindex)+fracmass
                    ENDIF
                 ENDIF
              ENDIF

 20        CONTINUE

           IF(ABS(tmass-pmass(i))/pmass(i).GT.1.e-4) THEN
              CALL berror(' mass assignment error in halovel ')
           ENDIF

 30     CONTINUE

        DO 100 i=1,nbin
           xmbin(i)=xmbin(i-1)+gridmass(i)
 100    CONTINUE

        xmbin(nbin+1)=xmbin(nbin)

        IF(halotype.EQ.'LH') THEN

           DO 130 i=1,nbin

              ri=i*dr
              ria=ri/ahalo
              rhoi=halomass/(2.*pi*ahalo**3*ria*(1.+ria)**3)

              DO 125 j=i,nbin
                 rj=j*dr+0.5*dr
                 rja=rj/ahalo
                 rhoj=halomass/(2.*pi*ahalo**3*rja*(1.+rja)**3)
                 xintbin(i)=xintbin(i)+rhoj*0.5*(xmbin(j)+
     &                      xmbin(j+1))*dr/(rj*rj)
 125          CONTINUE

              xintbin(i)=xintbin(i)/rhoi

 130       CONTINUE
   
        ELSE

           DO 137 i=1,nbin

              ri=i*dr
              rhoi=halomass*alphhalo*EXP(-ri**2/rthalo**2)/
     &             (2.*pi*SQRT(pi)*rthalo*(ri*ri+gamhalo**2))

              DO 135 j=i,nbin
                 rj=j*dr+0.5*dr
                 rhoj=halomass*alphhalo*EXP(-rj**2/rthalo**2)/
     &                (2.*pi*SQRT(pi)*rthalo*(rj*rj+gamhalo**2))
                 xintbin(i)=xintbin(i)+rhoj*0.5*(xmbin(j)+
     &                      xmbin(j+1))*dr/(rj*rj)
 135          CONTINUE

              xintbin(i)=xintbin(i)/rhoi

 137       CONTINUE

        ENDIF

        xintbin(nbin+1)=0.0

        DO 140 i=nbodies-nhalo+1,nbodies
           r=SQRT(x(i)**2+y(i)**2+z(i)**2)
           ir=r/dr
           ir=MIN(ir,nbin)
           sigr(i)=((xintbin(ir+1)-xintbin(ir))*r+xintbin(ir)*
     &             (ir+1)*dr-xintbin(ir+1)*ir*dr)/dr
           sigr(i)=SQRT(sigr(i))
 140    CONTINUE

C   Initialize velocities isotropically so that the distribution of
C   speeds is proportional to v^2 EXP[-v^2/(2*sigma_r^2)].  Limit
C   speed to the local escape speed.
C   ---------------------------------------------------------------

        CALL halopot
C            -------

        ntest=nbodies-nhalo+1

        one=1.0
        twopi=2.0*4.0*ATAN(one)

        vkludge=0.95

 160    CONTINUE

        vescape=SQRT(2.0*ABS(pot(ntest)))
        
        potrmax= -(halomass+bulgmass+diskmass)/rmaxhalo
        vmax=vescape**2+2.0*potrmax

        IF(vmax.GT.0.0) THEN
           vmax=SQRT(vmax)

           xv=dxrand()*vescape

           vmaxsig=vkludge*vmax/(SQRT(2.0)*sigr(ntest))
           sigfix=1.0-erfcc(vmaxsig)-8.0*vmaxsig*(0.75+0.5*vmaxsig**2)*
     &            EXP(-vmaxsig**2)/(3.0*SQRT(pi))
           sigfix=sigr(ntest)/SQRT(sigfix)

           gspeed=xv*xv*EXP(-xv*xv/(2.0*sigfix**2))

           gspeed=gspeed/(2.0*sigfix*sigfix*EXP(-1.0))

           p=dxrand()
        ELSE
           gspeed=0.0
           p=0.0
        ENDIF

        IF(p.LE.gspeed.OR.vmax.LE.0.0) THEN
           cth=2.0*(dxrand()-0.5)
           sth=SQRT(1.0-cth*cth)
           signs=2.0*(dxrand()-0.5)
           cth=signs*cth/ABS(signs)
           phi=twopi*dxrand()
           vx(ntest)=xv*sth*COS(phi)
           vy(ntest)=xv*sth*SIN(phi)
           vz(ntest)=xv*cth

           IF(SQRT(vx(ntest)**2+vy(ntest)**2+vz(ntest)**2).GT.
     &        vkludge*vmax.AND.vmax.GT.0.0) GO TO 160
              
           ntest=ntest+1
        ELSE
           GO TO 160
        ENDIF

        IF(ntest.LE.nbodies) GO TO 160

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE inbmass
C
C
C***********************************************************************
C
C
C     Subroutine to initialize arrays, etc. having to do with bulge
C     masses.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 mtrunc

C=======================================================================

        DO 10 i=nbodies+1-nbulge,nbodies
           pmass(i)=bulgmass/nbulge
 10     CONTINUE

        IF(.NOT.axibulge) THEN
           bulgmass=bulgmass/(rmaxbulg*rmaxbulg/(rmaxbulg+abulge)**2)
        ELSE
           mtrunc=SQRT(rmaxbulg**2/abulge**2+zmaxbulg**2/cbulge**2)
           bulgmass=bulgmass/(mtrunc**2/(1.0+mtrunc)**2)
        ENDIF

        WRITE(6,*) ' inbmass>> corrected bulge mass = ',bulgmass

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE indmass
C
C
C***********************************************************************
C
C
C     Subroutine to initialize arrays, etc. having to do with disk
C     masses.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 dgasmass

C=======================================================================

        xgasmass=gasmass-gasmass*((1.-(1.+rmingas/h)*
     &           EXP(-rmingas/h))/(1.-(1.+rmaxgas/h)*
     &           EXP(-rmaxgas/h)))
        dgasmass=gasmass-xgasmass
        if(dgasmass.LT.0.0) dgasmass=0.0
        gasmass=gasmass-dgasmass

        xgasmass=gasmass

        IF(.NOT.selfggas) gasmass=0.0

        DO 10 i=ndgas+1,nbodies
           pmass(i)=(diskmass-gasmass)/ndstars
 10     CONTINUE
	

	print*,'diskmass',diskmass
	print*,'ndstars',ndstars
	print*,'pmass in indmass',pmass(1)
        DO 20 i=1,ndgas
           pmass(i)=gasmass/ndgas
 20     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE inhmass
C
C
C***********************************************************************
C
C
C     Subroutine to initialize arrays, etc. having to do with halo
C     masses.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i

C=======================================================================

        DO 10 i=nbodies+1-nhalo,nbodies
           pmass(i)=halomass/nhalo
 10     CONTINUE

        IF(halotype.EQ.'LH') THEN

           halomass=halomass/(rmaxhalo*rmaxhalo/(rmaxhalo+ahalo)**2)

        ELSE

           IF(rmaxhalo.GE.rhalo(ntabhalo)) THEN
              halomass=halomass*halomass/xmhalo(ntabhalo)
           ELSE
              DO 20 i=2,ntabhalo
                 IF(rmaxhalo.GT.rhalo(i-1).AND.rmaxhalo.LE.rhalo(i))
     &              halomass=halomass*halomass/xmhalo(i)
 20           CONTINUE
           ENDIF

        ENDIF

        WRITE(6,*) ' inhmass>> corrected halo mass = ',halomass

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE initbulg
C
C
C***********************************************************************
C
C
C     Subroutine to initialize density structure of the bulge.  Bulges
C     are modeled by the density profile analyzed by Hernquist (1990;
C     Ap. J. 356, 359).  Namely,
C
C        rho[B](r) = M_bulge * a_bulge /(r * (a_bulge + r)**3),
C
C     where rho[B](r) is the volume mass density of the bulge at
C     spherical coordinate (r), a_bulge is the scale-length of the 
C     bulge and M_bulge is the bulge mass.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C=======================================================================

        CALL readbulg
C            --------

        IF(usebulge.AND.selfgbul) THEN

           CALL inbmass
C               -------
           CALL setbulge
C               --------
           CALL cmbulge
C               -------
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE initdisk
C
C
C***********************************************************************
C
C
C     Subroutine to initialize density structure of the disk.  Disks are
C     modeled as isothermal sheets (Spitzer 1942, Ap.J. 95, 329):
C
C        rho[D](R,z) = rho[D](0,0) exp(-R/h) [sech(z/z0)]**2 ,
C
C     where rho[D](R,z) is the volume mass density of the disk at
C     cylinderical coordinates (R,z), h is the  exponential scale
C     length of the disk and z0 is the z scale height.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C=======================================================================

        CALL readdisk
C            --------
        CALL indmass
C            -------
        CALL setdisk
C            -------
        CALL cmdisk
C            ------
        CALL surfden
C            -------

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE inithalo
C
C
C***********************************************************************
C
C
C     Subroutine to initialize density structure of the halo.  Halos are
C     modeled either according to the density profile
C
C        rho[H](r)   EXP[-r*r/(r_c*r_c)]
C       --------- = -------------------
C        rho[H](0)   (1 +  [r/gamma]**2)
C 
C     where rho[H](r) is the volume mass density of the halo at
C     spherical radius r, gamma is the halo mass scale length,
C     and r_c is a cut-off radius; or according to the model analyzed
C     by Hernquist (1990; Ap. J. 356; 359).  Namely,
C
C        rho[H](r) = M_halo * a_halo /(r * (a_halo + r)**3),
C
C     where rho[H](r) is the volume mass density of the halo at
C     spherical coordinate (r), a_halo is the scale-length of the
C     halo and M_halo is the halo mass.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C=======================================================================

        CALL readhalo
C            --------

        IF(usehalo.AND.selfghal) THEN

           CALL inhmass
C               -------
           CALL sethalo
C               -------
           CALL cmhalo
C               ------

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE initsys
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the state of the system.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*1 yesno
        INTEGER iseed

C=======================================================================


        WRITE(6,*) ' Welcome to buildgal .... '
        WRITE(6,*) ' '

        WRITE(6,10)
 10     FORMAT(' Input iseed ',$)
        READ(5,*) iseed

C   Initialize random number generator.
C   -----------------------------------
c--        CALL xraninit(iseed)
C            --------

        WRITE(6,20)
 20     FORMAT(' Output softening lengths (y/n) ? ',$)
        READ(5,'(a)') yesno

        IF(yesno.EQ.'y'.OR.yesno.EQ.'Y') THEN
           outpteps=.TRUE.
        ELSE
           outpteps=.FALSE.
        ENDIF

        h=1.0
        G=1.0
        diskmass=1.0

        nbodies=0

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE insmass
C
C
C***********************************************************************
C
C
C     Subroutine to initialize arrays, etc. having to do with satellite
C     masses.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i

C=======================================================================

        DO 10 i=nbodies+1-nsat,nbodies
           pmass(i)=satmass/nsat
 10     CONTINUE

        satmass=satmass/(rmaxsat*rmaxsat/(rmaxsat+asat)**2)

        WRITE(6,*) ' insmass>> Corrected satellte mass = ',satmass

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE meanrot
C
C
C***********************************************************************
C
C
C   Subroutine to compute the mean circular velocity for each particle
C   according to the C function and Q under the van der Kruit and 
C   Searle radial dispersion profile, according to the collisionless 
C   Boltzmann moment equation.
C
C              2       2      2                           2
C         Vmean - Vcirc = < pi > ( 1. - (R/h) - R dln(< pi >) - C(R/h))
C                                                 -----------
C                                                     dR
C
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,nerror,nout,ierror,imax
        REAL*8 errormax,errorsum,rsuma,rsumb,radmax,xrad,delC,error,
     &         term1,term2,sigratio,errmean,rmeana,cfunc,sigmax,t1max,
     &         t2max,square,rmeanb

C=======================================================================

        DO 5 i=1,ndgas
           rotmean(i)=SQRT(rotcirc(i)**2)
 5      CONTINUE

        nerror=0
        errormax=0.
        errorsum=0.
        nout=0
        rsuma=0.
        rsumb=0.
        radmax=0.

        DO 10 i=ndgas+1,ndisk
           xrad=radcyl(i)/h
           sigratio=cfunc(rotcirc(i),radcyl(i),kappa(i))

           IF(sigratio.GT.1.)THEN
              rsuma=rsuma+radcyl(i)
              nerror=nerror+1
              delC=sigratio-1.
              error=delC*100./sigratio
              errorsum=errorsum+error

              IF(error.GE.errormax) THEN
                 errormax=error
                 ierror=i
              END IF

              sigratio=1.
           END IF 

           term1=rotcirc(i)*rotcirc(i)/(sigr(i)*sigr(i))+1-xrad-
     &           radcyl(i)*radcyl(i)/(h*SQRT(radcyl(i)*radcyl(i)+
     &           2.*acorr*acorr))
           term2=1-xrad-radcyl(i)*radcyl(i)/(h*SQRT(radcyl(i)*
     &           radcyl(i)+2.*acorr*acorr))

           IF(term1.LE.sigratio.OR.term2.GE.sigratio) THEN
              nout=nout+1
              rsumb=rsumb+radcyl(i)

              IF(radcyl(i).ge.radmax)THEN
                 radmax=radcyl(i)
                 imax=i
                 sigmax=sigratio
                 t1max=term1
                 t2max=term2
              END IF

              square=0.

           ELSE

              square=rotcirc(i)*rotcirc(i)-sigr(i)*sigr(i)*
     &               (sigratio+xrad-1.+radcyl(i)*radcyl(i)/(h*
     &               SQRT(radcyl(i)*radcyl(i)+2.*acorr*acorr)))

           ENDIF

           rotmean(i)=SQRT(square)

  10    CONTINUE

        IF(nerror.NE.0) errmean=errorsum/FLOAT(nerror)
        IF(nerror.NE.0) rmeana=rsuma/FLOAT(nerror)

        WRITE(6,*) 'meanrot>> ' 
        WRITE(6,*) 'meanrot>> ',nerror,' particles had C gt 1'
        WRITE(6,*) 'meanrot>> mean error  : ',errmean 
        WRITE(6,*) 'meanrot>> mean radius : ',rmeana 
        WRITE(6,*) 'meanrot>> '
        WRITE(6,*) 'meanrot>> ',nout,' particles out of limits on C'

        IF(nout.NE.0) rmeanb=rsumb/FLOAT(nout)

        WRITE(6,*) 'meanrot>> mean radius = ',rmeanb

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE meshacc(nphi,xmesh,ymesh,rmesh,tabbuff)
C
C
C***********************************************************************
C
C
C   Subroutine to compute the acceleration at each mesh point by a 
C   direct sum over the disk particles and by evaluating the halo
C   contribution at each point.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n,i,j,k,nphi,ihalo
        REAL*8 xmesh(1),ymesh(1),rmesh(1),axm,aym,drhalo,sum,delx,dely,
     &         delz,term,term2,aradd,xmhr,aradh,aradb
        REAL*8 tabbuff(1)

C=======================================================================

        n=0
        drhalo=rhalo(5)-rhalo(4)

        DO 10 i=1,maxtab 
           sum=0.

           DO 20 j=1,nphi
              n=n+1
              axm=0.
              aym=0.

              DO 30 k=1,ndisk
                 delx=xmesh(n)-x(k)
                 dely=ymesh(n)-y(k)
                 delz=z(k)
                 term=delx*delx + dely*dely + delz*delz + epsdisk2
                 term2=term*SQRT(term)
                 axm=axm-pmass(k)*delx/term2
                 aym=aym-pmass(k)*dely/term2
 30           CONTINUE

              aradd=(axm*xmesh(n)+aym*ymesh(n))/rmesh(n)

              IF(usebulge.AND.axibulge.AND.selfgbul) THEN

                 axm=0.0
                 aym=0.0

                 DO 35 k=ndisk+1,ndisk+nbulge
                    delx=xmesh(n)-x(k)
                    dely=ymesh(n)-y(k)
                    delz=z(k)
                    term=delx*delx + dely*dely + delz*delz + epsbulge**2
                    term2=term*SQRT(term)
                    axm=axm-pmass(k)*delx/term2
                    aym=aym-pmass(k)*dely/term2
 35              CONTINUE

                 aradb=(axm*xmesh(n)+aym*ymesh(n))/rmesh(n)

                 aradd=aradd+aradb

              ENDIF

              aradh=0.0

              IF(usehalo) THEN

                 IF(halotype.EQ.'IS') THEN
                    ihalo=rmesh(n)/drhalo
                    ihalo=ihalo+1

                    IF(rmesh(n).GT.rhalo(ntabhalo)) THEN
                       xmhr=xmhalo(ntabhalo)
                    ELSE
                       IF(rmesh(n).LE.rhalo(2)) THEN
                          xmhr=xmhalo(2)*rmesh(n)**3/rhalo(2)**3
                       ELSE
                          xmhr=(rmesh(n)-rhalo(ihalo))*xmhalo(ihalo+1)/
     &                         drhalo-(rmesh(n)-rhalo(ihalo+1))*
     &                         xmhalo(ihalo)/drhalo
                       ENDIF
                    ENDIF

                    aradh=-xmhr/(rmesh(n)**2+1.e-10)

                 ELSE
                    aradh= -halomass/(rmesh(n)+ahalo)**2
                 ENDIF

              ENDIF

              sum=sum+aradd+aradh

              IF(usebulge.AND.(.NOT.axibulge)) THEN
                 sum=sum-bulgmass/(rmesh(n)+abulge)**2
              ENDIF

 20        CONTINUE

           tabbuff(i)=sum/FLOAT(nphi)

 10     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE obsigma(xR,xz,sigR,sigz,sigp,abulge,cbulge,bmass,ns,
     &                     rmaxbulg,zmaxbulg)
C
C
C***********************************************************************
C
C
C     Subroutine to compute velocity dispersions for an axisymmetric
C     mass distribution, according to:
C
C                                          psi
C                                           /
C                                   1       |
C       sigma_R^2 = sigma_z^2 = ----------  | rho(psi,R) d psi ,
C                               rho(psi,R)  |
C                                           /
C                                           0
C
C                                             psi
C                                              /
C                                      1       |   d(rho(psi,R))
C       sigma_phi^2 = sigma_R^2 +  ----------  | R ------------- d psi,
C                                  rho(psi,R)  |        dR
C                                              /
C                                              0
C
C
C     where psi = - potential.
C     
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        INTEGER ns,nsimp
        REAL*8 xR,xz,sigR,sigz,pot,psi,sigt,rho,G,zero,one,M,a,c,
     &         R,z,sigp,abulge,cbulge,bmass,rmaxbulg,zmaxbulg,mRz,pi

        COMMON/parcom/G,zero,one,M,a,c
        COMMON/potcom/R,z
        COMMON/simpcom/nsimp

C=======================================================================

        R=xR
        z=xz

        a=abulge
        c=cbulge
        M=bmass

        zero=0.0
        one=1.0
        G=1.0

        nsimp=ns

        CALL oblatep(xR,xz,pot)
C            -------

        psi=-pot

        CALL obsigt(psi,xR,xz,sigt)
C            ------
        CALL obsigp(psi,xR,xz,sigp)
C            ------

        pi=4.0*ATAN(one)

        mRz=SQRT(xR**2/a**2+xz**2/c**2)
        rho=M/(2.*pi*a*a*c*mRz*(1.+mRz)**3)

        sigR=sigt/rho
        sigp=sigR+sigp/rho
        sigR=SQRT(sigR)

        sigz=sigR

        sigp=SQRT(sigp)

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE oblatep(R,z,pot)
C
C
C***********************************************************************
C
C
C     Subroutine to determine the negative potential psi from the
C     integral:
C
C                       infinity
C                          /
C                    G M   |              1
C       psi (R,z) = -----  | du  -------------------- ,
C                     2    |     delta(u) (1+m(u))**2
C                          /
C                          0
C
C
C     where:
C
C                    R^2         z^2
C       m^2(u) =  --------- + ---------
C                  a^2 + u     c^2 + u
C
C
C       delta (u) = (a^2 + u) SQRT(c^2 + u).
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        INTEGER nscoef
        PARAMETER(nscoef=1000)

        INTEGER i,nsimp
        LOGICAL firstc
        REAL*8 a,c,M,R,z,pot,G,zero,one,pot1,simpcoef(nscoef),
     &         xarg(nscoef),a2,c2,xarg6(nscoef),xarg2(nscoef)

        COMMON/parcom/G,zero,one,M,a,c
        COMMON/simpcom/nsimp

        SAVE firstc,simpcoef,xarg,a2,c2,xarg6,xarg2

        DATA firstc/.TRUE./

C=======================================================================

        IF(firstc) THEN
           firstc=.FALSE.

           IF(2*nsimp+1.GT.nscoef) CALL berror(' overflow in oblatep ')
C                                       ------

           simpcoef(1)=1.0

           DO 10 i=2,2*nsimp,2
              simpcoef(i)=4.0
              simpcoef(i+1)=2.0
 10        CONTINUE

           simpcoef(2*nsimp+1)=1.0

           DO 20 i=1,2*nsimp+1
              xarg(i)=(i-1.0)/(2.0*nsimp)
              xarg2(i)=xarg(i)**2
              xarg6(i)=xarg2(i)**3
 20        CONTINUE

           a2=a*a
           c2=c*c

        ENDIF

        pot=0.0
        pot1=0.0

        DO 30 i=1,2*nsimp+1
           pot=pot+simpcoef(i)/((1.0+SQRT(R*R/(a2+xarg(i))+z*z/
     &         (c2+xarg(i))))**2 * (a2+xarg(i)) * SQRT(c2+xarg(i)))
           pot1=pot1+simpcoef(i)*xarg2(i)/((1.0+SQRT(R*R*xarg6(i)/
     &          (xarg6(i)*a2+1)+z*z*xarg6(i)/(xarg6(i)*c2+1)))**2*
     &          (a2*xarg6(i)+1.0)*SQRT(c2*xarg6(i)+1.0))
 30     CONTINUE

        pot = -0.5*G*M*(pot+6.*pot1)/(2.0*nsimp*3.0)

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE obsigt(psi,R,z,sigt)
C
C
C***********************************************************************
C
C
C     Subroutine to determine the integral:
C
C           psi
C            /
C            |
C       I =  | d-psi rho(psi,R).
C            |
C            /
C            0
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        INTEGER nscoef
        PARAMETER(nscoef=1000)

        INTEGER i,nsimp
        LOGICAL firstc
        REAL*8 a,c,M,R,z,psi,G,zero,one,sigt,simpcoef(nscoef),
     &         xarg(nscoef),rhosave,zsave,rho,znew

        COMMON/parcom/G,zero,one,M,a,c
        COMMON/simpcom/nsimp
        COMMON/savecom/rhosave(nscoef),zsave(nscoef)

        SAVE firstc,simpcoef,xarg

        DATA firstc/.TRUE./

C=======================================================================

        IF(firstc) THEN
           firstc=.FALSE.

           IF(2*nsimp+1.GT.nscoef) CALL berror(' overflow in oblatep ')
C                                       ------

           simpcoef(1)=1.0

           DO 10 i=2,2*nsimp,2
              simpcoef(i)=4.0
              simpcoef(i+1)=2.0
 10        CONTINUE

           simpcoef(2*nsimp+1)=1.0

           DO 20 i=1,2*nsimp+1
              xarg(i)=(i-1.0)/(2.0*nsimp)
 20        CONTINUE

        ENDIF

        sigt=0.0

        DO 30 i=1,2*nsimp+1
           CALL rhopsiR(R,psi*xarg(i),rho,znew)
C               -------
           sigt=sigt+simpcoef(i)*rho
           rhosave(i)=rho
           zsave(i)=znew
 30     CONTINUE

        sigt=sigt*psi/(2.0*nsimp*3.0)

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE obsigp(psi,R,z,sigp)
C
C
C***********************************************************************
C
C
C     Subroutine to determine the integral:
C
C           psi
C            /
C            |         partial-rho(psi,R)
C       I =  | d-psi R ------------------ .
C            |             partial-R
C            /
C            0
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        INTEGER nscoef
        PARAMETER(nscoef=1000)

        INTEGER i,nsimp
        LOGICAL firstc
        REAL*8 a,c,M,R,z,psi,G,zero,one,sigp,simpcoef(nscoef),
     &         rhosave,zsave,mRz,drhodR,drhodz,dpotdz,dpotdR

        COMMON/parcom/G,zero,one,M,a,c
        COMMON/simpcom/nsimp
        COMMON/savecom/rhosave(nscoef),zsave(nscoef)

        SAVE firstc,simpcoef

        DATA firstc/.TRUE./

C=======================================================================

        IF(firstc) THEN
           firstc=.FALSE.

           IF(2*nsimp+1.GT.nscoef) CALL berror(' overflow in oblatep ')
C                                       ------

           simpcoef(1)=1.0

           DO 10 i=2,2*nsimp,2
              simpcoef(i)=4.0
              simpcoef(i+1)=2.0
 10        CONTINUE

           simpcoef(2*nsimp+1)=1.0

        ENDIF

        sigp=0.0

        DO 30 i=1,2*nsimp+1
           mRz=SQRT(R**2/a**2+zsave(i)**2/c**2)

           drhodR=-rhosave(i)*(1.+4.*mRz)*R/(mRz**2*a**2*(1.+mRz))
           drhodz=-rhosave(i)*(1.+4.*mRz)*zsave(i)/(mRz**2*c**2*
     &            (1.+mRz))

           CALL obdpotdz(R,zsave(i),dpotdz)
C               --------
           CALL obdpotdR(R,zsave(i),dpotdR)
C               --------

           drhodR=drhodR-drhodz*dpotdR/dpotdz

           sigp=sigp+simpcoef(i)*R*drhodR

 30     CONTINUE

        sigp=sigp*psi/(2.0*nsimp*3.0)

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE rhopsiR(xR,xpsi,rho,znew)
C
C
C***********************************************************************
C
C
C     Subroutine to compute rho (psi,R), given psi = psi(R,z) and 
C     rho = rho(R,z).  From input values of R and psi, the relation
C     psi(R,z) is inverted numerically to provide z = z(psi,R).
C     The density is then computed from R and z.
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        INTEGER nit
        REAL*8 a,c,M,R,z,pot,G,zero,one,xR,xpsi,rho,zinit,zold,znew,fz,
     &         dfdz,mRz,pi,dpotdz,pottest,zlow,zhigh,potlow,pothigh,
     &         ztest,potmax

        COMMON/parcom/G,zero,one,M,a,c
        COMMON/potcom/R,z

C=======================================================================

        IF(xpsi.EQ.0.0) THEN
           rho=0.0
           znew=1.e10
           RETURN
        ENDIF

        CALL oblatep(xR,zero,potmax)
C            -------

        IF(xpsi.EQ.-potmax) THEN
           znew=zero
           RETURN
        ENDIF

        IF(xpsi.GT.-potmax) THEN
           WRITE(6,*) ' range error in rhopsiR '
           WRITE(6,*) ' psimax = ',-potmax
           STOP
        ENDIF

        pi=4.0*ATAN(one)

        zinit=(G*M/xpsi-a)**2-xR**2

        IF(zinit.LE.1.e-4) THEN
           zinit=c
        ENDIF

        zinit=SQRT(zinit)

        zold=zinit
        znew=zinit

        nit=0

 10     CONTINUE

        CALL oblatep(xR,zold,pot)
C            -------

        fz= -pot - xpsi

        CALL obdpotdz(xR,zold,dpotdz)
C            --------

        dfdz = -dpotdz

        znew = zold - fz/dfdz

        IF(ABS(zold-znew)/zold.GT.1.e-6) THEN

           zold=znew

           nit=nit+1
 
           IF(nit.GT.20) GO TO 15

           GO TO 10

        ENDIF

 15     CALL oblatep(xR,znew,pottest)
C            -------

        IF(ABS(-pottest-xpsi)/xpsi.GT.1.e-6) THEN

           zlow=0.0
           zhigh=1000.

           CALL oblatep(xR,zlow,potlow)
C               -------
           CALL oblatep(xR,zhigh,pothigh)
C               -------

           IF((xpsi+potlow)*(xpsi+pothigh).GE.0.0) THEN
              WRITE(6,*) ' limit error in rhopsiR '
              STOP
           ENDIF

           nit=0

 16        CONTINUE

           ztest=0.5*(zlow+zhigh)

           CALL oblatep(xR,ztest,pottest)
C               -------

           IF(ABS(pottest+xpsi)/xpsi.LE.1.e-6) GO TO 20

           IF(xpsi.LE.-pottest) THEN
              zlow=ztest
           ELSE
              zhigh=ztest
           ENDIF

           nit=nit+1

           IF(nit.GT.100) THEN
              WRITE(6,*) ' nit1 error in rhopsiR '
              STOP
           ENDIF

           GO TO 16

 20        znew=ztest

        ENDIF

        mRz=SQRT(xR**2/a**2+znew**2/c**2)

        rho=M/(2.*pi*a*a*c*mRz*(1.+mRz)**3)

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE obdpotdz(R,z,dpotdz)
C
C
C***********************************************************************
C
C
C     Subroutine to compute partial-psi / partial-z.
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        INTEGER nscoef
        PARAMETER(nscoef=1000)

        INTEGER i,nsimp
        LOGICAL firstc
        REAL*8 a,c,M,R,z,dpotdz,G,zero,one,dpotdz1,simpcoef(nscoef),
     &         xarg(nscoef),a2,c2,xarg6(nscoef),xarg2(nscoef),mofu,
     &         xarg3(nscoef),xarg5(nscoef)

        COMMON/parcom/G,zero,one,M,a,c
        COMMON/simpcom/nsimp

        SAVE firstc,simpcoef,xarg,a2,c2,xarg6,xarg2,xarg3,xarg5

        DATA firstc/.TRUE./

C=======================================================================

        IF(firstc) THEN
           firstc=.FALSE.

           IF(2*nsimp+1.GT.nscoef) CALL berror(' overflow in oblatep ')
C                                       ------

           simpcoef(1)=1.0

           DO 10 i=2,2*nsimp,2
              simpcoef(i)=4.0
              simpcoef(i+1)=2.0
 10        CONTINUE

           simpcoef(2*nsimp+1)=1.0

           DO 20 i=1,2*nsimp+1
              xarg(i)=(i-1.0)/(2.0*nsimp)
              xarg2(i)=xarg(i)**2
              xarg3(i)=xarg2(i)*xarg(i)
              xarg5(i)=xarg2(i)*xarg3(i)
              xarg6(i)=xarg2(i)**3
 20        CONTINUE

           a2=a*a
           c2=c*c

        ENDIF

        dpotdz=0.0
        dpotdz1=0.0

        DO 30 i=1,2*nsimp+1
           mofu=SQRT(R*R/(a2+xarg(i))+z*z/(c2+xarg(i)))
           dpotdz=dpotdz+simpcoef(i)*z/(((1.0+mofu)**3 * 
     &            (a2+xarg(i)) * SQRT(c2+xarg(i)))*mofu*(c2+xarg(i)))
           mofu=SQRT(R*R/(xarg6(i)*a2+1.0)+z*z/(xarg6(i)*c2+1.0))
           dpotdz1=dpotdz1+simpcoef(i)*xarg5(i)*z/((1.0+xarg3(i)*
     &             mofu)**3*(a2*xarg6(i)+1.0)*SQRT(c2*xarg6(i)+1.0)*
     &             mofu*(c2*xarg6(i)+1.0))
 30     CONTINUE

        dpotdz = -0.5*G*M*(-2.*dpotdz-12.*dpotdz1)/(2.0*nsimp*3.0)

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE obdpotdR(R,z,dpotdR)
C
C
C***********************************************************************
C
C
C     Subroutine to compute partial-psi / partial-R.
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        INTEGER nscoef
        PARAMETER(nscoef=1000)

        INTEGER i,nsimp
        LOGICAL firstc
        REAL*8 a,c,M,R,z,dpotdR,G,zero,one,dpotdR1,simpcoef(nscoef),
     &         xarg(nscoef),a2,c2,xarg6(nscoef),xarg2(nscoef),mofu,
     &         xarg3(nscoef),xarg5(nscoef)

        COMMON/parcom/G,zero,one,M,a,c
        COMMON/simpcom/nsimp

        SAVE firstc,simpcoef,xarg,a2,c2,xarg6,xarg2,xarg3,xarg5

        DATA firstc/.TRUE./

C=======================================================================

        IF(firstc) THEN
           firstc=.FALSE.

           IF(2*nsimp+1.GT.nscoef) CALL berror(' overflow in oblatep ')
C                                  ------

           simpcoef(1)=1.0

           DO 10 i=2,2*nsimp,2
              simpcoef(i)=4.0
              simpcoef(i+1)=2.0
 10        CONTINUE

           simpcoef(2*nsimp+1)=1.0

           DO 20 i=1,2*nsimp+1
              xarg(i)=(i-1.0)/(2.0*nsimp)
              xarg2(i)=xarg(i)**2
              xarg3(i)=xarg2(i)*xarg(i)
              xarg5(i)=xarg2(i)*xarg3(i)
              xarg6(i)=xarg2(i)**3
 20        CONTINUE

           a2=a*a
           c2=c*c

        ENDIF

        dpotdR=0.0
        dpotdR1=0.0

        DO 30 i=1,2*nsimp+1
           mofu=SQRT(R*R/(a2+xarg(i))+z*z/(c2+xarg(i)))
           dpotdR=dpotdR+simpcoef(i)*R/(((1.0+mofu)**3 * 
     &            (a2+xarg(i)) * SQRT(c2+xarg(i)))*mofu*(a2+xarg(i)))
           mofu=SQRT(R*R/(xarg6(i)*a2+1.0)+z*z/(xarg6(i)*c2+1.0))
           dpotdR1=dpotdR1+simpcoef(i)*xarg5(i)*R/((1.0+xarg3(i)*
     &             mofu)**3*(a2*xarg6(i)+1.0)*SQRT(c2*xarg6(i)+1.0)*
     &             mofu*(a2*xarg6(i)+1.0))
 30     CONTINUE

        dpotdR = -0.5*G*M*(-2.*dpotdR-12.*dpotdR1)/(2.0*nsimp*3.0)

        RETURN
        END





C***********************************************************************
C
C
        SUBROUTINE outbods(filename)
C
C
C***********************************************************************
C
C
C   Output particle realization.
C
c++     This routine modified to output in different formats, including
c++     binary data. (cmb 21.11.1999). 
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER *40 filename
        INTEGER ndim,ns,i
        REAL*8 t,xt,yt,zt,vxt,vyt,vzt,one,pi,costhet1,sinthet1,cosphi1,
     &         sinphi1,costhet2,sinthet2,cosphi2,sinphi2
        character*1   answer, binmode, outformat 
        character*2   nbodyform 
C=======================================================================

c++     Renormalise gas-particle masses - if any. 

        DO 10 i=1,ndgas
           IF(.NOT.selfggas) pmass(i)=xgasmass/ndgas
 10     CONTINUE

c++     Inquire output format - 

      write (*,*) ' '
      write (*,FMT='($,a)') 'Data format - ascii (a) or binary (b) ? '
      read  (*,'(a)') binmode 
      write (*,FMT='($,a)') 
     &      'Data - linear (LH), nbody6 (nb6)  or superbox (sup)? '
      read  (*,'(a)') nbodyform 

      write (*,*) ' '
      write (*,FMT='($,a)') 'Input model name (< 40 characters) '
      read  (*,'(a)') filename 

      if( binmode .eq. 'a' .or. binmode .eq. 'A' ) then 
          call LH_ascii( filename,nbodyform  ) ! Hernquist ascii mode on 
      else   
          call LH_binary( filename,nbodyform  ) ! Hernquist binary mode on 
      endif 

      RETURN
      END



c++     subroutine to output in the original LH format 
c       cmb 22.11.1999. 

      subroutine LH_ascii( filename,binmode ) 

        INCLUDE 'buildgal.h'

        character*2 binmode 
        CHARACTER *40 filename, formstat 
        INTEGER ndim,ns,i
        REAL*8 t,xt,yt,zt,vxt,vyt,vzt,one,pi,costhet1,sinthet1,cosphi1,
     &         sinphi1,costhet2,sinthet2,cosphi2,sinphi2

        ndim=3
        t=0.
        ns=0

        OPEN(10,file=filename,form='formatted')

c++     Start with one-line details of model

        IF(ndgas.EQ.0) THEN 
              IF(.NOT.addmods) THEN
                 WRITE(10,334) nbodies,ndim,t
              ELSE
                 WRITE(10,334) 2*nbodies,ndim,t
              ENDIF              
        ELSE
           IF(.NOT.addmods) THEN
              WRITE(10,334) nbodies,ndim,t
           ELSE
              WRITE(10,334) 2*nbodies,ndim,t
           ENDIF
        ENDIF
c
        WRITE(6,*) ' %%%%%%%%%% '
        IF(.NOT.addmods) THEN
           WRITE(6,*) nbodies,ndim,t,ndgas
        ELSE
           WRITE(6,*) 2*nbodies,ndim,t,2*ndgas
        ENDIF
        WRITE(6,*) ' %%%%%%%%%% '

c++     This mode only possible if .not. nbody6 or superbox data (cmb)

        if( binmode .eq. 'L' .or. binmode .eq. 'l' ) then 

        DO 20 i=1,ndgas
           WRITE(10,444) pmass(i)
 20     CONTINUE

        IF(addmods) THEN
           DO 30 i=1,ndgas
              WRITE(10,444) pmass(i)
 30        CONTINUE
        ENDIF

        endif 

c++     mass of particles, positions : select mode. 

        if( binmode .eq. 'n' .or. binmode .eq. 'N' ) then 

c++     Output data in Nbody6 format : 

        do i = ndgas+1, nbodies 
           write( 10, FMT='(I6,1x,7(1x,f12.6))' ) 
     &     i, pmass(i), x(i),y(i),z(i),vx(i),vy(i),vz(i)
        end do 

        return 
        
        endif 

        if( binmode .eq. 's' .or. binmode .eq. 'S' ) then 

c++     Output data in SUPERBOX  format : 

           call superbox_out( filename,'ascii' )  
        endif 

c++     Output data in LH format : (the old algorithm) 

        if( binmode .eq. 'L' .or. binmode .eq. 'l' ) then 

        DO 40 i=ndgas+1,nbodies
           WRITE(10,444) pmass(i)
 40     CONTINUE

        IF(addmods) THEN
           DO 50 i=ndgas+1,nbodies
              WRITE(10,444) pmass(i)
 50        CONTINUE
        ENDIF

        IF(.NOT.addmods) THEN
           xmod1=0.0
           ymod1=0.0
           zmod1=0.0
           vxmod1=0.0
           vymod1=0.0
           vzmod1=0.0
           thetmod1=0.0
           phimod1=0.0
           costhet1=1.0
           sinthet1=0.0
           cosphi1=1.0
           sinphi1=0.0
        ELSE
           one=1.0
           pi=4.0*ATAN(one)
           thetmod1=2.0*pi*thetmod1/360.0
           phimod1=2.0*pi*phimod1/360.0
           thetmod2=2.0*pi*thetmod2/360.0
           phimod2=2.0*pi*phimod2/360.0
           costhet1=COS(thetmod1)
           sinthet1=SIN(thetmod1)
           cosphi1=COS(phimod1)
           sinphi1=SIN(phimod1)
           costhet2=COS(thetmod2)
           sinthet2=SIN(thetmod2)
           cosphi2=COS(phimod2)
           sinphi2=SIN(phimod2)
        ENDIF

        DO 60 i=1,ndgas
           xt= x(i)*cosphi1+y(i)*sinphi1*costhet1+z(i)*sinphi1*sinthet1
           yt= -x(i)*sinphi1+y(i)*cosphi1*costhet1+z(i)*cosphi1*sinthet1
           zt= -y(i)*sinthet1+z(i)*costhet1
           WRITE(10,444) xt+xmod1,yt+ymod1,zt+zmod1
 60     CONTINUE

        IF(addmods) THEN
           DO 70 i=1,ndgas
              xt= x(i)*cosphi2+y(i)*sinphi2*costhet2+z(i)*sinphi2*
     &            sinthet2
              yt= -x(i)*sinphi2+y(i)*cosphi2*costhet2+z(i)*cosphi2*
     &            sinthet2
              zt= -y(i)*sinthet2+z(i)*costhet2
              WRITE(10,444) xt+xmod2,yt+ymod2,zt+zmod2
 70        CONTINUE
        ENDIF

        DO 80 i=ndgas+1,nbodies
           xt= x(i)*cosphi1+y(i)*sinphi1*costhet1+z(i)*sinphi1*sinthet1
           yt= -x(i)*sinphi1+y(i)*cosphi1*costhet1+z(i)*cosphi1*sinthet1
           zt= -y(i)*sinthet1+z(i)*costhet1
           WRITE(10,444) xt+xmod1,yt+ymod1,zt+zmod1
 80     CONTINUE

        IF(addmods) THEN
           DO 90 i=ndgas+1,nbodies
              xt= x(i)*cosphi2+y(i)*sinphi2*costhet2+z(i)*sinphi2*
     &            sinthet2
              yt= -x(i)*sinphi2+y(i)*cosphi2*costhet2+z(i)*cosphi2*
     &            sinthet2
              zt= -y(i)*sinthet2+z(i)*costhet2
              WRITE(10,444) xt+xmod2,yt+ymod2,zt+zmod2
 90        CONTINUE
        ENDIF

        DO 100 i=1,ndgas
           vxt= vx(i)*cosphi1+vy(i)*sinphi1*costhet1+vz(i)*sinphi1*
     &          sinthet1
           vyt= -vx(i)*sinphi1+vy(i)*cosphi1*costhet1+vz(i)*cosphi1*
     &          sinthet1
           vzt= -vy(i)*sinthet1+vz(i)*costhet1
           WRITE(10,444) vxt+vxmod1,vyt+vymod1,vzt+vzmod1
 100    CONTINUE

        IF(addmods) THEN
           DO 110 i=1,ndgas
              vxt= vx(i)*cosphi2+vy(i)*sinphi2*costhet2+vz(i)*sinphi2*
     &            sinthet2
              vyt= -vx(i)*sinphi2+vy(i)*cosphi2*costhet2+vz(i)*cosphi2*
     &            sinthet2
              vzt= -vy(i)*sinthet2+vz(i)*costhet2
              WRITE(10,444) vxt+vxmod2,vyt+vymod2,vzt+vzmod2
 110       CONTINUE
        ENDIF

        DO 120 i=ndgas+1,nbodies
           vxt= vx(i)*cosphi1+vy(i)*sinphi1*costhet1+vz(i)*sinphi1*
     &          sinthet1
           vyt= -vx(i)*sinphi1+vy(i)*cosphi1*costhet1+vz(i)*cosphi1*
     &          sinthet1
           vzt= -vy(i)*sinthet1+vz(i)*costhet1
           WRITE(10,444) vxt+vxmod1,vyt+vymod1,vzt+vzmod1
 120    CONTINUE

        IF(addmods) THEN
           DO 130 i=ndgas+1,nbodies
              vxt= vx(i)*cosphi2+vy(i)*sinphi2*costhet2+vz(i)*sinphi2*
     &            sinthet2
              vyt= -vx(i)*sinphi2+vy(i)*cosphi2*costhet2+vz(i)*cosphi2*
     &            sinthet2
              vzt= -vy(i)*sinthet2+vz(i)*costhet2
              WRITE(10,444) vxt+vxmod2,vyt+vymod2,vzt+vzmod2
 130       CONTINUE
        ENDIF

        IF(outpteps) THEN

           DO 160 i=ndgas+1,ndisk
              WRITE(10,444) epsdisk
 160       CONTINUE

           DO 170 i=1,nbulge
              WRITE(10,444) epsbulge
 170       CONTINUE

           DO 180 i=1,nhalo
              WRITE(10,444) epshalo
 180       CONTINUE

           DO 190 i=1,nsat
              WRITE(10,444) epssat
 190       CONTINUE

           IF(addmods) THEN

              DO 200 i=ndgas+1,ndisk
                 WRITE(10,444) epsdisk
 200          CONTINUE

              DO 210 i=1,nbulge
                 WRITE(10,444) epsbulge
 210          CONTINUE

              DO 220 i=1,nhalo
                 WRITE(10,444) epshalo
 220          CONTINUE

           ENDIF

        ENDIF

        DO 140 i=1,ndgas
           WRITE(10,444) gastemp
 140    CONTINUE

        IF(addmods) THEN
           DO 150 i=1,ndgas
              WRITE(10,444) gastemp
 150       CONTINUE
        ENDIF

        CLOSE(10)
        return 

        end if  ! condition on LH format .. 

 333    FORMAT(1x,1i6,1i6,1i6,1i6)
 334    FORMAT(1x,1i6,1i6,1f12.6,1f12.6)
 444    FORMAT(1x,10(1pe14.6))

        RETURN
        END


c++     subroutine to output in the original LH format 
c       cmb 22.11.1999. 

      subroutine LH_binary( filename,binmode ) 

        INCLUDE 'buildgal.h'

        character*1 binmode 
        CHARACTER *40 filename, formstat 
        INTEGER ndim,ns,i
        REAL*8 t,xt,yt,zt,vxt,vyt,vzt,one,pi,costhet1,sinthet1,cosphi1,
     &         sinphi1,costhet2,sinthet2,cosphi2,sinphi2

        ndim=3
        t=0.
        ns=0

c++     mass of particles, positions : select mode. 

        if( binmode .eq. 'n' .or. binmode .eq. 'N' ) then 

        OPEN(10,file=filename,form='unformatted')

c++     Output data in Nbody6 format : (ndgas shoudl = 0, always) 

        WRITE(10) nbodies,ndgas,t

        do i = ndgas+1, nbodies 
        write( 10 ) 
     &  i, pmass(i), x(i),y(i),z(i),vx(i),vy(i),vz(i)
        end do 

        close( 10 ) 
        return 
        
        endif 


        if( binmode .eq. 's' .or. binmode .eq. 'S' ) then 

c++     Output data in SUPERBOX  format : 

           call superbox_out( filename, 'binary' )  
        endif 

c++     Output data in LH format : (the old algorithm) 

      if( binmode .eq. 'L' .or. binmode .eq. 'l' ) then 

        OPEN(10,file=filename,form='unformatted')

        IF(ndgas.EQ.0) THEN 
              IF(.NOT.addmods) THEN
                 WRITE(10) nbodies,ndim,t
              ELSE
                 WRITE(10) 2*nbodies,ndim,t
              ENDIF              
        ELSE
           IF(.NOT.addmods) THEN
              WRITE(10) nbodies,ndim,t
           ELSE
              WRITE(10) 2*nbodies,ndim,t
           ENDIF
        ENDIF
c
        WRITE(6,*) ' %%%%%%%%%% '
        IF(.NOT.addmods) THEN
           WRITE(6,*) nbodies,ndim,t,ndgas
        ELSE
           WRITE(6,*) 2*nbodies,ndim,t,2*ndgas
        ENDIF
        WRITE(6,*) ' %%%%%%%%%% '

        DO 20 i=1,ndgas
           WRITE(10) pmass(i)
 20     CONTINUE

        IF(addmods) THEN
           DO 30 i=1,ndgas
              WRITE(10) pmass(i)
 30        CONTINUE
        ENDIF

        DO 40 i=ndgas+1,nbodies
           WRITE(10) pmass(i)
 40     CONTINUE

        IF(addmods) THEN
           DO 50 i=ndgas+1,nbodies
              WRITE(10) pmass(i)
 50        CONTINUE
        ENDIF

        IF(.NOT.addmods) THEN
           xmod1=0.0
           ymod1=0.0
           zmod1=0.0
           vxmod1=0.0
           vymod1=0.0
           vzmod1=0.0
           thetmod1=0.0
           phimod1=0.0
           costhet1=1.0
           sinthet1=0.0
           cosphi1=1.0
           sinphi1=0.0
        ELSE
           one=1.0
           pi=4.0*ATAN(one)
           thetmod1=2.0*pi*thetmod1/360.0
           phimod1=2.0*pi*phimod1/360.0
           thetmod2=2.0*pi*thetmod2/360.0
           phimod2=2.0*pi*phimod2/360.0
           costhet1=COS(thetmod1)
           sinthet1=SIN(thetmod1)
           cosphi1=COS(phimod1)
           sinphi1=SIN(phimod1)
           costhet2=COS(thetmod2)
           sinthet2=SIN(thetmod2)
           cosphi2=COS(phimod2)
           sinphi2=SIN(phimod2)
        ENDIF

        DO 60 i=1,ndgas
           xt= x(i)*cosphi1+y(i)*sinphi1*costhet1+z(i)*sinphi1*sinthet1
           yt= -x(i)*sinphi1+y(i)*cosphi1*costhet1+z(i)*cosphi1*sinthet1
           zt= -y(i)*sinthet1+z(i)*costhet1
           WRITE(10) xt+xmod1,yt+ymod1,zt+zmod1
 60     CONTINUE

        IF(addmods) THEN
           DO 70 i=1,ndgas
              xt= x(i)*cosphi2+y(i)*sinphi2*costhet2+z(i)*sinphi2*
     &            sinthet2
              yt= -x(i)*sinphi2+y(i)*cosphi2*costhet2+z(i)*cosphi2*
     &            sinthet2
              zt= -y(i)*sinthet2+z(i)*costhet2
              WRITE(10) xt+xmod2,yt+ymod2,zt+zmod2
 70        CONTINUE
        ENDIF

        DO 80 i=ndgas+1,nbodies
           xt= x(i)*cosphi1+y(i)*sinphi1*costhet1+z(i)*sinphi1*sinthet1
           yt= -x(i)*sinphi1+y(i)*cosphi1*costhet1+z(i)*cosphi1*sinthet1
           zt= -y(i)*sinthet1+z(i)*costhet1
           WRITE(10) xt+xmod1,yt+ymod1,zt+zmod1
 80     CONTINUE

        IF(addmods) THEN
           DO 90 i=ndgas+1,nbodies
              xt= x(i)*cosphi2+y(i)*sinphi2*costhet2+z(i)*sinphi2*
     &            sinthet2
              yt= -x(i)*sinphi2+y(i)*cosphi2*costhet2+z(i)*cosphi2*
     &            sinthet2
              zt= -y(i)*sinthet2+z(i)*costhet2
              WRITE(10) xt+xmod2,yt+ymod2,zt+zmod2
 90        CONTINUE
        ENDIF

        DO 100 i=1,ndgas
           vxt= vx(i)*cosphi1+vy(i)*sinphi1*costhet1+vz(i)*sinphi1*
     &          sinthet1
           vyt= -vx(i)*sinphi1+vy(i)*cosphi1*costhet1+vz(i)*cosphi1*
     &          sinthet1
           vzt= -vy(i)*sinthet1+vz(i)*costhet1
           WRITE(10) vxt+vxmod1,vyt+vymod1,vzt+vzmod1
 100    CONTINUE

        IF(addmods) THEN
           DO 110 i=1,ndgas
              vxt= vx(i)*cosphi2+vy(i)*sinphi2*costhet2+vz(i)*sinphi2*
     &            sinthet2
              vyt= -vx(i)*sinphi2+vy(i)*cosphi2*costhet2+vz(i)*cosphi2*
     &            sinthet2
              vzt= -vy(i)*sinthet2+vz(i)*costhet2
              WRITE(10) vxt+vxmod2,vyt+vymod2,vzt+vzmod2
 110       CONTINUE
        ENDIF

        DO 120 i=ndgas+1,nbodies
           vxt= vx(i)*cosphi1+vy(i)*sinphi1*costhet1+vz(i)*sinphi1*
     &          sinthet1
           vyt= -vx(i)*sinphi1+vy(i)*cosphi1*costhet1+vz(i)*cosphi1*
     &          sinthet1
           vzt= -vy(i)*sinthet1+vz(i)*costhet1
           WRITE(10) vxt+vxmod1,vyt+vymod1,vzt+vzmod1
 120    CONTINUE

        IF(addmods) THEN
           DO 130 i=ndgas+1,nbodies
              vxt= vx(i)*cosphi2+vy(i)*sinphi2*costhet2+vz(i)*sinphi2*
     &            sinthet2
              vyt= -vx(i)*sinphi2+vy(i)*cosphi2*costhet2+vz(i)*cosphi2*
     &            sinthet2
              vzt= -vy(i)*sinthet2+vz(i)*costhet2
              WRITE(10) vxt+vxmod2,vyt+vymod2,vzt+vzmod2
 130       CONTINUE
        ENDIF

        IF(outpteps) THEN

           DO 160 i=ndgas+1,ndisk
              WRITE(10) epsdisk
 160       CONTINUE

           DO 170 i=1,nbulge
              WRITE(10) epsbulge
 170       CONTINUE

           DO 180 i=1,nhalo
              WRITE(10) epshalo
 180       CONTINUE

           DO 190 i=1,nsat
              WRITE(10) epssat
 190       CONTINUE

           IF(addmods) THEN

              DO 200 i=ndgas+1,ndisk
                 WRITE(10) epsdisk
 200          CONTINUE

              DO 210 i=1,nbulge
                 WRITE(10) epsbulge
 210          CONTINUE

              DO 220 i=1,nhalo
                 WRITE(10) epshalo
 220          CONTINUE

           ENDIF

        ENDIF

        DO 140 i=1,ndgas
           WRITE(10) gastemp
 140    CONTINUE

        IF(addmods) THEN
           DO 150 i=1,ndgas
              WRITE(10) gastemp
 150       CONTINUE
        ENDIF

        CLOSE(10)

        return 
      end if    ! Condition on LH format 
 
        RETURN
        END


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++     Subroutine superbox_out - to produce outputs in 
c       format suitable for superbox, including headers. 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine superbox_out( fname,binmode )

      include 'buildgal.h' 

      character*1 binmode 
      character*40 fname, lfname, backupname 
      real*4      mtot, rlim, mrlim, mmtot, tcr, gp, gs, gi,gm,go 
      real*4      mrpl,mp,cms(6),hstar(6),fh(60), alphhalo 
      real*4      star(6,nbodsmax), pi, dt, kmspcy
      real*4      a,b,c,d, cv, massref, radref, rc, vta,vtp
      real*8      erfcc, qh 
      integer     ntot,numgal,ih(60),ioff,idum,stno,unitl,ijk,model

c     intialize constants and variables : note that G = 2 in model units. 

      pi    = 4.*atan(1.)
      gs     = 2.

      if( binmode .eq. 'a' .or. binmode .eq. 'A' ) then 
         
      ! Ascii output format ... 

      write( 6,* ) ' Superbox_out - can not output ascii .. '
      stop

      endif 
         
c     input of grid-data :

      write( 6,* ) ' ' 
      write (*,FMT='($,a)') 'input radius grid 1 (in m.u.) : '
      read (*,*) gi
      write (*,FMT='($,a)') 'input radius grid 2 (in m.u.) : '
      read (*,*) gm
      if( gm .le. gi ) write( 6,* ) ' Warnung! grid 2 too small '
      write (*,FMT='($,a)') 'input radius grid 3 (in m.u.) : '
      read (*,*) go
      if( go .le. gm ) write( 6,* ) ' Warnung! grid 3 too small '     
      if( go .le. 1. ) write( 6,* ) ' Warnung! grids too small '
      write(6,* ) ' ' 

      write(*,*) 'Conversion factors for scaling - done once only ' 
      write (*,FMT='($,a)') 'Input mass of this galaxy [Msun]   '
      write (*,FMT='($,a)') '(0 for defining crossing time) '
      read  (*,*) mtot

c       Renormalise component masses to this amount; remember discmass=1
c       by default. 

      if( .not. usehalo ) halomass = 0. 
      if( .not. usebulge ) bulgmass = 0. 

c++   Total model mass = 

      mmtot = 1. + halomass + bulgmass 

      if (mtot .eq. 0.0) then
         write (*,FMT='($,a)') 'input crossing time   [Myrs] :  '
         read  (*,*) tcr
      endif 
      write (*,*) '  '        

9090  write (*,FMT='(a,$)') 'Give a unit of length: pc (0) or kpc (1)  '
      read  (*,*) unitl 
      write (*,FMT='(a,$)') 
     &  'In this unit, what is the size of the galaxy? '
      read  (*,*) rlim 

c++     Determine the model limiting radius - 

      mrlim = max( rmaxbulg, rmaxhalo ) 
      if( halotype .eq. 'IS' .or. halotype .eq. 'is' ) then 
         mrlim = max( mrlim, rthalo ) 
      endif 

c     In pc, years and solar masses,   G = 4.89 10^-3  Myrs^-2Pc^3Sm^-1
c     In kpc, Myears and solar masses, G = 4.89 10^-12 Myrs^-2Kpc^3Sm^-1
c                                        = 4.89 10^-24  yrs^-2Kpc^3Sm^-1
      if( unitl .lt. 0.5 ) then 
              gp    = 4.89e-3
      elseif( unitl .gt. 0.5 .and. unitl .lt. 1.5 ) then 
              gp    = 4.89e-12 
      else 

      write(6,* ) ' Choice out of range ' 
      go to 9090 

      endif 

c       Transformation coefficient - in model units for SUPERBOX 
c       we have 
c       Gs = 2 ; Ms = 1 = Rs   (s = superbox)  
c       I write the physical quantities G, M, R and t such that 
c       Gs = a x G  ; Ms = b x M ; Rs = c x R and ts = d x t 
c
c       The physical unit of mass is the solar mass, time is in years
c       and length has been defined above. Then by construction 
c         
c       [ts]^-2 = GsMs/Rs^3 = ab/c^3 [t]^-2 => [ts] = sqrt(c^3/(ab)) [t]
c       [Vs]^2  = GsMs/Rs   = ab/c GM/R   => [Vs] = sqrt(ab/c) [V] 
c
c       Also 
c       [Ms] = ts^-2Rs^3/Gs = d^-2c^3/a  t^-2R/G => [Ms] = d^-2c^3/a [M]
c       and of course d = sqrt(c^3/(ab)). 
c 
c       Find the conversion coefficients below : store in fh(54:60) 


c     Conversion factor from km/sec -> pc/Myear 

      kmspcy   = (pi/3.)

c     if unit of length is kpc, store km/sec -> kpc/Myear instead 
      if( unitl .gt. 0.5 ) kmspcy = kmspcy / 1000. 
 
      gp= gp 
      a = gs / gp 
      b = mmtot / mtot 
      c = mrlim /rlim 
      d = sqrt( c**3/a/b ) 
      cv= c/d * kmspcy 

c++      else 
c
c      mtot = mmtot/b 
c      rlim = mrlim/c 
c      rpl  = rpl/c 
c
c++      endif 

c++   If there is no bulge, no halo but only a disk : set 
c     'ff' to a dynamical time already worked out by the code. 

      if( usebulge ) then 

c++   For a Hernquist profile a particle a r(0) = a reaches the
c     origin in about 

      ff = 2.572 * ( abulge**3 / gs / bulgmass )**(0.5) 

c++   The (smaller) bulge would provide a shorter timescale - 

      elseif( usehalo .and. .not. usebulge ) then 

        if( halotype .eq. 'LH' .or. halotype .eq. 'lh' ) then       
      
          ff = 2.572 * ( ahalo**3 / gs / halomass )**(0.5) 
        else 
c++     For an isothermal sphere: let the central density define 
c       the dynamical time as in BT 1987, p.38-39 : 
         qh = gamhalo/rthalo 
         alphhalo=1./( 1.-SQRT(pi)*qh*EXP(qh*qh)*
     &                 erfcc(qh) )
         ff = halomass*alphhalo/2/pi**1.5/rthalo**3 / qh**2 
         ff = 1./sqrt(2.) * sqrt( 3.*pi /16./gs/ff ) 
        endif 

      endif 

c     calculating missing input-parameter (mass or crossing-time)
c     The total mass is in solar masses, the time in years. 

      if (mtot .eq. 0.0) then 
         mtot = mmtot * (tcr/ff)**2 * a / c**3  
         write(6,*) ' Total mass = ', mtot, ' solar masses '  
      else
         tcr  = ff / d 
         write(6,*) 
     &       ' Crossing time = ', tcr, ' million years ' 
      endif


c     Compute a (temporary) timestep - 

      write(6,*) ' Note: crossing time = ', tcr*d, ' in m.u. '
      write(6,*) ' Tentative time increament = 1/100 the above ' 
      dt   = tcr*d / 100. 

      ! Binary output mode - default superbox format 

      write(6,'( a,$ )' ) 
     &         ' Data save in physical (0) or model (1) units ? ' 
      read( 5,* ) model  

*       If data need be written out in physical units, rescales data arrays: 
c++ 
      if( model .lt. 0.5 ) then 

         mmtot = mtot 
         tcr   = tcr / d 
         dt    = dt / d 
         mrlim = rlim 
         gi  = gi/c
         gm  = gm/c
         go  = go/c 

c!!!!!    we only write BULGE and HALO particles, copied in each 3D 
c         cuadrant-->we copy each paticle in the another 7 quadrants
c         making a parity transf. We can write DISK particles, \
c         but without copying because it has rotation.
c!!!!!

         fact=sqrt(8.)

         do i = 1,ndisk  

            star(1,i) = x(i)/c 
            star(2,i) = y(i)/c 
            star(3,I) = z(i)/c 
            do j = 1,6  
               cms(j) = 0.
            end do 

            star(4,i) = vx(i)/cv*sqrt(2.)  
            star(5,i) = vy(i)/cv*sqrt(2.)  
            star(6,i) = vz(i)/cv*sqrt(2.)  
         end do 
         
         do i=ndisk+1,nbodies
            star(1,i) = x(i)/c 
            star(2,i) = y(i)/c 
            star(3,I) = z(i)/c
            star(4,i) = vx(i)/cv*sqrt(2.) 
            star(5,i) = vy(i)/cv*sqrt(2.) 
            star(6,i) = vz(i)/cv*sqrt(2.)
         enddo

         j=ndisk
         do i=nbodies+1,(2*nbodies-ndisk)
            j=j+1
            star(1,i) = -x(j)/c 
            star(2,i) = y(j)/c 
            star(3,I) = z(j)/c
            star(4,i) = -vx(j)/cv*sqrt(2.)
            star(5,i) = vy(j)/cv*sqrt(2.)
            star(6,i) = vz(j)/cv*sqrt(2.)
         enddo

         j=ndisk
         do i=2*nbodies-ndisk+1,(3*nbodies-2*ndisk)
            j=j+1
            star(1,i) = x(j)/c 
            star(2,i) = -y(j)/c 
            star(3,I) = z(j)/c
            star(4,i) = vx(j)/cv*sqrt(2.)
            star(5,i) = -vy(j)/cv*sqrt(2.)
            star(6,i) = vz(j)/cv*sqrt(2.)
        enddo

         j=ndisk
         do i=3*nbodies-2*ndisk+1,(4*nbodies-3*ndisk)
            j=j+1
            star(1,i) = x(j)/c 
           star(2,i) = y(j)/c 
            star(3,I) = -z(j)/c
            star(4,i) = vx(j)/cv*sqrt(2.)
            star(5,i) = vy(j)/cv*sqrt(2.)
           star(6,i) = -vz(j)/cv*sqrt(2.)
        enddo

         j=ndisk
         do i=4*nbodies-3*ndisk+1,(5*nbodies-4*ndisk)
            j=j+1
            star(1,i) = -x(j)/c 
            star(2,i) = -y(j)/c 
            star(3,I) = z(j)/c 
            star(4,i) = -vx(j)/cv*sqrt(2.)  
            star(5,i) = -vy(j)/cv*sqrt(2.)  
            star(6,i) = vz(j)/cv*sqrt(2.)
         enddo

         j=ndisk
         do i=5*(nbodies)-4*ndisk+1,6*(nbodies)-5*ndisk
            j=j+1
            star(1,i) = -x(j)/c 
            star(2,i) = y(j)/c 
            star(3,I) = -z(j)/c
            star(4,i) = -vx(j)/cv*sqrt(2.)
            star(5,i) = vy(j)/cv*sqrt(2.)
            star(6,i) = -vz(j)/cv*sqrt(2.)
         enddo

         j=ndisk
         do i=6*(nbodies)-5*ndisk+1,7*(nbodies)-6*ndisk
            j=j+1
            star(1,i) = x(i)/c 
            star(2,i) = -y(i)/c 
            star(3,I) = -z(i)/c
            star(4,i) = vx(i)/cv*sqrt(2.)
            star(5,i) = -vy(i)/cv*sqrt(2.)
            star(6,i) = -vz(i)/cv*sqrt(2.)
         enddo

         j=ndisk
         do i=7*(nbodies)-6*ndisk+1,8*(nbodies)-7*ndisk
            j=j+1
            star(1,i) = -x(j)/c 
            star(2,i) = -y(j)/c 
            star(3,I) = -z(j)/c
            star(4,i) = -vx(j)/cv*sqrt(2.)
            star(5,i) = -vy(j)/cv*sqrt(2.)
            star(6,i) = -vz(j)/cv*sqrt(2.)
         enddo

 

       else 
        
          fact=sqrt(8.)
          do i = 1,ndisk  

             star(1,i) = x(i) 
             star(2,i) = y(i) 
             star(3,I) = z(i) 
             do j = 1,6 
                cms(j) = 0.
             end do 

             star(4,i) = vx(i)* sqrt(2.)
             star(5,i) = vy(i)* sqrt(2.)  
             star(6,i) = vz(i)* sqrt(2.)  
          end do 
         
          do i=ndisk+1,nbodies
             star(1,i) = x(i) 
             star(2,i) = y(i) 
             star(3,I) = z(i)
             star(4,i) = vx(i)* sqrt(2.)  
             star(5,i) = vy(i)* sqrt(2.)  
             star(6,i) = vz(i)* sqrt(2.)
          enddo
          
          j=ndisk
          do i=nbodies+1,(2*nbodies-ndisk)
             j=j+1
             star(1,i) = -x(j)
             star(2,i) = y(j) 
             star(3,I) = z(j)
             star(4,i) = -vx(j)* sqrt(2.)  
             star(5,i) = vy(j)* sqrt(2.)  
             star(6,i) = vz(j)* sqrt(2.)
          enddo

          j=ndisk
          do i=2*nbodies-ndisk+1,(3*nbodies-2*ndisk)
             j=j+1
             star(1,i) = x(j) 
             star(2,i) = -y(j) 
             star(3,I) = z(j)
             star(4,i) = vx(j)* sqrt(2.)
             star(5,i) = -vy(j)* sqrt(2.)
             star(6,i) = vz(j)* sqrt(2.)
          enddo
         
          j=ndisk
          do i=3*nbodies-2*ndisk+1,(4*nbodies-3*ndisk)
             j=j+1
             star(1,i) = x(j) 
             star(2,i) = y(j) 
             star(3,I) = -z(j)
             star(4,i) = vx(j)* sqrt(2.)
             star(5,i) = vy(j)* sqrt(2.)
             star(6,i) = -vz(j)* sqrt(2.)
          enddo

          j=ndisk
          do i=4*nbodies-3*ndisk+1,(5*nbodies-4*ndisk)
             j=j+1
             star(1,i) = -x(j) 
            star(2,i) = -y(j) 
             star(3,I) = z(j) 
             star(4,i) = -vx(j)* sqrt(2.)
             star(5,i) = -vy(j)* sqrt(2.)
             star(6,i) = vz(j)* sqrt(2.)
          enddo
          
          j=ndisk
          do i=5*(nbodies)-4*ndisk+1,6*(nbodies)-5*ndisk
             j=j+1
             star(1,i) = -x(j) 
             star(2,i) = y(j) 
             star(3,I) = -z(j)
             star(4,i) = -vx(j)* sqrt(2.)
             star(5,i) = vy(j)* sqrt(2.) 
             star(6,i) = -vz(j)* sqrt(2.)
         enddo
          
          j=ndisk
          do i=6*(nbodies)-5*ndisk+1,7*(nbodies)-6*ndisk
             j=j+1
            star(1,i) = x(j) 
             star(2,i) = -y(j) 
             star(3,I) = -z(j)
             star(4,i) = vx(j)* sqrt(2.)
             star(5,i) = -vy(j)* sqrt(2.)
             star(6,i) = -vz(j)* sqrt(2.)
          enddo

          j=ndisk
          do i=7*(nbodies)-6*ndisk+1,8*(nbodies)-7*ndisk
             j=j+1
            star(1,i) = -x(j) 
             star(2,i) = -y(j) 
            star(3,I) = -z(j)
             star(4,i) = -vx(j)* sqrt(2.)
             star(5,i) = -vy(j)* sqrt(2.)
             star(6,i) = -vz(j)* sqrt(2.)
          enddo


       endif 
c++     

C     setting up header parameter (INTEGER)
      
      
      ih(1) = 1                 ! number of galaxies
      ih(2) = 8*nbodies-7*ndisk           ! number of particles (simulation)
      ih(3) = 0                 ! current integration step
      ih(4) = 5                 ! last integration step
      ih(5) = model             ! physical units (0) or model (1) 
      ih(6) = 0                 ! center of mass
      ih(7) = 10                ! save data every n-steps, n = ih(7) 
      ih(8) = 0                 ! no sticky particle code
      ih(9) = 0                 ! no black hole code : = 0 no BH. 
      ih(10)= 0 
      ih(11)= 0                 ! flag for comoving (1) or fixed (0) coord.
      ih(12)= 0                 ! constant (0) or adjustable (1) dt 
      ih(13)= unitl             ! unit length = pc (0) or kpc (1) 

      do i = 14,60
         ih(i) = 0
      enddo

C     setting up header parameter (REAL) : all data in model units. 
      
      fh(1) = gi                ! grid radius inner grid
      fh(2) = gm                ! grid radius middle grid
      fh(3) = go                ! grid radius outer grid
      fh(4) = dt                ! dt
      fh(5) = mmtot             ! total mass of galaxy
      fh(6) = real(8*nbodies-7*ndisk)     ! number of particles (this galaxy) !w/o disk!!!!!
      fh(7) = real(8*nbodies-7*ndisk)     ! number of particles left          !w/o disk!!!!! 
      do j = 1,6
         fh( 7+j) = cms(j)      ! center of mass and velocity
         fh(13+j) = cms(j)      ! focus of grids (center of mass)
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

c     Basic properties of the galaxy as given here - 
      
      fh(51) = tcr*d 
      fh(52) = mmtot 
      fh(53) = mrlim 

c     Store scaling factors for conversion to physical units 

      fh(54) = gp 
      fh(55) = a 
      fh(56) = b 
      fh(57) = c
      fh(58) = d 
      fh(59) = cv 
      fh(60) = kmspcy 

c     save data onto disk : 

      call makeext (fname,lfname,'.CONT') ! creating filename
      call makeext (fname,backupname,'.STARTUP') ! creating filename
      
      write (*,*) '  '
      write (*,FMT='($,a)') 'writing data to disc ...'
            
      open (1,FILE=lfname,ACCESS='DIRECT',RECL=240)
      open (2,FILE=backupname,ACCESS='DIRECT',RECL=240)
      
      write (1,rec = 1) ih      ! write INTEGER header
      write (1,rec = 2) fh      ! write REAL header

      close (1)
           
      write (2,rec = 1) ih      ! write INTEGER header
      write (2,rec = 2) fh      ! write REAL header

      close (2)
      
c     reopen file with new recl

      open (1,FILE=lfname,ACCESS='DIRECT',RECL=24)
      open (2,FILE=backupname,ACCESS='DIRECT',RECL=24)

c     write data of all stars          !!!!!!with disk stars!!!!!!!!!!!

      ioff = 20
                                
      k=0
       do i = 1,8*(nbodies-ndisk)+ndisk             ! ih(3) = TOTAL number of stars
          k=k+1
          do j = 1,6
            hstar(j) = star(j,i)
         end do
         write (1,rec = ioff+k) hstar
         write (2,rec = ioff+k) hstar

      end do
         
   
      write(6,*) ' '
      write(6,*) '-------------------------------------'
      write(6,*) 'particle # in simulation',8*nbodies-7*ndisk,'  :) :D '
      WRITE(6,*) '--------------------------------------' 
      
      close (1)
      close (2)

      write (*,'(a,/)' ) 'done'


c++      only one galaxy, please .. 
c++      end do                    ! Galaxy number - ijk loop. 


c     creamos un output con las caracteristicas de la simulacion

      call makeext (fname,backupname,'.read') ! creating filename
      open(20,file=backupname)
      write(20,*) '   disc'
      write(20,*) '   ----'
      write(20,*) 'N particles in disc',ndisk 
      write(20,*) 'thickness',z0
      write(20,*) 'rsolar',rsolar
      write(20,*) 'Q solar',qsolar
      write(20,*) ' '
      write(20,*) '  bulge'
      write(20,*) '  -----'
      write(20,*) 'N particles in bulge (total)',8*nbulge
      write(20,*) 'bulge scale length',abulge
      write(20,*) 'rmax bulge',rmaxbulg
      write(20,*) 'bulge mass',bulgmass
      write(20,*) ' '
      write(20,*) '  halo'
      write(20,*) '  ----'
      write(20,*) 'N particles in halo (total)',8*nhalo
      write(20,*) 'halo rcore (planar)',gamhalo2,'   old',gamhalo
      write(20,*) 'halo rcore (perpend)',b
      write(20,*) 'r tidal',rthalo
      write(20,*) 'rmax halo',rmaxhalo
      write(20,*) 'halo mass',halomass
      write(20,*) ' '
      write(20,*) '  tcrossing (t=r/v)'
      write(20,*) '  ---------'
      adisc=1.
      write(20,*) 't disc:',adisc
      write(20,*) 't bulge:',(abulge**3./bulgmass)**0.5
      write(20,*) 't halo:',(gamhalo**3./halomass)**0.5

      return 
      end


c******************************************************************

      SUBROUTINE makeext (f1,f2,ext)
      
      character*(*) f1,f2,ext
      
      do i = 1,40
         f2(i:i) = ' '
      enddo

      do i = 1,len(f1)
         if (f1(i:i) .eq. ' ') goto 10
         f2(i:i) = f1(i:i)
      enddo
      
 10   continue

      do j = 1,len(ext)
         if (ext(j:j) .eq. ' ') goto 20
         f2(i-1+j:i-1+j) = ext(j:j)
      enddo

      f2(i-1+j:i-1+j) = '\n'

 20   continue

      do k = i-1+j,len(f2)
         f2(k:k) = ' '
      enddo

      return

      end


C***********************************************************************
C
C
        SUBROUTINE radacc
C
C
C***********************************************************************
C
C
C   Subroutine to compute radial acceleration on each disk particle.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 smallnum

C=======================================================================

        smallnum=1.e-07

        DO 10 i=1,ndisk
           IF(radcyl(i).GT.smallnum) THEN
              aradcyl(i)=(x(i)*ax(i)+y(i)*ay(i))/radcyl(i)
           ELSE
              WRITE(6,*) 'radacc>> particle :',i,' had rad =',radcyl(i)
              aradcyl(i)=0.
           END IF
 10     CONTINUE

        WRITE(6,*) 'radacc>> Radial acceleration components computed <<' 

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE readbulg
C
C
C***********************************************************************
C
C
C     Subroutine to read-in parameters associated with the bulge.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*1 yesno

C=======================================================================


        WRITE(6,*) ' '
        WRITE(6,*) ' Time to input bulge parameters '
        WRITE(6,*) ' '

        WRITE(6,10)
 10     FORMAT( ' Do you want to include a bulge (y/n)? ',$)
        READ(5,'(a)') yesno

        IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN

           usebulge=.FALSE.
           nbulge=0
           bulgmass=0.0

        ELSE

           usebulge=.TRUE.

           WRITE(6,20)
 20        FORMAT(' Bulge mass? (Disk mass = 1) ',$)

           READ(5,*) bulgmass

           WRITE(6,30)
 30        FORMAT(' Bulge scale-length? (Disk h = 1) ',$)

           READ(5,*) abulge

           WRITE(6,40)
 40        FORMAT(' Include bulge self-gravity (y/n)? ',$)

           READ(5,'(a)') yesno
 
           IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN
              selfgbul=.FALSE.
              nbulge = 0 
              axibulge=.FALSE.
              bulgerot=.FALSE.
              brotfrac=0.0
           ELSE
              selfgbul=.TRUE.
              WRITE(6,50)
 50           FORMAT(' Number of bulge particles? ',$)

              READ(5,*) nbulge

              WRITE(6,52)
 52           FORMAT(' Maximum radius for bulge particles? ',$)

              READ(5,*) rmaxbulg

              WRITE(6,55)
 55           FORMAT(' Softening length for bulge particles? ',$)

              READ(5,*) epsbulge

              WRITE(6,60)
 60           FORMAT(' Do you want a non-spherical bulge (y/n)? ',$)

              READ(5,'(a)') yesno

              IF(yesno.EQ.'Y'.OR.yesno.EQ.'y') THEN

                 axibulge=.TRUE.

 62              WRITE(6,65)
 65              FORMAT(' Input c-bulge: ',$)

                 READ(5,*) cbulge

                 IF(cbulge.GT.abulge) THEN
                    WRITE(6,*) ' c must be <= a '
                    GO TO 62
                 ENDIF

                 WRITE(6,70)
 70              FORMAT(' Input zmaxbulg (suggest c*rmaxbulg/a) ',$)

                 READ(5,*) zmaxbulg

                 WRITE(6,75)
 75              FORMAT(' Input nsimpson ',$)

                 READ(5,*) nsimpson

                 WRITE(6,80)
 80              FORMAT(' Include bulge rotation (y/n)? ',$)

                 READ(5,'(a)') yesno

                 IF(yesno.EQ.'Y'.OR.yesno.EQ.'y') THEN
                    bulgerot=.TRUE.

 82                 CONTINUE

                    WRITE(6,85)
 85                 FORMAT(' Fraction of reversed particles ? ',$)

                    READ(5,*) brotfrac

                    IF(brotfrac.LT.0.5.OR.brotfrac.GT.1.0) THEN
                       WRITE(6,*) ' Frac. must be between 0.5 and 1 '
                       WRITE(6,*) ' '
                       GO TO 82
                    ENDIF

                 ELSE
                    bulgerot=.FALSE.
                    brotfrac=0.0
                 ENDIF

              ELSE
                 axibulge=.FALSE.
                 bulgerot=.FALSE.
                 brotfrac=0.0
              ENDIF

           ENDIF

        ENDIF

        nbodies=nbodies+nbulge
        IF(nbodies.GT.nbodsmax) CALL berror(' nbulge in readbulg ')
C                                    ------

        RETURN
        END


C***********************************************************************
C
C
        SUBROUTINE readdisk
C
C
C***********************************************************************
C
C
C     Subroutine to read-in parameters associated with the disk.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*1 yesno

C=======================================================================


        WRITE(6,*) ' '
        WRITE(6,*) ' Time to input disk parameters '
        WRITE(6,*) ' '

        WRITE(6,10)
   10	FORMAT(' Number of disk star-particles? ',$)
	
        READ(5,*) ndstars
	
        nbodies=nbodies+ndstars
        IF(nbodies.GT.nbodsmax) CALL berror(' ndstars in readdisk ')
C                                    ------

        WRITE(6,20)
 20     FORMAT(' z0? (Typically, z0 = 0.2) ',$)

        READ(5,*) z0

        WRITE(6,30)
 30     FORMAT(' Solar radius? (Suggest: 8.5/3.5 = 2.428571429) ',$)

        READ(5,*) rsolar

        WRITE(6,40)
 40     FORMAT(' Q-solar? (Typically, Q = 1.5) ',$)

        READ(5,*) qsolar

        WRITE(6,50)
 50     FORMAT(' Softening length for disk particles? ',$)

        READ(5,*) epsdisk

        epsdisk2=epsdisk*epsdisk

        WRITE(6,60)
 60     FORMAT(' zmax? (Suggest 10 * z0) ',$)

        READ(5,*) zmax

        WRITE(6,70)
 70     FORMAT(' rmax? (Suggest 15 * h = 15) ',$)

        READ(5,*) rmax


        WRITE(6,80)
 80     FORMAT(' Do you want to include gas in the disk (y/n)? ',$)

        READ(5,'(a)') yesno

        IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN

           usegas=.FALSE.
           ndgas=0
           gasmass=0.0
           zmaxgas=0.0
           rmaxgas=rmax
           rmingas=0.0

        ELSE

           usegas=.TRUE.

           WRITE(6,90)
 90        FORMAT(' Number of disk gas-particles? ',$)

           READ(5,*) ndgas

           WRITE(6,100)
 100       FORMAT(' Total mass of gas? (Total disk mass = 1 ) ',$)

           READ(5,*) gasmass

           WRITE(6,110)
 110       FORMAT(' Temperature of gas? (Suggest 10^4) ',$)

           READ(5,*) gastemp

           WRITE(6,120)
 120       FORMAT(' z0gas? ',$)

           READ(5,*) z0gas

           WRITE(6,130)
 130       FORMAT(' zmaxgas? (Suggest 10 * z0gas) ',$)

           READ(5,*) zmaxgas

           WRITE(6,140)
 140       FORMAT(' rmaxgas? (Suggest rmax) ',$)

           READ(5,*) rmaxgas

           WRITE(6,150)
 150       FORMAT(' rmingas? ',$)

           READ(5,*) rmingas

           WRITE(6,160)
 160       FORMAT(' Include gas self-gravity? ',$)

           READ(5,'(a)') yesno

           IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN
              selfggas=.FALSE.
           ELSE
              selfggas=.TRUE.
           ENDIF

        ENDIF

        nbodies=nbodies+ndgas
        IF(nbodies.GT.nbodsmax) CALL berror(' ndgas in readdisk ')
C                                    ------
        ndisk=ndgas+ndstars

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE readhalo
C
C
C***********************************************************************
C
C
C     Subroutine to read-in parameters associated with the halo.  For
C     the 'IS' halo, a table is set up containing its properties,
C     defined by:
C
C
C                     M alpha          EXP(-r**2/rt**2)
C        rho(r) = ---------------      ----------------- ,
C                 2 pi**1.5 r_t**3     q**2 + r**2/rt**2
C
C
C     where q = gamhalo / rthalo, and alphhalo = 1/[1 - sqrt(pi) * q *
C     EXP(q**2) * (1 - erf(q)) ].  The cumulative mass and potential
C     for this mass density are:
C
C
C                           r/rt
C                            /
C                 2 M alpha  | x**2 EXP(-x**2)
C        M(r) =   ---------  | ---------------  dx
C                  sqrt(pi)  |   x**2 + q**2
C                            /
C                            0
C
C
C                    GM(r)     G M alpha
C        Phi(r) = -  ----- +  -----------  Ei[ - (r/rt)**2 - q**2 ],
C                      r      sqrt(pi) rt
C
C
C     where Ei is an exponential integral.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*1 yesno
        INTEGER i,j
        REAL*8 one,pi,rootpi,drhalo,qhalo,alphhalo,xj,deltax,expint,
     &         erfcc

C=======================================================================


        WRITE(6,*) ' '
        WRITE(6,*) ' Time to input halo parameters '
        WRITE(6,*) ' '

        WRITE(6,10)
 10     FORMAT(' Do you want to include a halo (y/n)? ',$)

        READ(5,'(a)') yesno

        IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN

           usehalo=.FALSE.
           nhalo=0
           halomass=0.0

        ELSE

           usehalo=.TRUE.

           WRITE(6,20)
 20        FORMAT(' Include halo self-gravity (y/n)? ',$)

           READ(5,'(a)') yesno
 
           WRITE(6,35)
 35        FORMAT(' Maximum radius of halo? ',$)

           READ(5,*) rmaxhalo       

           IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN
              selfghal=.FALSE.
              nhalo = 0 

           ELSE
              selfghal=.TRUE.
              WRITE(6,30)
 30           FORMAT(' Number of halo particles? ',$)

              READ(5,*) nhalo

              WRITE(6,37)
 37           FORMAT(' Softening length for halo particles? ',$)

              READ(5,*) epshalo

           ENDIF

           WRITE(6,40)
 40        FORMAT(' Type of halo (LH or IS)? ',$)

           READ(5,'(a)') halotype

           WRITE(6,50)
 50        FORMAT(' Halo mass? (Disk mass = 1) ',$)

           READ(5,*) halomass

           IF(halotype.EQ.'IS') THEN

              WRITE(6,60)
 60           FORMAT(' Halo core radius (gamma)? (Disk h = 1) ',$)

              READ(5,*) gamhalo

              WRITE(6,70)
 70           FORMAT(' Halo tidal radius (rthalo)? (Disk h = 1) ',$)
              READ(5,*) rthalo

              one=1.0
              pi=4.0*ATAN(one)
              rootpi=SQRT(pi)
              ntabhalo=maxtabh
              print*,'ntabhalo',ntabhalo
              drhalo=rmaxhalo/(ntabhalo-1)

              qhalo=gamhalo/rthalo
              alphhalo=1./(1.-SQRT(pi)*qhalo*EXP(qhalo*qhalo)*
     &                 erfcc(qhalo))
              print*,'alpha',alphhalo

              rhalo(1)=1.e-10

              DO 75 i=2,ntabhalo
                 rhalo(i)=(i-1)*drhalo
 75           CONTINUE

              xmhalo(1)=0.0

              DO 85 i=2,ntabhalo
                 xmhalo(i)=xmhalo(i-1)

                 DO 80 j=1,100
                    deltax=(rhalo(i)-rhalo(i-1))/(100.0*rthalo)
                    xj=rhalo(i-1)/rthalo+j*deltax-0.5*deltax
                    xmhalo(i)=xmhalo(i)+xj**2*EXP(-xj**2)*deltax/
     &                        (xj**2+qhalo**2)
 80              CONTINUE

 85           CONTINUE

              DO 90 i=1,ntabhalo
                 xmhalo(i)=xmhalo(i)*2.0*halomass*alphhalo/rootpi
 90           CONTINUE

              DO 95 i=1,ntabhalo
                 uhalo(i)=-xmhalo(i)/rhalo(i)-halomass*alphhalo*
     &                    EXP(qhalo**2)*expint(rhalo(i)**2/rthalo**2+
     &                    qhalo**2)/(rootpi*rthalo)
 95           CONTINUE
              print*,'xmhalo',xmhalo(ntabhalo)
              IF(ABS(xmhalo(ntabhalo)-halomass)/halomass.GT.1.e-3)
     &           CALL berror(' Halo mass error in readhalo ')
C                     ------
           ELSE

              IF(halotype.EQ.'LH') THEN

                 WRITE(6,110)
 110             FORMAT(' Halo scale-length? (Disk h = 1) ',$)
                 READ(5,*) ahalo

              ELSE
                 CALL berror(' Halo type error in readhalo ')
C                     ------
              ENDIF
           ENDIF
        ENDIF

        nbodies=nbodies+nhalo
        IF(nbodies.GT.nbodsmax) CALL berror(' nhalo in readhalo ')
C                                    ------

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE readsat
C
C
C***********************************************************************
C
C
C     Subroutine to read-in parameters associated with an optional
C     satellite.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*1 yesno
        INTEGER i
        REAL*8 menclose,one,pi,fx,fpx,xold,xnew,rhobar,acoef,rtidal

C=======================================================================


        WRITE(6,*) ' '
        WRITE(6,*) ' Time to input satellite parameters '
        WRITE(6,*) ' '

        WRITE(6,10)
 10     FORMAT( ' Do you want to include a satellite (y/n)? ',$)
        READ(5,'(a)') yesno

        IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN

           usesat=.FALSE.
           nsat=0
           satmass=0.0

        ELSE

           usesat=.TRUE.

           WRITE(6,20)
 20        FORMAT(' Satellite mass? (Disk mass = 1) ',$)

           READ(5,*) satmass

           WRITE(6,30)
 30        FORMAT(' Satellite scale-length? (Disk h = 1) ',$)

           READ(5,*) asat

           WRITE(6,40)
 40        FORMAT(' Satellite coordinates? ',$)

           READ(5,*) xsat,ysat,zsat

C   Determine approximate tidal radius of satellite.
C   ------------------------------------------------
           radsat=SQRT(xsat**2+ysat**2+zsat**2)

           menclose=0.0

           DO 50 i=1,nbodies
              IF(radsph(i).LE.radsat) menclose=menclose+pmass(i)
 50        CONTINUE

           IF(usebulge.AND.(.NOT.selfgbul)) menclose=menclose+
     &        bulgmass*radsat**2/(radsat+abulge)**2

           IF(usehalo.AND.(.NOT.selfghal)) THEN

              IF(halotype.EQ.'LH') menclose=menclose+halomass*
     &           radsat**2/(radsat+ahalo)**2

              IF(halotype.EQ.'IS') THEN

                 DO 60 i=2,ntabhalo
                    IF(radsat.GE.rhalo(i-1).AND.radsat.LT.rhalo(i))
     &                 menclose=menclose+xmhalo(i-1)
 60              CONTINUE

                 IF(radsat.GE.rhalo(ntabhalo)) menclose=menclose+
     &              xmhalo(ntabhalo)

              ENDIF

           ENDIF

           one=1.0
           pi=4.0*ATAN(one)

           rhobar=3.0*menclose/(4.0*pi*radsat**3)
           acoef=9.0*satmass/(4.0*pi*asat**3*rhobar)

           xnew=1.0

 70        xold=xnew
           fx=xold**3+2.0*xold**2+xold-acoef
           fpx=3.0*xold**2+4.0*xold+1.0

           xnew=xold-fx/fpx

           WRITE(6,*) ' xnew = ',xnew

           IF(ABS((xold-xnew)/xold).GT.1.E-6) GO TO 70

           rtidal=xnew*asat

           WRITE(6,*) ' rtidal = ',rtidal

           WRITE(6,80)
 80        FORMAT(' Maximum satellite radius? (Disk h = 1) ',$)

           READ(5,*) rmaxsat

           WRITE(6,90)
 90        FORMAT(' Include satellite self-gravity (y/n)? ',$)

           READ(5,'(a)') yesno
 
           IF(yesno.EQ.'N'.OR.yesno.EQ.'n') THEN
              selfgsat=.FALSE.
              nsat = 0 

              CALL berror(' Rigid satellite option not implemented ')
C                  ------

           ELSE
              selfgsat=.TRUE.
              WRITE(6,100)
 100          FORMAT(' Number of satellite particles? ',$)

              READ(5,*) nsat

              WRITE(6,105)
 105          FORMAT(' Softening length of satellite particles? ',$)

              READ(5,*) epssat
           ENDIF

        ENDIF

        nbodies=nbodies+nsat
        IF(nbodies.GT.nbodsmax) CALL berror(' nsat in readsat ')
C                                    ------

        RETURN
        END
C***********************************************************************
C
C
        FUNCTION rgauss(xmean,sigma)
C
C
C***********************************************************************
C
C
C   Function to produce a Gaussianly distributed random number.
C
C
C=======================================================================

C   Declaration of local variables.
C   -------------------------------

        REAL*8 sigrange,a,x,b,p,rgauss,xmean,sigma
        REAL*8 dxrand

C=======================================================================

        sigrange=5.

  10    CONTINUE

        a=dxrand()
        x=2.*(a-0.5)*sigma*sigrange + xmean
        b=dxrand()
        p=EXP(-(x-xmean)*(x-xmean)/(2.*sigma*sigma))

        IF(b.gt.p) GO TO 10

        rgauss=x

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE satvel
C
C
C***********************************************************************
C
C
C   Subroutine to initialize satellite bulk velocity.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        REAL*8 smallnum,aradsat,vcircsat,vescsat

C=======================================================================

        CALL forceds
C            -------
        IF(usebulge) CALL forcebs
C                         -------
        IF(usehalo) CALL forcehs
C                        -------

        smallnum=1.e-07

        IF(radsat.GT.smallnum) THEN
           aradsat=(xsat*axsat+ysat*aysat+zsat*azsat)/radsat
        ELSE
           WRITE(6,*) 'satvel>> satellite had rad =',radsat
           aradsat=0.
        END IF

        vcircsat=SQRT(radsat*ABS(aradsat))
        vescsat=SQRT(2.0*ABS(potsat))

        WRITE(6,*) ' Satellite circular velocity = ',vcircsat
        WRITE(6,*) ' Satellite escape velocity = ',vescsat

        WRITE(6,35)
 35     FORMAT(' Satellite velocities? ',$)

        READ(5,*) vxsat,vysat,vzsat

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE setbulge
C
C
C***********************************************************************
C
C
C   Establish a bulge corresponding to model analyzed by Hernquist
C   (1991).   This subroutine initializes only the spatial coordinates;
C   velocities are set in bulgevel.f.  The variable const represents 
C   the maximum value of r**2 * v**2 * f(q), which occurs at 
C   r/r_0 = 0.63798179 and v/v_g = v_circ.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ntot
        REAL*8 one,twopi,const,xr,xv,e,p,fq,pn,q,cth,sth,signs,phi,
     &         cthv,sthv,phiv,Rtest,ztest,mtest,ftest
        REAL*8 dxrand

C=======================================================================

        ntot=0

        one=1.0
        twopi=2.*4.*ATAN(one)
        const=1.20223157581242

        IF(.NOT.axibulge) THEN

 10        xr=dxrand()*rmaxbulg
           xv=dxrand()*SQRT(2.*bulgmass/abulge)
           e=0.5*xv*xv-bulgmass/(xr+abulge)
           IF(e.GT.0.0.OR.e.LT.-bulgmass/abulge) GO TO 10

           q=SQRT(-abulge*e/bulgmass)

           p=dxrand()

           fq=(3.*ASIN(q)+q*SQRT(1.-q*q)*(1.-2.*q*q)*(8.*q**4-
     &        8.*q**2-3.))/(1.-q**2)**2.5
           pn=(xr*xr/(abulge*abulge))*(xv*xv/(bulgmass/abulge))*fq/const

           IF(pn.GT.1.0) THEN
              WRITE(6,*) ' pn error in setbulge '
              WRITE(6,66) pn
 66           FORMAT(1x,1pe16.8)
              pn=1.0
           ENDIF

           IF(p.LE.pn) THEN

              cth=2.*(dxrand()-.5)
              sth=SQRT(1.-cth*cth)
              signs=2.*(dxrand()-.5)
              cth=signs*cth/ABS(signs)
              phi=twopi*dxrand()

              cthv=2.*(dxrand()-.5)
              sthv=SQRT(1.-cthv*cthv)
              signs=2.*(dxrand()-.5)
              cthv=signs*cthv/ABS(signs)
              phiv=twopi*dxrand()

              ntot=ntot+1

              x(ntot+nbodies-nbulge)=xr*sth*COS(phi)
              y(ntot+nbodies-nbulge)=xr*sth*SIN(phi)
              z(ntot+nbodies-nbulge)=xr*cth
              vx(ntot+nbodies-nbulge)=xv*sthv*COS(phiv)
              vy(ntot+nbodies-nbulge)=xv*sthv*SIN(phiv)
              vz(ntot+nbodies-nbulge)=xv*cthv

           ELSE
              GO TO 10
           ENDIF

           IF(ntot.LT.nbulge) GO TO 10

        ELSE

 100       CONTINUE

           Rtest=rmaxbulg*dxrand()
           ztest=zmaxbulg*dxrand()

           p=dxrand()
           IF(p.LE.0.5) ztest = -ztest

           mtest=SQRT(Rtest**2/abulge**2+ztest**2/cbulge**2)
           ftest=Rtest/(abulge*abulge*cbulge*mtest*(1.0+mtest)**3)
           ftest=ftest*abulge*cbulge

           p=dxrand()

           IF(p.LE.ftest) THEN

              ntot=ntot+1

              z(ntot+nbodies-nbulge)=ztest

              p=dxrand()

              x(ntot+nbodies-nbulge)=Rtest*COS(twopi*p)
              y(ntot+nbodies-nbulge)=Rtest*SIN(twopi*p)

              vx(ntot+nbodies-nbulge)=0.0
              vy(ntot+nbodies-nbulge)=0.0
              vz(ntot+nbodies-nbulge)=0.0

           ELSE
              GO TO 100
           ENDIF

           IF(ntot.LT.nbulge) GO TO 100

        ENDIF

        WRITE(6,*) 'setbulge>> Bulge established <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE setdisk
C
C
C***********************************************************************
C
C
C   Establish an exponental surface density distribution with scale 
C   length h and constant scale height z0.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n,m
        REAL*8 pi,one,phi,p,r,r2,anum,zp,zps,bnum,pz
        REAL*8 dxrand

C=======================================================================

        one=1.0
        pi=4.0*ATAN(one)

        n=0

 10     CONTINUE

        r2=dxrand()*rmaxgas*rmaxgas

        IF(r2.LT.rmingas**2) GO TO 10

        phi=dxrand()*2.*pi
        r=SQRT(r2)
        p=EXP(-r/h)
        anum=dxrand()

        IF(anum.gt.p) GO TO 10

        n=n+1

        x(n)=r*COS(phi)
        y(n)=r*SIN(phi)

        IF(n.LT.ndgas) GO TO 10

 15     CONTINUE
       
        r2=dxrand()*rmax*rmax
        phi=dxrand()*2.*pi
        r=SQRT(r2)
        p=EXP(-r/h)
        anum=dxrand()
	
        IF(anum.GT.p) GO TO 15

        n=n+1
        x(n)=r*COS(phi)
        y(n)=r*SIN(phi)
        
        IF(n.LT.ndisk) GO TO 15

        m=0

 20     CONTINUE

        IF(m.LT.ndgas) THEN
           zp=2.*(dxrand()-0.5)*zmaxgas
           zps=zp/z0gas
        ELSE
           zp=2.*(dxrand()-0.5)*zmax
           zps=zp/z0
        ENDIF

        bnum=dxrand()

C   Assign pz according to sech^2 distribution.
C   -------------------------------------------

        pz=4./((EXP(zps)+EXP(-zps))*(EXP(zps)+EXP(-zps)))

        IF(bnum.GT.pz) GO TO 20

        m=m+1
        z(m)=zp

        IF(m.LT.ndisk) GO TO 20

        WRITE(6,*) 'setdisk>> Disk established <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE sethalo
C
C
C***********************************************************************
C
C
C   Establish a halo corresponding either to the model analyzed by 
C   Hernquist (1991) or a truncated, non-singular isothermal sphere.
C   This subroutine initializes only the spatial coordinates; velocities
C   are set in halovel.f.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ntot
        REAL*8 one,twopi,const,xr,xv,e,p,fq,pn,q,cth,sth,signs,phi,
     &         cthv,sthv,phiv
        REAL*8 dxrand

C=======================================================================

        IF(halotype.EQ.'LH') THEN

C   Establish a halo corresponding to model analyzed by Hernquist
C   (1991).  The variable const represents the maximum value of 
C   r**2 * v**2 * f(q), which occurs at r/r_0 = 0.63798179 and 
C   v/v_g = v_circ.

           ntot=0

           one=1.0
           twopi=2.*4.*ATAN(one)
           const=1.20223157581242

 10        xr=dxrand()*rmaxhalo
           xv=dxrand()*SQRT(2.*halomass/ahalo)
           e=0.5*xv*xv-halomass/(xr+ahalo)
           IF(e.GT.0.0.OR.e.LT.-halomass/ahalo) GO TO 10

           q=SQRT(-ahalo*e/halomass)

           p=dxrand()

           fq=(3.*ASIN(q)+q*SQRT(1.-q*q)*(1.-2.*q*q)*(8.*q**4-
     &        8.*q**2-3.))/(1.-q**2)**2.5
           pn=(xr*xr/(ahalo*ahalo))*(xv*xv/(halomass/ahalo))*fq/const

           IF(pn.GT.1.0) THEN
              WRITE(6,*) ' pn error in sethalo '
              WRITE(6,66) pn
 66           FORMAT(1x,1pe16.8)
              pn=1.0
           ENDIF

           IF(p.LE.pn) THEN

              cth=2.*(dxrand()-.5)
              sth=SQRT(1.-cth*cth)
              signs=2.*(dxrand()-.5)
              cth=signs*cth/ABS(signs)
              phi=twopi*dxrand()

              cthv=2.*(dxrand()-.5)
              sthv=SQRT(1.-cthv*cthv)
              signs=2.*(dxrand()-.5)
              cthv=signs*cthv/ABS(signs)
              phiv=twopi*dxrand()

              ntot=ntot+1

              x(ntot+nbodies-nhalo)=xr*sth*COS(phi)
              y(ntot+nbodies-nhalo)=xr*sth*SIN(phi)
              z(ntot+nbodies-nhalo)=xr*cth
              vx(ntot+nbodies-nhalo)=xv*sthv*COS(phiv)
              vy(ntot+nbodies-nhalo)=xv*sthv*SIN(phiv)
              vz(ntot+nbodies-nhalo)=xv*cthv

           ELSE
              GO TO 10
           ENDIF

           IF(ntot.LT.nhalo) GO TO 10

        ELSE

C   Establish a halo corresponding to a truncated isothermal
C   model.  The variable const represents the maximum value of 
C   r**2 * rho(r), which occurs at r**2 = 0.5*(-gamhalo**2+
C   SQRT(gamhalo**4 + 4 * r_t**2 * gamhalo**2)).

           ntot=0

           one=1.0
           twopi=2.*4.*ATAN(one)
           const=SQRT(0.5*(-gamhalo**2+SQRT(gamhalo**4+4.0*gamhalo**2*
     &           rthalo**2)))
           const=const**2*EXP(-const**2/rthalo**2)/(1.+const**2/
     &           gamhalo**2)

 30        xr=dxrand()*rmaxhalo

           p=dxrand()

           fq=xr**2*EXP(-xr**2/rthalo**2)/(1.+xr**2/gamhalo**2)

           pn=fq/const

           IF(pn.GT.1.0) THEN
              WRITE(6,*) ' pn error in sethalo '
              WRITE(6,66) pn
              pn=1.0
           ENDIF

           IF(p.LE.pn) THEN

              cth=2.*(dxrand()-.5)
              sth=SQRT(1.-cth*cth)
              signs=2.*(dxrand()-.5)
              cth=signs*cth/ABS(signs)
              phi=twopi*dxrand()

              ntot=ntot+1

              x(ntot+nbodies-nhalo)=xr*sth*COS(phi)
              y(ntot+nbodies-nhalo)=xr*sth*SIN(phi)
              z(ntot+nbodies-nhalo)=xr*cth
              vx(ntot+nbodies-nhalo)=0.0
              vy(ntot+nbodies-nhalo)=0.0
              vz(ntot+nbodies-nhalo)=0.0

           ELSE
              GO TO 30
           ENDIF

           IF(ntot.LT.nhalo) GO TO 30

        ENDIF

        WRITE(6,*) 'sethalo>> Halo established <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE setmesh(nphi,xmesh,ymesh,rmesh)
C
C
C***********************************************************************
C
C
C   Subroutine to set up a mesh of x and y points to evaluate the     
C   acceleration to compute the acceleration gradient.         
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j,n,nphi
        REAL*8 xmesh(1),ymesh(1),one,pi,delphi,delr,r,phi,rmesh(1)

C=======================================================================

        one=1.0
        pi=4.0*ATAN(one)

        delphi=2.0*pi/FLOAT(nphi)
        delr=MAX(rmax,rmaxgas)/FLOAT(maxtab)

        r=delr
        n=0

        DO 20 i=1,maxtab

           phi=0.

           DO 10 j=1,nphi
              n=n+1
              xmesh(n)=r*COS(phi)
              ymesh(n)=r*SIN(phi)
              rmesh(n)=r
              phi=phi+delphi
 10        CONTINUE

           r=r+delr

 20     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE setrot
C
C
C***********************************************************************
C
C
C   Subroutine to set up centripital, anticlockwise rotation according 
C   to the mean rotation curve.              
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 smallnum

C=======================================================================
           
        smallnum=1.e-07

        DO 10 i=1,ndisk
           IF(radcyl(i).GT.smallnum) THEN
              vx(i)=vx(i)-rotmean(i)*y(i)/radcyl(i)
              vy(i)=vy(i)+rotmean(i)*x(i)/radcyl(i)  
           ELSE
              WRITE(6,*) 'setrot>> particle:',i,' had rad = ',radcyl(i)
           ENDIF
 10     CONTINUE

        WRITE(6,*) 'setrot>> Disk rotation established <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE setsat
C
C
C***********************************************************************
C
C
C   Establish a satellite corresponding to model analyzed by Hernquist
C   (1991).  The variable const represents the maximum value of 
C   r**2 * v**2 * f(q), which occurs at r/r_0 = 0.63798179 and 
C   v/v_g = v_circ.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER ntot
        REAL*8 one,twopi,const,xr,xv,e,p,fq,pn,q,cth,sth,signs,phi,
     &         cthv,sthv,phiv
        REAL*8 dxrand

C=======================================================================

        ntot=0

        one=1.0
        twopi=2.*4.*ATAN(one)
        const=1.20223157581242

 10     xr=dxrand()*rmaxsat
        xv=dxrand()*SQRT(2.*satmass/asat)
        e=0.5*xv*xv-satmass/(xr+asat)
        IF(e.GT.0.0.OR.e.LT.-satmass/asat) GO TO 10

        q=SQRT(-asat*e/satmass)

        p=dxrand()

        fq=(3.*ASIN(q)+q*SQRT(1.-q*q)*(1.-2.*q*q)*(8.*q**4-8.*q**2-3.))/
     &     (1.-q**2)**2.5
        pn=(xr*xr/(asat*asat))*(xv*xv/(satmass/asat))*fq/const

        IF(pn.GT.1.0) THEN
           WRITE(6,*) ' pn error in setsat '
           WRITE(6,66) pn
 66        FORMAT(1x,1pe16.8)
           pn=1.0
        ENDIF

        IF(p.LE.pn) THEN
           cth=2.*(dxrand()-.5)
           sth=SQRT(1.-cth*cth)
           signs=2.*(dxrand()-.5)
           cth=signs*cth/ABS(signs)
           phi=twopi*dxrand()
           cthv=2.*(dxrand()-.5)
           sthv=SQRT(1.-cthv*cthv)
           signs=2.*(dxrand()-.5)
           cthv=signs*cthv/ABS(signs)
           phiv=twopi*dxrand()
           ntot=ntot+1
           x(ntot+nbodies-nsat)=xr*sth*COS(phi)
           y(ntot+nbodies-nsat)=xr*sth*SIN(phi)
           z(ntot+nbodies-nsat)=xr*cth
           vx(ntot+nbodies-nsat)=xv*sthv*COS(phiv)
           vy(ntot+nbodies-nsat)=xv*sthv*SIN(phiv)
           vz(ntot+nbodies-nsat)=xv*cthv
        ELSE
           GO TO 10
        ENDIF

        IF(ntot.LT.nsat) GO TO 10

        WRITE(6,*) 'setsat>> Satellite established <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE setsigma
C
C
C***********************************************************************
C
C
C   Subroutine to set up the velocity dispersion, correcting for the 
C   radial dispersion near the center.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 smallnum,vmean,gaussz,rgauss,gaussr,gaussphi

C=======================================================================

        smallnum=1.e-07

        vmean=0.

        DO 10 i=1,ndisk
           vx(i)=0.0
           vy(i)=0.0
           vz(i)=0.0
           IF(sigz(i).NE.0.0) THEN
              gaussz=rgauss(vmean,sigz(i))
              vz(i)=gaussz
           ELSE
              IF(i.GT.ndgas) THEN
                 WRITE(6,*) 'setv>> Particle :',i,' had sigz =',sigz(i)
              ENDIF
           ENDIF

           IF(radcyl(i).GT.smallnum) THEN
              IF(sigr(i).NE.0.0) THEN
                 gaussr=rgauss(vmean,sigr(i))
                 gaussphi=rgauss(vmean,sigphi(i))
                 vx(i)=gaussr*(x(i)/radcyl(i))-gaussphi*(y(i)/radcyl(i))
                 vy(i)=gaussr*(y(i)/radcyl(i))+gaussphi*(x(i)/radcyl(i))
              ELSE
                 IF(qsolar.NE.0.) THEN
                    IF(i.GT.ndgas) THEN
                       WRITE(6,*) 'setv>> particle:',i,' had sigr=',
     &                            sigr(i)
                    ENDIF
                 END IF
              END IF 
           ELSE
              WRITE(6,*) 'setv>>  particle :',i,' had rad = ',radcyl(i)
           END IF 

 10     CONTINUE

        WRITE(6,*) 'setsigma>> Velocity dispersion set up <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE sigalar
C
C
C***********************************************************************
C
C
C   Subroutine to compute the critical Toomre radial velocity     
C   dispersion at radius R for a given Q from
C
C                             pi Q g SIGMA(R)
C           sigma[R]_crit ~ -----------------
C                               KAPPA(R)
C
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 pi,one,term

C=======================================================================

        one=1.0
        pi=4.0*ATAN(one)
        term=qsolar*pi

        DO 10 i=1,ndisk
           IF(kappa(i).NE.0.) THEN
              sigt(i)=term*surfd(i)/kappa(i)
           ELSE
              WRITE(6,*) 'sigalar>> particle:',i,' had kappa=',kappa(i)
              sigt(i)=0.
           ENDIF
  10    CONTINUE

        WRITE(6,*) 'sigalar>> Toomre velocity dispersion computed <<'

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE sigcheck
C
C
C***********************************************************************
C
C
C   Subroutine to find the radius at which the mean motion becomes 
C   imaginary. This radius "a" is used to correct the radial velocity 
C   dispersion profile.  The radius where THETA_mean^2 goes negative 
C   (R=a) is taken as the scale length of the correcting function:
C
C        <pi^2>^(1/2)_new = <pi^2>(r=0)_old exp(-sqrt(r^2+2a^2)/2h).
C
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER nrings,i,nfound,j,isis
        REAL*8 delr,one,pi,r,fnsum,omega,fnmean,sign,asign

C=======================================================================

        delr=0.25

        acorrgas=0.0

        one=1.0
        pi=4.0*ATAN(one)
        r=rmax
        nrings=rmax/delr
 
        DO 10 i=1,nrings-1
           fnsum=0.
           nfound=0

           DO 20 j=ndgas+1,ndisk
              IF(radcyl(j).LE.r.AND.radcyl(j).GT.r-delr) THEN
                 nfound=nfound+1
                 omega=rotcirc(j)/radcyl(j)
                 fnsum=fnsum+1.-2.*(radcyl(j)/h)+(rotcirc(j)*
     &                 rotcirc(j)/(sigr(j)*sigr(j)))-0.25*
     &                 kappa(j)*kappa(j)/(omega*omega)
              END IF
 20        CONTINUE

           fnmean=0.

           IF(nfound.GT.0) fnmean=fnsum/FLOAT(nfound)

           IF(fnmean.GE.0.0.OR.nfound.EQ.0)THEN
              sign=1.
           ELSE
              sign=-1.
           ENDIF

           IF(i.eq.1)THEN
              asign=sign
           ELSE
              IF(asign.ne.sign)THEN
                 acorr=r-delr/2.
                 WRITE(6,*) 'sigcheck>> Critical radius acorr : ',acorr
                 GOTO 101
              ELSE
                 asign=sign
              ENDIF
           ENDIF

          r=r-delr

  10    CONTINUE 

        WRITE(6,*) 'sigcheck>> No critical radius found, acorr=0 <<'
        acorr=r-delr/2.

 101    CONTINUE

        DO 102 isis=ndgas+1,ndisk
           sigr(isis)=sigr0*EXP(-(SQRT(radcyl(isis)*radcyl(isis)+
     &                2.0*acorr*acorr)/(2.0*h)))
           sigz(isis)=SQRT(pi*g*z0*surfd0)*EXP(-(SQRT(radcyl(isis)*
     &                radcyl(isis)+2.0*acorr*acorr)/(2.0*h)))
 102    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE sigmap
C
C
C***********************************************************************
C
C
C   Subroutine to compute the azimuthal dispersion for each particle
C   given some functional form for the way the shape of the velocity 
C   ellipsoid (<theta^2>/<pi^2>) varies with radius (C(x), x=R/h).
C   The function C is constrained to be one at the origin and approach 
C   0.5 at large R where we approach the epicyclic limit of small 
C   perturbations on circular orbits for a flat rotation curve 
C   potential.
C
C   At present the C function is that given by epicyclic theory:
C
C           C = <theta^2>/<pi^2> = 0.25 KAPPA^2/OMEGA^2
C
C   where OMEGA is the angular rotational frequency in the disk and
C   halo field.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 sigratio,cfunc

C=======================================================================
        
        DO 10 i=1,ndisk
           sigratio=SQRT(cfunc(rotcirc(i),radcyl(i),kappa(i)))
           IF(sigratio.GT.1.0) THEN
C***          This condition is reported in meanrot  ***
              sigratio=1.
           ENDIF
           sigphi(i)=sigratio*sigr(i)
 10     CONTINUE

        WRITE(6,*) 'sigmap>> Phi velocity dispersion computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE sigmar
C
C
C***********************************************************************
C
C
C   Subroutine to compute the radial velocity dispersion as a function
C   of radius. The radial dispersion drops like exp(-R/(2*h)) in 
C   accordance with the van der Kruit and Searle conjecture and direct 
C   observations of the disk of the Milky Way by Freeman and Lewis.  
C   The dispersion is normalized to that required by the Toomre formula
C   for axisymmetric stability at the solar radius for a given value 
C   of Q.    
C
C   Note : This dispersion profile will be modified near the center of 
C          the disk to insure no net radial motion and a real mean 
C          motion.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,nsum
        REAL*8 rtoll,rminsol,rmaxsol,rsum,sigmean,sigzero

C=======================================================================

        rtoll=0.25
 
        rminsol=rsolar-rtoll
        rmaxsol=rsolar+rtoll

        rsum=0.
        nsum=0

        DO 10 i=1,ndisk
           IF(radcyl(i).GT.rminsol.AND.radcyl(i).LE.rmaxsol) THEN
              nsum=nsum+1
              rsum=rsum+sigt(i)
           ENDIF
 10     CONTINUE
	print*,'nsum',nsum
	print*,'rsum',rsum
        sigmean=rsum/FLOAT(nsum)
        sigzero=sigmean/EXP(-rsolar/(2.0*h))
        sigr0=sigzero

        WRITE(6,*) 'sigmar>> Reference dispersion : ',sigmean
        WRITE(6,*) 'sigmar>> Central Toomre dispersion : ',sigzero

        DO 20 i=1,ndisk
           IF(i.LE.ndgas) THEN
              sigr(i)=0.0
           ELSE
              sigr(i)=sigzero*exp(-radcyl(i)/(2.0*h))
           ENDIF
 20     CONTINUE

        WRITE(6,*) 'sigmar>> Radial velocity dispersion computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE sigmaz
C
C
C***********************************************************************
C
C
C   Subroutine to compute the z velocity dispersion.  For the 
C   isothermal sheet, the z dispersion is given by :
C
C           sigma[z](z=0) = sqrt [ pi g SIGMA(R) z0 ].
C
C   NOTE: The Z dispersion is modified by the same function as the 
C         radial dispersion profile so the ratio of Z to R dispersion 
C         is everywhere equal.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 pi,one,term

C=======================================================================

        one=1.0
        pi=4.0*ATAN(one)
        term=pi*z0

        DO 10 i=1,ndisk
           IF(i.LE.ndgas) THEN
              sigz(i)=0.0
           ELSE
              sigz(i)=SQRT(term*surfd(i))
           ENDIF
 10     CONTINUE

        WRITE(6,*) 'sigmaz>> Z velocity dispersion computed <<'

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE stackmod
C
C
C***********************************************************************
C
C
C     Subroutine to optionally stack together two DBH models on a
C     parabolic orbit.  By convention, the orbital plane is the z = 0
C     plane.
C
C
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*1 yesno
        INTEGER i
        REAL*8 rp,rsep,rmod1,rmod2,mtotal,vsep

C=======================================================================

        addmods=.FALSE.

        IF(usesat) RETURN

        WRITE(6,10)
 10     FORMAT( ' Do you want to add two models (y/n)? ',$)
        READ(5,'(a)') yesno

        IF(yesno.EQ.'N'.OR.yesno.EQ.'n') RETURN

        IF(usehalo.AND.(.NOT.selfghal)) THEN
           WRITE(6,*) ' Sorry, self-gravitating halo required '
           RETURN
        ENDIF

        IF(usebulge.AND.(.NOT.selfgbul)) THEN
           WRITE(6,*) ' Sorry, self-gravitating bulge required '
           RETURN
        ENDIF

        addmods=.TRUE.

        WRITE(6,20)
 20     FORMAT( ' Pericenter distance for parabolic orbit? ',$)
        READ(5,*) rp

        WRITE(6,30)
 30     FORMAT( ' Initial center of mass separation? ',$)
        READ(5,*) rsep

        rmod1=rsep/2.0
        rmod2=rsep/2.0

        mtotal=0.0

        DO 40 i=1,nbodies
           mtotal=mtotal+pmass(i)
 40     CONTINUE

        xmod1=rmod1-rp
        ymod1=SQRT(2.0*rmod1*rp-rp*rp)
        zmod1=0.0

        vxmod1= -SQRT(mtotal*(2.0*rmod1-rp))/(2.0*rmod1)
        vymod1= -SQRT(mtotal*rp)/(2.0*rmod1)
        vzmod1=0.0

        xmod2=-xmod1
        ymod2=-ymod1
        zmod2=0.0

        vxmod2= -vxmod1
        vymod2= -vymod1
        vzmod2=0.0

        vsep=SQRT((vxmod1-vxmod2)**2+(vymod1-vymod2)**2+
     &            (vzmod1-vzmod2)**2)

        WRITE(6,*) ' vsep = ',vsep

        WRITE(6,50)
 50     FORMAT( ' Rotation angles (degrees) theta, phi for disk 1? ',$)
        READ(5,*) thetmod1,phimod1

        WRITE(6,60)
 60     FORMAT( ' Rotation angles (degrees) theta, phi for disk 2? ',$)
        READ(5,*) thetmod2,phimod2

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE surfden
C
C
C***********************************************************************
C
C                                                                 
C   Subroutine to compute the surface density corresponding to each
C   particle.
C
C                                                                 
C=======================================================================

        INCLUDE 'buildgal.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 one,pi,corr

C=======================================================================

        one=1.0
        pi=4.0*ATAN(one)

        corr=1.-EXP(-rmax/h)*(1.+rmax/h) 

        surfd0=diskmass/(2.*pi*h*h*corr)

        DO 10 i=1,ndisk
           surfd(i)=surfd0*EXP(-radcyl(i)/h)
 10     CONTINUE

        WRITE(6,*) 'surfden>> Disk surface density computed <<'

        RETURN
        END
c
        function dxrand()

        INCLUDE 'buildgal.h'
	real*8 dxrand
        real*4 ran1
        integer iseed
        dxrand=dble(ran1(iseed))
         return
        end

      function ran1(idum)
c************************************************************
c                                                           *
c  Uniform random deviate (from Numerical Recipes 2)        *
c  This rountine must be called with a negative seed        *
c     and then a POSITIVE number from then on               *
c                                                           *
c************************************************************
 
      implicit none
 
      integer idum, ia, im, iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
 
      parameter (ia=16807, im=2147483647, am=1./im, iq=127773, ir=2836,
     &     ntab=32, ndiv=1+(im-1)/ntab, eps=1.2e-7, rnmx=1.-eps)
 
      integer j,k,iv(ntab),iy
      save iv, iy
      data iv /ntab*0/, iy /0/
 
      if (idum.le.0 .or. iy.eq.0) then
         idum=max(-idum, 1)
         do j = ntab + 8, 1, -1
            k = idum/iq
            idum = ia*(idum - k*iq) - ir*k
            if (idum.lt.0) idum = idum + im
            if (j.le.ntab) iv(j) = idum
         end do
         iy = iv(1)
      endif
      k = idum/iq
      idum = ia*(idum - k*iq) - ir*k
      if (idum.lt.0) idum = idum + im
      j = 1 + iy/ndiv
      iy = iv(j)
      iv(j) = idum
      ran1 = min(am*iy, rnmx)
 
      RETURN
      END





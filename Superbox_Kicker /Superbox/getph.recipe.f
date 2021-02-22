
c**********************************************************************
c
c     G E T P H . REC. F.
c     Version mit komplett neuer fft - the numerical recipes version 
c     coded by mike fellhauer, ari summer 1999. 
c
c**********************************************************************
c       The purpose of this file is unknown to me (cmb 17.11.1999)
c**********************************************************************
c
      subroutine getph
c
      include 'super.cb'
c
      integer i,j,k,ii
      real    daten(n),eps,soft,qj
c
c     (1.) get the 1/r - numbers in the active grid :
c     ===============================================
c
      soft = 0.0   ! hp(0,0,0) = 1./eps soft could be 0.25

      do i = 0,n-1
         ii = i*i
         do j = 0,n-1
            qj = j*j + ii + soft
            do k = 0,n-1
               if ((qj+k*k).ne.0.0) then
                  hp(i,j,k) = 1. / sqrt(real(qj + k*k))
               else
                  hp(0,0,0) = 1.0 / eps
               endif
            enddo
         enddo
      enddo
c
c     (2.1) 1. Fourier-transformation of hp (z-direction) :
c     ====================================================
c
      do i = 0,n-1
         do j = 0,n-1
c
            do k = 0,n-1
               daten(k+1) = hp(i,j,k)
            enddo
c
            call realft (daten,n,1) 
c
            hp(i,j,0)   = daten(1)
            hp(i,j,nd2) = daten(2)
c
            do k = 1,(nd2-1)
               hp(i,j,    k) = daten(2*k + 1)
               hp(i,j,nd2+k) = daten(2*k + 2)
            enddo
c
         enddo
      enddo
c
c     (2.2) 2. Fourier-transformation of hp (y-direction) :
c     ====================================================
c
      do i = 0,n-1
         do k = 0,n-1
c
            do j = 0,n-1
               daten(j+1) = hp(i,j,k)
            enddo
c
            call realft (daten,n,1)
c
            hp(i,0,k)   = daten(1)
            hp(i,nd2,k) = daten(2)
c
            do j = 1,nd2-1
               hp(i,    j,k) = daten(2*j + 1)
               hp(i,nd2+j,k) = daten(2*j + 2)
            enddo
c
         enddo
      enddo
c
c     (2.3) 3. Fourier-transformation of hp (x-direction) :
c     ====================================================
c
      do j = 0,n-1
         do k = 0,n-1
c
            do i = 0,n-1
               daten (i+1) = hp(i,j,k)
            enddo
c
            call realft (daten,n,1)
c
            hp(0,j,k)   = daten(1)
            hp(nd2,j,k) = daten(2)
c
            do i = 1,nd2-1
               hp(    i,j,k) = daten(2*i + 1)
               hp(nd2+i,j,k) = daten(2*i + 2)
            enddo
c
         enddo
      enddo
c
      return
c     
      end
c     
c******************************************************************
c     end of file
c******************************************************************

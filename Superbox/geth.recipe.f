c**********************************************************************
c
c     G E T H . REC. F.
c     Version mit komplett neuer fft - the numerical recipes version 
c     coded by mike fellhauer, ari summer 1999. 
c
c**********************************************************************
c
c     This file contents the transformation routine for the
c     Greens-function
c
      subroutine geth
c
c     subroutine 'geth' : calculates the fourier-transform of the 
c                         1/r - numbers
c                         it's only necessary to calculate once
c
c     output generated in h(i,j,k) for i,j,k = 0 to n
c
      include 'super.cb'
c
      integer i,j,k,ii, nm2, nd2
      real    daten(2*n),soft,qj
c
c     (1.) get the 1/r - numbers in the active n*n*n grid :
c     =====================================================
c
      soft = 0.0               
      nm2  = 2*n 
      nd2  = n / 2

      ! h(0,0,0) = 1./eps; soft could be 0.25
                              
      do i = 0,n
         ii = i*i
         do j = 0,n
            qj = j*j + ii + soft
            do k = 0,n
               if ((qj+k*k).ne.0.0) then
                  h(i,j,k) = 1. / sqrt(real(qj + k*k))
               else
                  h(0,0,0) = 1.0 / eps
               endif
            enddo
         enddo
      enddo
c
c     (2.1) 1. Fourier-transformation of h (z-direction) :
c     ====================================================
c
      do i = 0,n
         do j = 0,n
c
            do k = 0,n
               daten(k+1) = h(i,j,k)
            enddo
            do k = 1,n-1
               daten(n+k+1) = h(i,j,n-k)
            enddo
c
            call realft (daten,nm2,1) ! see: fft.f
c
            h(i,j,0) = daten(1)
            h(i,j,n) = daten(2)
c
            do k = 1,(n-1)
               h(i,j,k) = daten(2*k + 1)
            enddo
c
         enddo
      enddo
c
c     (2.2) 2. Fourier-transformation of h (y-direction) :
c     ====================================================
c
      do i = 0,n
         do k = 0,n
c
            do j = 0,n
               daten(j+1) = h(i,j,k)
            enddo
            do j = 1,n-1
               daten(n+j+1) = h(i,n-j,k)
            enddo
c
            call realft (daten,nm2,1) ! see: fft.f
c
            h(i,0,k) = daten(1)
            h(i,n,k) = daten(2)
c
            do j = 1,n-1
               h(i,j,k) = daten(2*j + 1)
            enddo
c
         enddo
      enddo
c
c     (2.3) 3. Fourier-transformation of h (x-direction) :
c     ====================================================
c
      do j = 0,n
         do k = 0,n
c
            do i = 0,n
               daten(i+1) = h(i,j,k)
            enddo
            do i = 1,n-1
               daten(n+i+1) = h(n-i,j,k)
            enddo
c
            call realft (daten,nm2,1) ! see: fft.f
c
            h(0,j,k) = daten(1)
            h(n,j,k) = daten(2)
c
            do i = 1,n-1
               h(i,j,k) = daten(2*i + 1)
            enddo
c
         enddo
      enddo
c
      return
c     
      end
c     
c**********************************************************************
c
      subroutine getph
c
      include 'super.cb'
c
      integer i,j,k,ii
      real    daten(n),soft,qj
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
            hp(i,j,n/2) = daten(2)
c
            do k = 1,(n/2-1)
               hp(i,j,    k) = daten(2*k + 1)
               hp(i,j,n/2+k) = daten(2*k + 2)
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
            hp(i,n/2,k) = daten(2)
c
            do j = 1,n/2-1
               hp(i,    j,k) = daten(2*j + 1)
               hp(i,n/2+j,k) = daten(2*j + 2)
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
            hp(n/2,j,k) = daten(2)
c
            do i = 1,n/2-1
               hp(    i,j,k) = daten(2*i + 1)
               hp(n/2+i,j,k) = daten(2*i + 2)
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

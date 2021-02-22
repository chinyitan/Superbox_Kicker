c**********************************************************************
c
c     G E T H . F
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
c     output generated in h(i,j,k) for i,j,k = 0 to 64
c
      include 'super.cb'
c
      integer nm2,nd2,i,j,k,ii
      real    y1(n),y2(n),data(n*2),soft,qj
c
      nm2   =  n*2
      nd2   =  n/2
c     
c     (1.) get the 1/r - numbers in the active 32*32*32 grid :
c     ========================================================
c
c++     This line commented out by cmb 03.06.1999 -> eps in super.cb
c++      eps  = 0.75               ! now there is no additional softening
      soft = 0.0                ! and you should have a value for
                                ! h(0,0,0) = 1./eps
                                ! soft could be 0.25
      do 12 i = 0,n
         ii = i*i
         do 11 j = 0,n
            qj = j*j + ii + soft
            do 10 k = 0,n
               if ((qj+k*k).ne.0.0) then
                  h(i,j,k) = 1. / sqrt(real(qj + k*k))
               else
                  h(0,0,0) = 1.0 / eps
               endif
 10         continue
 11      continue
 12   continue
c
c     (2.1) 1. Fourier-transformation of h (z-direction) :
c     ====================================================
c
      do 30 i = 0,n
         do 25 j = 0,n
c
            do 15 k = 0,(nd2-1)
               y1(k+1) = h(i,j,2*k)
               y2(k+1) = h(i,j,2*k + 1)
 15         continue
c
            y1(nd2+1) = h(i,j,n)
            y2(nd2+1) = h(i,j,n-1)
c
            do 20 k = (nd2+1),(n-1)
               y1(k+1) = h(i,j,nm2 - 2*k)
               y2(k+1) = h(i,j,nm2 - 2*k - 1)
 20         continue
c
            call fft (y1,y2,data,1) ! see: fft.f
c
            h(i,j,0) = data(1)
            h(i,j,n) = data(2)
c
            do 22 k = 1,(n-1)
               h(i,j,k) = data(nm2 - 2*k + 1)
 22         continue
c
 25      continue
 30   continue
c
c     (2.2) 2. Fourier-transformation of h (y-direction) :
c     ====================================================
c
      do 50 i = 0,n
         do 45 k = 0,n
c
            do 35 j = 0,(nd2-1)
               y1(j+1) = h(i,2*j,k)
               y2(j+1) = h(i,2*j + 1,k)
 35         continue
c
            y1(nd2+1) = h(i,n,k)
            y2(nd2+1) = h(i,n-1,k)
c
            do 40 j = (nd2+1),(n-1)
               y1(j+1) = h(i,nm2 - 2*j,k)
               y2(j+1) = h(i,nm2 - 2*j - 1,k)
 40         continue
c
            call fft (y1,y2,data,1) ! see: fft.f
c
            h(i,0,k) = data(1)
            h(i,n,k) = data(2)
c
            do 42 j = 1,n-1
               h(i,j,k) = data(nm2 - 2*j + 1)
 42         continue
c
 45      continue
 50   continue
c
c     (2.3) 3. Fourier-transformation of h (x-direction) :
c     ====================================================
c
      do j = 0,n
         do k = 0,n
c
            do i = 0,(nd2-1)
               y1(i+1) = h(2*i,j,k)
               y2(i+1) = h(2*i+1,j,k)
            enddo
c
            y1(nd2+1) = h(n,j,k)
            y2(nd2+1) = h(n-1,j,k)
c
            do i = (nd2+1),(n-1)
               y1(i+1) = h(nm2 - 2*i,j,k)
               y2(i+1) = h(nm2 - 2*i - 1,j,k)
            enddo
c
            call fft (y1,y2,data,1)  ! see: fft.f
c
            h(0,j,k) = data(1)
            h(n,j,k) = data(2)
c
            do i = 1,n-1
               h(i,j,k) = data(nm2 - 2*i + 1)
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

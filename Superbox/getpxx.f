c**********************************************************************
c
c     G E T X X . F
c
c**********************************************************************
c
      subroutine getpxx (grnr)
c
c     This routine calculates the potential function in each
c     mesh-point of grid grnr via fast fourier transformation
c
c     input : mass-density-array rho, Greensfunction-array h
c     output: potential-array    rho
c
      include 'super.cb'
c
      real    r11(0:n-1,0:2*n-1,0:2*n-1),data(2*n),z(0:2*n-1)
      real    odnm2,y1(n),y2(n)
      integer grnr,nm2,nd2,i,j,k,im2,jm2,km2
c
      nm2   = n*2
      nd2   = n/2
      odnm2 = 1.0 / nm2
c
c     (1.) Fourier-transformation of rho (z-direction) :
c     ==================================================
c
      do 15 i = 0,n-1
         do 14 j = 0,n-1
c
            do 10 k = 0,(nd2-1)
               y1(k + 1)       = rho(i,j,2*k,grnr)
               y2(k + 1)       = rho(i,j,2*k + 1,grnr)
               y1(nd2 + 1 + k) = 0.0
               y2(nd2 + 1 + k) = 0.0
 10         continue
c
            call fft(y1,y2,data,1) ! see: fft.f
c
            r11(i,j,0) = data(1)
            r11(i,j,n) = data(2)
c
            do 12 k = 1,n-1
               r11(i,j,k)     =  data(nm2 - 2*k + 1)
               r11(i,j,n + k) = -data(2*k + 2)
 12         continue
c
 14      continue
 15   continue
c
c     (2.) Fourier-transformation of rho (y-direction) :
c     ==================================================
c
      do 20 i = 0,n-1
         do 19 k = 0,(nm2-1)
c
            do 16 j = 0,(nd2-1)
               y1(j + 1)       = r11(i,2*j,k)
               y2(j + 1)       = r11(i,2*j + 1,k)
               y1(nd2 + 1 + j) = 0.0
               y2(nd2 + 1 + j) = 0.0
 16         continue
c
            call fft(y1,y2,data,1) ! see: fft.f
c
            r11(i,0,k) = data(1)
            r11(i,n,k) = data(2)
c
            do 18 j = 1,n-1
               r11(i,j,k)     =  data(nm2 - 2*j + 1)
               r11(i,n + j,k) = -data(2*j + 2)
 18         continue
c
 19      continue
 20   continue
c
c     (3.) Fourier-transformation of rho (x-direction) :
c     ==================================================
c
      do 50 j = 0,(nm2-1)
         do 49 k = 0,(nm2-1)
c
            do 22 i = 0,(nd2-1)
               y1(i+1)         = r11(2*i,j,k)
               y2(i+1)         = r11(2*i + 1,j,k)
               y1(nd2 + 1 + i) = 0.0
               y2(nd2 + 1 + i) = 0.0
 22         continue
c
            call fft(y1,y2,data,1) ! see: fft.f
c
            z(0) = data(1)
            z(N) = data(2)
c
            do 26 i = 1,n-1
               z(i)     =  data(nm2 - 2*i + 1)
               z(n + i) = -data(2*i + 2)
 26         continue
c
c     (4.) Multiplication with h (= Fourier-transformed Greensfkt.) :
c     ===============================================================
c
            if ((j.gt.n).and.(k.le.n)) goto 30
            if ((j.le.n).and.(k.gt.n)) goto 34
            if ((j.gt.n).and.(k.gt.n)) goto 38
c
            do 28 i = 1,n-1
               z(i)     = z(i)     * h(i,j,k)
               z(n + i) = z(n + i) * h(n - i,j,k)
 28         continue
c
            z(0) = z(0) * h(0,j,k)
            z(n) = z(n) * h(n,j,k)
c
            goto 42
c
 30         continue
c
            do 32 i = 1,n-1
               z(i)     = z(i)     * h(i,nm2 - j,k)
               z(n + i) = z(n + i) * h(n - i,nm2 - j,k)
 32         continue
c
            z(0) = z(0) * h(0,nm2 - j,k)
            z(n) = z(n) * h(n,nm2 - j,k)
c
            goto 42
c
 34         continue
c
            do 36 i = 1,n-1
               z(i)     = z(i)     * h(i,j,nm2 - k)
               z(n + i) = z(n + i) * h(n - i,j,nm2 - k)
 36         continue
c
            z(0) = z(0) * h(0,j,nm2 - k)
            z(n) = z(n) * h(n,j,nm2 - k)
c
            goto 42
c
 38         continue
c
            do 40 i = 1,n-1
               z(i)     = z(i)     * h(i,nm2 - j,nm2 - k)
               z(n + i) = z(n + i) * h(n - i,nm2 - j,nm2 - k)
 40         continue
c
            z(0) = z(0) * h(0,nm2 - j,nm2 - k)
            z(n) = z(n) * h(n,nm2 - j,nm2 - k)
c
 42         continue
c
c     (5.) Re-transformation of rho (x-direction) :
c     =============================================
c
            do 44 i = 1,(nd2-1)
               im2       = i*2
               y1(i + 1)       = z(im2)         + z(nm2 - im2)
               y2(i + 1)       = z(im2 + 1)     + z(nm2 - im2 - 1)
               y1(nd2 + i + 1) = z(n - im2)     - z(n + im2)
               y2(nd2 + i + 1) = z(n - im2 - 1) - z(n + im2 + 1)
 44         continue
c
            y1(1)       = z(0)
            y2(1)       = z(1)     + z(nm2 - 1)
            y1(nd2 + 1) = z(n)
            y2(nd2 + 1) = z(n - 1) - z(n + 1)
c
            call fft(y1,y2,data,-1) ! see: fft.f
c
            r11(0,j,k) = data(1)
c
            do 46 i = 1,n-1
               r11(i,j,k) = odnm2 * (data(nm2 - 2*i + 1) -
     $              data(nm2 - 2*i + 2))
 46         continue
c
 49      continue
 50   continue
c
c     (6.) Re-transformation of rho (y-direction) :
c     =============================================
c
      do 59 i = 0,n-1
         do 58 k = 0,nm2-1
c
            do 52 j = 1,(nd2-1)
               jm2             = j*2
               y1(j + 1)       = r11(i,jm2,k)
     $                         + r11(i,nm2 - jm2,k)
               y2(j + 1)       = r11(i,jm2 + 1,k)
     $                         + r11(i,nm2 - jm2 - 1,k)
               y1(nd2 + j + 1) = r11(i,n - jm2,k)
     $                         - r11(i,n + jm2,k)
               y2(nd2 + j + 1) = r11(i,n - jm2 - 1,k)
     $                         - r11(i,n + jm2 + 1,k)
 52         continue
c
            y1(1)       = r11(i,0,k)
            y2(1)       = r11(i,1,k)     + r11(i,nm2 - 1,k)
            y1(nd2 + 1) = r11(i,n,k)
            y2(nd2 + 1) = r11(i,n - 1,k) - r11(i,n + 1,k)
c
            call fft(y1,y2,data,-1) ! see: fft.f
c
            r11(i,0,k) = data(1)
            r11(i,n,k) = data(2)
c
            do 54 j = 1,n-1
               r11(i,j,k)     = odnm2 * (data(nm2 - 2*j + 1) -
     $              data(nm2 - 2*j + 2))
               r11(i,n + j,k) = odnm2 * (data(2*j + 1) +
     $              data(2*j + 2))
 54         continue
c
 58      continue
 59   continue
c
c     (7.) Re-transformation of rho (z-direction) :
c     =============================================
c
      do 67 i = 0,n-1
         do 66 j = 0,n-1
c
            do 60 k = 1,(nd2-1)
               km2             = k*2
               y1(k + 1)       = r11(i,j,km2)        
     $                         + r11(i,j,nm2 - km2)
               y2(k + 1)       = r11(i,j,km2 + 1)    
     $                         + r11(i,j,nm2 - km2 - 1)
               y1(nd2 + k + 1) = r11(i,j,n - km2)    
     $                         - r11(i,j,n + km2)
               y2(nd2 + k + 1) = r11(i,j,n - km2 - 1)
     $                         - r11(i,j,n + km2 + 1)
 60         continue
c
            y1(1)       = r11(i,j,0)
            y2(1)       = r11(i,j,1)     + r11(i,j,nm2 - 1)
            y1(nd2 + 1) = r11(i,j,n)
            y2(nd2 + 1) = r11(i,j,n - 1) - r11(i,j,n + 1)
c
            call fft(y1,y2,data,-1) ! see:fft.f
c
            rho(i,j,0,grnr) = data(1)
c
            do 62 k = 1,n-1
               rho(i,j,k,grnr) = odnm2 * (data(nm2 - 2*k + 1) -
     $              data(nm2 - 2*k + 2))
 62         continue
c
 66      continue
 67   continue
c
      return
c
      end

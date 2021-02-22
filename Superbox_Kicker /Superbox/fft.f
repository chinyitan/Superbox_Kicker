c**********************************************************************
c
c     F F T . S U P E R. F : the original superbox fast-fourier transform 
C     routine. 
c
c**********************************************************************
c
c     This routine makes a discrete fast fourier transformation
c     the origin of this routine is:
c     H.Werner, R.Schaback : Praktische Mathemathik II
c     Springer Berlin, Heidelberg, New York 1979
c
c       This is the original algorithm used in Fellhauer et al. 1999. 
c       last update: cmb 26.10.1999 
c
      subroutine fft (xreal,ximag,data,isign)
c
      include 'super.cb'
c
      integer ip,nh,i,l2,iadr,l,j1max,j2,j1min,j1,j3,j4,j
      integer i1,i2,i3,i4
      real  w1,w2,yr1,yi1,h1r,h1i,h2r,h2i
      real  xreal(n),ximag(n),yr(n*2),yi(n*2),data(n*2)
c
      ip   = m
      nh   = n/2
c      
      do i=1,n
         yr(i) = xreal(i)
         yi(i) = ximag(i)
      enddo
c
      do i = (n+1),(n*2)
         yr(i) = 0.
         yi(i) = 0.
      enddo
c
      l2   = 1
      iadr = 0
c
      do 50 l = 1,ip
c
         j1max = iadr
         iadr  = n - iadr
         j2    = iadr
c
         do 40 j = 1,nh,l2
c
            j1min = j1max + 1
            j1max = j1max + l2
            w1    = wr(j)
            w2    = wi(j)
c
            do 30 j1 = j1min,j1max
c
               j2     = j2 + 1
               j3     = j1 + nh
               j4     = j2 + l2
               yr(j2) = yr(j1) + yr(j3)
               yi(j2) = yi(j1) + yi(j3)
               yr1    = yr(j1) - yr(j3)
               yi1    = yi(j1) - yi(j3)
               yr(j4) = yr1*w1 - yi1*w2
               yi(j4) = yr1*w2 + yi1*w1
c
 30         continue
c
            j2 = j2 + l2
c
 40      continue
c
         l2 = 2*l2
c
 50   continue
c
      if (iadr.eq.0) goto 70
c
      do 60 j = 1,n
         l     = n + j
         yr(j) = yr(l)
         yi(j) = yi(l)
 60   continue
c
 70   continue
c
      do i = 1,n
         data(i*2-1) = yr(i)    ! xreal(i)=yr(i)
         data(i*2)   = yi(i)    ! ximag(i)=yi(i)
      enddo
c     
      do 11 i = 2,(nh+1)
c
         i2  = 2*i
         i1  = i2   - 1
         i3  = n2p3 - i2
         i4  = i3   + 1
c
         h1r =  c1 * (data(i1) + data(i3))
         h1i =  c1 * (data(i2) - data(i4))
         h2r = -c2 * (data(i2) + data(i4))
         h2i =  c2 * (data(i1) - data(i3))
c
         data(i1) =  h1r + wr1(i)*h2r - wi1(i)*h2i
         data(i2) =  h1i + wr1(i)*h2i + wi1(i)*h2r
         data(i3) =  h1r - wr1(i)*h2r + wi1(i)*h2i
         data(i4) = -h1i + wr1(i)*h2i + wi1(i)*h2r
c
 11   continue
c
      h1r = data(1)
c
      data(1) = h1r + data(2)
      data(2) = h1r - data(2)
c
      return
c
      end
c     
c******************************************************************
c     end of file
c******************************************************************














c      
c******************************************************************
c
c     Calculate all eigenvalues and eigenvectors of a real,
c     symmetric matrix
c
c     source code from : Numerical Recipies In The Physical
c     Sciences 2nd edition
c
      subroutine jacobi(a,n,np,d,v,nrot)
c
      integer n,np,nrot,nmax
      real    a(np,np),d(np),v(np,np)
c
      parameter (nmax = 500)
c
      integer i,ip,iq,j
      real    c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)
c     
      do 12 ip = 1,n
         do 11 iq = 1,n
            v(ip,iq) = 0.
 11      continue
         v(ip,ip) = 1.
 12   continue
c
      do 13 ip = 1,n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.
 13   continue
c
      nrot = 0
c
      do 24 i = 1,50
c
         sm = 0.
c
         do 15 ip = 1,n-1
            do 14 iq = ip+1,n
               sm = sm + abs(a(ip,iq))
 14         continue
 15      continue
c
         if (sm .eq. 0.0) return
c
         if (i .lt. 4) then
            tresh = 0.2 * sm / n**2
         else
            tresh = 0.0
         endif
c
         do 22 ip = 1,n-1
            do 21 iq = ip+1,n
c
               g = 100. * abs(a(ip,iq))
c
               if ((i .gt. 4) .and. (abs(d(ip))+g .eq. abs(d(ip)))
     $              .and. (abs(d(iq))+g .eq. abs(d(iq)))) then
                  a(ip,iq) = 0.0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h        = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t     = a(ip,iq) / h
                  else
                     theta = 0.5 * h / a(ip,iq)
                     t     = 1. / (abs(theta) + sqrt(1.+theta**2))
                     if (theta .lt. 0.0) t = -t
                  endif
c
                  c     = 1. / sqrt(1+t**2)
                  s     = t * c
                  tau   = s / (1.+c)
                  h     = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.
c
                  do 16 j = 1,ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
 16               continue
c
                  do 17 j = ip+1,iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
 17               continue
c
                  do 18 j = iq+1,n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
 18               continue
c
                  do 19 j = 1,n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
 19               continue
c
                  nrot = nrot + 1
c
               endif
c
 21         continue
c
 22      continue
c
         do 23 ip = 1,n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.
 23      continue
c
 24   continue
c
      pause 'too many iterations in jacobi'
c
      return
c
      end

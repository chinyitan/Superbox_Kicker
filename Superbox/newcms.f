c      
c******************************************************************
c
c     Routine to determine center of density of galaxy gn
c
      subroutine newcms (step,gn,istart,istop)
c      
      include 'super.cb'
c      
      integer step,gn,istart,istop,count,hc
      integer i
      double precision  hx(6),hxn(6),r,x,y,z,dist,rmax,maxdist
c      
      rmax    = dble(rcore(gn) * rcore(gn))
      maxdist = 0.1d0 * rmax
c
      do i = 1,6
         hx(i) = 0.d0
      enddo
c
      hc = 0
c      
      do i = istart,istop
c
         if (cflag(i) .eq. 1) then
            hx(1) = hx(1) + dble(star(1,i))
            hx(2) = hx(2) + dble(star(2,i))
            hx(3) = hx(3) + dble(star(3,i))
            hx(4) = hx(4) + dble(star(4,i))
            hx(5) = hx(5) + dble(star(5,i))
            hx(6) = hx(6) + dble(star(6,i))
            hc    = hc + 1
            cflag(i) = 0
         endif
c
      enddo
c      
      if (hc .gt. 0) then
         do i = 1,6
            hx(i) = hx(i) / dble(hc)
         enddo
      else  
         write (99,'(a)') '----------------------------------------'
         write (99,'(a)') 'no stars to determine center of density'
         write (99,'(a)') 'saving data and stop program'
         write (99,'(a)') '----------------------------------------'
         call wdata
         close (99)
         stop
      endif
c
      count = 0
c      
 1    continue
c
      count = count + 1
c
      do i = 1,6
         hxn(i) = 0.d0
      enddo
c
      hc = 0
c      
      do i = istart,istop
c
         x = dble(star(1,i)) - hx(1)
         y = dble(star(2,i)) - hx(2)
         z = dble(star(3,i)) - hx(3)
         r = x*x + y*y + z*z
c
         if (r .lt. rmax) then
            cflag(i) = 1
            hxn(1) = hxn(1) + dble(star(1,i))
            hxn(2) = hxn(2) + dble(star(2,i))
            hxn(3) = hxn(3) + dble(star(3,i))
            hxn(4) = hxn(4) + dble(star(4,i))
            hxn(5) = hxn(5) + dble(star(5,i))
            hxn(6) = hxn(6) + dble(star(6,i))
            hc     = hc + 1
         else
            cflag(i) = 0
         endif
c
      enddo
c
      if (hc .gt. 0) then
         do i = 1,6
            hxn(i) = hxn(i) / dble(hc)
         enddo
      endif
c
      dist = dsqrt((hx(1)-hxn(1))**2 + (hx(2)-hxn(2))**2 +
     $     (hx(3)-hxn(3))**2)
c      
      if (dist .gt. maxdist) then
         do i = 1,6
            hx(i) = hxn(i)
         enddo
         if (count .lt. 5) goto 1
      endif
c      
      do i = 1,6
         dgcms(i,gn) = real(hxn(i))
      enddo
c      
      return
c
      end

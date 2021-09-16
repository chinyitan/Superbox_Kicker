c
c******************************************************************
c
      subroutine getflags (hx)
c      
      include 'super.cb'
c      
      real    hx(6,mg),x,y,z,r,rmax
      integer istart,istop,k,i,icount
c      
      do k = 1,gnum
c
         if (k .eq. 1) then
            istart = 1
            istop  = gstno(1)
         else
            istart = 0
            do i = 1,k-1
               istart = istart + gstno(i)
            enddo
            istop  = istart + gstno(k)
            istart = istart + 1
         endif
c
         rmax   = rcore(k) * rcore(k)
         icount = 0
c
         do i = istart,istop
c
            x = star(1,i) - hx(1,k)
            y = star(2,i) - hx(2,k)
            z = star(3,i) - hx(3,k)
            r = x*x + y*y + z*z
            if (r .lt. rmax) then
               cflag(i) = 1
               icount   = icount + 1
            else
               cflag(i) = 0
            endif
c
         enddo
c
      enddo
c
      return
c      
      end

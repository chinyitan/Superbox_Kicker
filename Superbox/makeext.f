c      
c******************************************************************
c
c     Concanate 'fname' with 'ext'
c
      subroutine makeext (f1,f2,ext)
c      
      character*(*) f1,f2,ext
c      
      do i = 1,100
         f2(i:i) = ' '
      enddo
c
      do i = 1,len(f1)
         if (f1(i:i) .eq. ' ') goto 10
         f2(i:i) = f1(i:i)
      enddo
c
 10   continue
c
      do j = 1,len(ext)
         if (ext(j:j) .eq. ' ') goto 20
         f2(i-1+j:i-1+j) = ext(j:j)
      enddo
c
      f2(i-1+j:i-1+j) = '\n'
c
 20   continue
c
      do k = i-1+j,len(f2)
         f2(k:k) = ' '
      enddo
c
      return
c
      end

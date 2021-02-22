c
c     (C) Copr. 1986-92 Numerical Recipes Software >).
c     
c******************************************************************
c
c     this routine creates the names of the output-files
c     (without ending)
c
      subroutine createna (fname,lfname,ign,ierr)
c      
      character*100 fname,lfname
      integer       ign,izero,ierr,i
c      
      izero  =    0
c      
      do i = 1,100
         lfname(i:i) = ' '
      enddo
c
      do i = 1,100
         if (fname(i:i) .eq. ' ') goto 10
         lfname(i:i) = fname(i:i)
      enddo
      ierr = 1
      goto 11
c
 10   lfname(i:i) = '-'
      lfname(i+1:i+1) = 'g'
c
      if (ign .lt. 10) then
         write (lfname(i+2:i+2),'(i1)') izero
         write (lfname(i+3:i+3),'(i1)') ign
      else
         write (lfname(i+2:i+3),'(i2)') ign
      endif
c
      lfname(i+4:i+4) = '.'
      lfname(i+5:i+5) = ' '
c
      ierr = 0
c
 11   continue
c     
      end

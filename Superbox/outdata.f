c      
c******************************************************************
c
c     Write data (energy, cumulative mass...of all galaxies
c     to disc. ngnum = number of galaxies :
c     =====================================================
c
      subroutine outdata (ngnum)
c      
      include 'super.cb'
c      
      integer   ngnum,ipos,ihlen,i,j
      real      hoh(hlen),fpos,fhlen
      logical   ex
      character lfname*100
c      
      do i = 1,100
         lfname(i:i) = ' '
      enddo
c
      call makeext (fname,lfname,'.HEAD')
c
c     Name of this output-data-file is xxx.HEAD :
c     ===========================================
c
      inquire (file=lfname,exist=ex)
c      
      if (ex .eqv. .true.) then ! file already exists
         open (1,file=lfname,access='direct',recl=lori,status='old')
         read (1,rec = 1) fpos
         read (1,rec = 2) fhlen
         ipos  = int(fpos)
         ihlen = int(fhlen)
         close (1)
         open (1,file=lfname,access='direct',recl=ihlen*lori,
     $        status='old')
      else                      ! create a new file
         ihlen = hlen
         open (1,file=lfname,access='direct',recl=ihlen*lori,
     $        status='new')
         ipos  = 1
         ihlen = hlen
         write (1,rec = 1) hoh
      endif
c      
      do i = 1,ngnum
         do j = 1,ihlen
            hoh(j) = outhead(j,i)
         enddo
         write (1,rec = ipos+i) hoh
      enddo
c      
      hoh(1) = real(ipos+ngnum)
      hoh(2) = real(ihlen)
      hoh(3) = real(ngnum)
c
      write (1,rec = 1) hoh
c
      close (1)
c
      return
c
      end

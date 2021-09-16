c     
c******************************************************************
c
c     Write data (x,y,z,u,v,w) of a set of stars to disc:
c     ===================================================
c
c     step = current integrationstep, ign = which galaxy
c
      subroutine wsdata (step,ign)
c      
      include 'super.cb'
c      
      real      scv1,scl1,hstar(6)
      integer   step,ign,istart,istop,irec,gstart,i,ierr,j
      character lfname*100
c          
c     Creating name of star-data-file :
c     =================================
c
      do i = 1,100
         lfname(i:i) = ' '
      enddo
c
      call makename (fname,lfname,step,ign,ierr)
c      
      open (1,file=lfname,access='DIRECT',recl=17*lori)!!!!!careful with number of columns
c      
c     Finding first star of galaxy ign :
c     ==================================
c
      gstart = 0
c      
      if (ign.gt.1) then
         do i = 1,ign-1
            gstart = gstart + gstno(i)
         enddo
      endif
c      
c     Writing data of the desired stars in the output file :
c     ======================================================
c
      istart = tstart(ign)      ! int (fh(22,ign))
      istop  = tstop(ign)       ! int (fh(23,ign))
      irec = 1
c
      do i = gstart+istart,gstart+istop
c
         do j = 1,3
            hstar(j) = star(j,i)/scl
         enddo
         do j = 4,6
            hstar(j) = star(j,i)/scv
         enddo
         write (1,rec = irec) hstar(1),hstar(2),hstar(3),hstar(4),
     &        hstar(5),hstar(6),
     &        fh(14,ign),fh(15,ign),fh(16,ign),
     &        fh(17,ign),fh(18,ign),fh(19,ign),
     &        de(1,i),de(2,i),tm(i), de(4,i),dv2(i)

         irec = irec + 1
c
      enddo
c
      close (1)
c
      return
c
      end

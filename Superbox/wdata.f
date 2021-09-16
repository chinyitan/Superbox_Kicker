c      
c******************************************************************
c
c     Write all data to disc :
c     ========================
c
      subroutine wdata
c      
      include 'super.cb'
c      
      real      hfh(60),hstar(6)
      integer   i,j,iso,ioff
      character lfname*100
c      
      do i = 1,100
         lfname(i:i) = ' '
      enddo
c
      if (model .eq. 0) call convert (1) ! Convert back to physical units.. 

      call makeext (fname,lfname,'.CONT')
c      
      open (1,file=lfname,access='direct',recl=60*lori,iostat=iso)
c      
c     First, write INTEGER-header :
c     =============================
c
      write (1,rec = 1) ih
c
c     Next, write REAL-header for every galaxy :
c     ==========================================
c
      do i = 1,ih(1)            ! ih(1) = number of galaxies
         do j = 1,60
            hfh(j) = fh(j,i)
         enddo
         write (1,rec = i+1) hfh
      enddo
c      
      close (1)
c      
c     Reopen file with new recl :
c     ===========================
c
      open (1,file=lfname,access='direct',recl=6*lori,iostat=iso)
c
c     Write data of all stars :
c     =========================
c
      ioff = 10 + ih(1)*10        ! skip header data
c      
      do i = 1,ih(2)            ! ih(2) = TOTAL number of stars
         do j = 1,6
            hstar(j) = star(j,i)
         enddo
         write (1,rec = ioff+i) hstar
      enddo
c      
      close (1)
c      
c     if you wish to have additional backup files namely
c     xxx.CONT0 and xxx.CONT1 you should set the backup-flag
c     backup = 1
c      
      if (backup.eq.0) goto 31
c
      do i = 1,100
         lfname(i:i) = ' '
      enddo
c
c     The name of the additional output-file switches between CONT0
c     and CONT1, so if a system breakdown happens while writing data
c     you still have an old xxx.CONTx file
c
      if (output.eq.0) then
         call makeext (fname,lfname,'.CONT0')
         output = 1
      else
         call makeext (fname,lfname,'.CONT1')
         output = 0
      endif
c      
      open (1,file=lfname,access='direct',recl=60*lori,iostat=iso)
      
c     First, write INTEGER-header :
c     =============================
c
      write (1,rec = 1) ih
c
c     Next, write REAL-header for every galaxy :
c     ==========================================
c
      do i = 1,ih(1)            ! ih(1) = number of galaxies
         do j = 1,60
            hfh(j) = fh(j,i)
         enddo
         write (1,rec = i+1) hfh
      enddo
c      
      close (1)
c      
c     Reopen file with new recl :
c     ===========================
c
      open (1,file=lfname,access='direct',recl=6*lori,iostat=iso)
c
c     Write data of all stars :
c     =========================
c
      ioff = 10 + ih(1)*10        ! skip header data
c      
      do i = 1,ih(2)            ! ih(2) = TOTAL number of stars
         do j = 1,6
            hstar(j) = star(j,i)
         enddo
         write (1,rec = ioff+i) hstar
      enddo
c      
      close (1)
c      
 31   continue
c
      if (model .eq. 0) call convert (0) ! .. convert back to model units 
      return
c
      end

c     this program reads the *.HEAD file
c     
      program rheader
c     
c     h( 1)        stepnumber
c     h( 2)        length of timestep
c     h( 3)        time of simulation
c     h( 4)..h( 9) center of mass (galaxy)
c     h(10)..h(15) center of density (galaxy)
c     h(16)..h(21) center of mass (whole system)
c     h(22)        number of particles left (galaxy)
c     h(23)        kinetic energy (galaxy)
c     h(24)        potential energy (galaxy)
c     h(25)        kinetic energy (of whole system)
c     h(26)        potential energy (of whole system)
c     h(27)..h(29) angular momentum (galaxy), in rectilinear reference frame 
c                  centered on the centre of mass of the galaxy 
c     h(30)..h(32) angular momentum in rectilinear reference frame
c                  centered on the centre of mass of the whole system
c     h(33)..h(41) Lagrange-radii of an individual galaxy 
c                  containing 10%..(10)..90% of the mass
c
      integer       galnum
      parameter     (galnum=20)   ! max.number of galaxies

      real          h(60)
      character*100 fname,lfname,fname1,lfname1
      logical       ex
      integer       wg,count(galnum)
c
      do i = 1,100
         fname(i:i)   = ' '
         lfname(i:i)  = ' '
         fname1(i:i)  = ' '
         lfname1(i:i) = ' '
      enddo
c
      write (6,'(a)') ' '
      write (*,*) ' program to read the HEAD - file '
      write (6,'(a)') ' '
c
      write (*,FMT='($,a)') 'input filename without extension : '
      read (*,'(a)') fname
c
      call makeext (fname,lfname,'.HEAD')
c     
      inquire (file=lfname,exist=ex)
      if (ex .eqv. .false.) then
         write (*,'(a)') 'file does not exist...'
         stop
      endif
c     
 42   do i = 1,20
         count(i) = 0
      enddo  
c     
      open (1,FILE=lfname,ACCESS='DIRECT',RECL=4,IOSTAT=iso,
     &     STATUS = 'OLD')
c
      read (1,rec = 1) fpos
      read (1,rec = 2) fhlen
      read (1,rec = 3) fgnum
      ipos  = int(fpos)
      ihlen = int(fhlen)
      ignum = int(fgnum)
c
      write (6,'(a)') ' '
      write (6,'(a,i6,a)') 'length of record   : ',ihlen,' floats'
      write (6,'(a,i6)') 'number of records  : ',(ipos-1)/ignum
c
      close (1)
c
      open (1,FILE=lfname,ACCESS='DIRECT',RECL=ihlen*4,IOSTAT=iso,
     &     STATUS = 'OLD')
c     
      write (*,*) '  '
c     
      write (6,'(a,i6)') 'number of galaxies : ',ignum
c
      if (ignum .gt. 1) then
         write (*,'(a)') ' '
 1       write (6,'($,a)') 'which galaxy (0 = all) : '
         read (5,*) wg
         if ((wg .lt. 0) .or. (wg .gt. ignum)) then
            write (6,'(a)')
     $           'your input is out of range... try iy again'
            goto 1
         else
            goto 2
         endif  
      else
         wg = 1
      endif
c 
 2    continue       
c
      write (*,*) ' '
c     
      write (6,'(a)') ' '
      write (6,'(a)') ' 1) energy '
      write (6,'(a)') ' 2) Lagrange radii (10%, 20% ..... 90%) '
      write (6,'(a)') ' 3) angular momentum '
      write (6,'(a)') ' 4) center of mass '
      write (6,'(a)') ' 5) center of density '
      write (6,'(a)') ' 6) center of simulation '
      write (6,'(a)') ' 7) number of particles left'
      write (6,'(a)') '99) exit program'
      write (6,'(a)') ' '
      write (6,'($,a)') 'expecting your command : '
      read (5,*) inum
c     
c     start input loop        
c     
 5    if ((inum .gt. 0) .and. (inum .lt. 8)) goto 10
c
      if (inum .eq. 99) stop
c
      write (6,'(a)') 'input is out of range...try it again'
      write (6,'($,a)') 'expecting your command : '
      read (5,*) inum
c
      goto 5    
c
 10   continue
c
      write (6,'(a)') ' '
c
c     stop input loop 
c      
c     filenames for output files 
c     energy for galaxy 1 : ene --> ene.g01
c      
      if (inum .eq. 1) fname1 = 'energy'
      if (inum .eq. 2) fname1 = 'lagr'
      if (inum .eq. 3) fname1 = 'angmom'
      if (inum .eq. 4) fname1 = 'xyz_m'
      if (inum .eq. 5) fname1 = 'xyz_d'
      if (inum .eq. 6) fname1 = 'tcms'
      if (inum .eq. 7) fname1 = 'rem'
c
      do k = 1,ignum
c
         count(k) = 0
c
         if ((wg .eq. 0) .or. (wg .eq. k)) then
c
            call makename (fname1,lfname1,0,k,ierr,ilen)
c
            open (2,FILE=lfname1)
c
            write (6,'($,a,a,a)') 
c
     &           'writing data in file ',lfname1(1:ilen),'.... '
c
            do i = 1,(ipos/ignum)-1
c
               ind = (i-1)*ignum+k+1
               read (1,rec = ind) h
c
               time = h(3)
c
c     energy : 
c     time [Myr], ekin of galaxy, epot of galaxy, etot of galaxy,
c     ekin of system, epot of system, etot of system; 
c     energies are in model units

               if (inum.eq.1) then 
                  write (2,'(I7,e11.3,6(1x,1p,e14.4))')
     $                INT(h(1)),time,
     $                h(23),h(24),h(23)+h(24),
     $                h(25),h(26),h(25)+h(26)
               endif

c     rad :
c     time and Lagrange radii 
c             [Myr]        [kpc]

               if (inum.eq.2) then 
                  write (2,'(I7,e11.3,9f12.4)') 
     $                 INT(h(1)),time,h(33),h(34),h(35),
     $                 h(36),h(37),h(38),h(39),h(40),h(41)
               endif
c
c     angmom :
c     time and angular momentum

               if (inum.eq.3) then 
c                  lg = sqrt(h(27)**2+h(28)**2+h(29)**2)
c                  ls = sqrt(h(30)**2+h(31)**2+h(32)**2)
c                  print*, lg,ls
                  write (2,'(I7,e11.3,1x,8f12.6)') INT(h(1)),time,h(27),
     &                 h(28),h(29),
     $                 sqrt(h(27)**2+h(28)**2+h(29)**2),h(30),h(31),
     &                 h(32),sqrt(h(30)**2+h(31)**2+h(32)**2)
               endif

c     xyz_m :
c     time and position and velocity of center of mass
c     [Myr]     [kpc]        [km/s]
c
               if (inum.eq.4) then
                  write(2,'(I7,e11.3,6f14.6)') 
     $                 INT(h(1)),time,h(04),h(05),h(06),h(07),
     $                 h(08),h(09)
               endif

c     xyz_d :
c     time [Myr] and position [kpc] and velocity [km/s] 
c     of center of density

               if (inum.eq.5) then
                  write(2,'(I7,e11.3,6f14.6)') 
     $                 INT(h(1)),time,h(10),h(11),h(12),h(13),
     $                 h(14),h(15)
               endif

c     tcms :
c     time [Myr], position [kpc] and velocity [km/s] of centre of mass 
c     of the whole system 
c
               if (inum.eq.6) then
                  write(2,'(e11.3,1x,6f12.6)') time,h(16),h(17),h(18),
     &                                         h(19), h(20),h(21)
               endif
c
c     rem :
c     time [Myr] and particles of a galaxy left in the computation
               if (inum .eq. 7) then
                  write (2,'(I7,e11.3,i10)')INT(h(1)),time,int(h(22))
               endif
c
               count(k) = count(k) + 1
c
            enddo
c
            close (2)
c
            write (6,'(a)') 'done.'
c
         endif
c
      enddo
c
      close (1)
c
      write (6,'(a)') ' '
c
      if (wg .eq. 0) then
         do i = 1,ignum
            write (6,'(a,i2,a,i6)')
     $           'number of records in galaxy ',i,' ',
     &           count(i)
         enddo 
      else
         write (6,'(a,i2,a,i6)')
     $        'number of records in galaxy ',wg,' ',
     &        count(wg)
      endif
c
      write (6,'(a)') ' '
c
      goto 42
c
      end
c     
c******************************************************************
c
c     concanate 'fname' with 'ext'
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
c
      enddo
c
      f2(i-1+j:i-1+j) = ' '
c
 20   continue
c
      end
c     
c******************************************************************
c
c     creating output file name
c
      subroutine makename (fname,lfname,step,ign,ierr,ilen)
c     
      character*100 fname,lfname
      integer       ind(5),step,ign,izero,istep,ierr,i,ilen
c     
      ind(1) = 10000
      ind(2) =  1000
      ind(3) =   100
      ind(4) =    10
      ind(5) =     1
      izero  =     0
      istep  = step
c     
      do i = 1,100
         lfname(i:i) = ' '
      enddo
c
      do i = 1,100
         if (fname(i:i) .eq. ' ') goto 10
         lfname(i:i) = fname(i:i)
      enddo
c
      ierr = 1
c
      goto 11
c
 10   lfname(i:i) = '.'
      lfname(i+1:i+1) = 'g'
c
      if (ign .lt. 10) then
         write (lfname(i+2:i+2),'(i1)') izero
         write (lfname(i+3:i+3),'(i1)') ign
      else
         write (lfname(i+2:i+3),'(i2)') ign
      endif
c
      ierr = 0
      ilen = i+3
c
 11   continue
c     
      end
c
c******************************************************************
c     end of file
c******************************************************************










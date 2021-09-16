       program itens
*
C Reads binary data from file.ITENS produced by superbox containing 
C rotational inertia data of a particular galaxy, and writes 
C here selected data (diagonal elements of tensor) to ASCII
C data file (itensdat.jnk), and if desired, the whole tensor and 
C eingenvalues to file itens_eigvect.jnk.
*
C FORMAT:
C FIRST RECORD: IPOS  = NUMBER OF RECORDS
C               IHLEN = NUMBER OF REAL'S
C
C TIME
C STEP
C RMAX
C NUMBER OF PARTICLES
*
*
*
C Time is stored in                   h(1)
C Step is in                          h(2)
C Number of stars with r<=rmax is in  h(3)
C rmax is in                          h(4)
C eigenvalues 1,2,3 are in            h(6,7,8)
C Eigenvectors are stored in h(9,12,15) for eigenvalue1
C                            h(10,13,16) for eigenvalue2
C                            h(11,14,17) for eigenvalue3
C See superbox routine "eval ": the relevant data are stored in 
C array hoh(k) such that here h(k) = hoh(k) 
C
C--------------checked it all again on 9.10.96!!!
*
*
*
        
       integer intstep, model, unitl, ilen
       parameter( ilen = 23 ) 
       real*4 time, ih5, ih13, kmspcy, scl,sct,scm,scv,gp  
       real h(ilen)
       character*100 fname,lfname,lfname1

       character    rpc*4, munit*5, tunit*5

      common /scales/ scl,scm,sct,scv, gp,kmspcy,
     &                model,unitl,tunit,munit,rpc  
*
       do i=1,100
          fname(i:i)=' '
       end do
       write (*,FMT='($,a)') 'name of model w/o extension: '
       read (*,'(a)')  fname
       write (*,FMT='($,a)') 'which galaxy number ? '
       read (*,'(i2)') ignum
       call createna (fname,lfname,ignum,ierr)
       call makeext  (lfname,lfname1,'ITENS')
       write(6,*)
       write(6,'(2a)')' Accessing file ',lfname1
*
       write(6,*)
       write(6,'($,a,a)')' write itens_eigvect.jnk for Nth ',
     &  'step ? (N<0 => not written) '
       read(5,*) intstep
       write(6,'(I6)') intstep
*
c open file with record length 4 (correct for REAL*4)
       open (1,FILE=lfname1,ACCESS='DIRECT',RECL=4,IOSTAT=iso,
     &         STATUS = 'OLD')
       read (1,rec = 1) fpos
       read (1,rec = 2) fhlen
       read (1,rec = 3) ih5   ! choice of model or physical units
       read (1,rec = 4) gp    ! Gravitational constant 
       read (1,rec = 6) scm
       read (1,rec = 7) scl
       read (1,rec = 8) sct
       read (1,rec = 9) scv
       read (1,rec = 10) kmspcy 
       read (1,rec = 11) ih13
       model = int(ih5) 
       unitl = int(ih13) 
       ipos  = int(fpos)
       ihlen = int(fhlen) + 6 
       close (1)
       write( *,* ) 'model and unit length = ', ih5, ih13 

*       Set scales if input/output not in model units - only cosmetic

       call set_scales 

* now open file with correct record length

       open (1,FILE=lfname1,ACCESS='DIRECT',RECL=ihlen*4,IOSTAT=iso,
     &         STATUS = 'OLD')

c+++
       print*, 'number of records: ',ipos-1,' hlen = ',ihlen

       open (unit=12,file='itensdat.jnk')
       do i = 1,ipos-1
        read (1,rec = i+1) h
        time = h(1)
*
*
C Write out: time  step rmax  count eigenvalues1,2,3 [Msun * kpc^2], axes, angles
        write (12,'(F8.4,F7.0,1x,F6.3,1x,F8.0,1x,9(x,E10.3))') 
     +  time,h(2),h(3),h(4),h(6),h(7),h(8),h(9),h(10),h(11),
     +   h(12),h(13),h(14)
c        write (*,'(F10.3,F8.0,F9.3,F10.0,3(1p,E9.3))') 
c     +           time,h(2),h(3),h(4),h(6),h(7),h(8)

*
C Write into file "eigenvect.jnk" the eigenvector:
        if (intstep.EQ.h(2)) then
          open(13,file='itens_eigenvect.jnk')
          write(13,'(a,F10.3,a,a,a,F8.0)')
     +           'time: ',time,' [ ',tunit,']; step: ',h(2)
          write(13,'(a,a,1x,a,a,3(F12.3,1x))')
     +    'LOG10( eigenvalues [ ',munit,rpc,'^2 ]): ',
     +              LOG10(h(6)),LOG10(h(7)),LOG10(h(8))
          write(13,'(a,3(F12.3,1x))')
     +   'x-comps:             ', h(15),h(16),h(17)
          write(13,'(a,3(F12.3,1x))')
     +   'y-comps:             ', h(18),h(19),h(20)
          write(13,'(a,3(F12.3,1x))')
     +   'z-comps:             ', h(21),h(22),h(23)
           close (13)
           write(6,*)
           write(6,*)' "itens_eigenvect.jnk" written'
           write(6,*)
        end if
       end do
*
*
*
       close (1)
       close (12)
*
       write(6,*)
       write(6,*)' "itensdat.jnk" written'
       write(6,*)
*
*
*
       end

C--------------------------------------------------------
C concanate 'fname' with 'ext'
C--------------------------------------------------------

      SUBROUTINE makeext (f1,f2,ext)

      character*(*) f1,f2,ext

      do i=1,100
         f2(i:i)=' '
      end do
      do i = 1,len(f1)
        if (f1(i:i) .eq. ' ') goto 10
        f2(i:i) = f1(i:i)
      end do
  10  continue
      do j = 1,len(ext)
         if (ext(j:j) .eq. ' ') goto 20
         f2(i-1+j:i-1+j) = ext(j:j)
      end do
      f2(i-1+j:i-1+j) = ' '
  20  continue
      end


       SUBROUTINE createna (fname,lfname,ign,ierr)

       character*100 fname,lfname
       integer       ign,izero,istep,ierr,i,j,ir

       izero  =    0
       istep  = step

        do i=1,100
           lfname(i:i)=' '
        end do
        do i = 1,100
         if (fname(i:i) .eq. ' ') goto 10
         lfname(i:i) = fname(i:i)
        end do
        ierr = 1
        goto 11
  10    lfname(i:i) = '-'
        lfname(i+1:i+1) = 'g'
        if (ign .lt. 10) then
         write (lfname(i+2:i+2),'(i1)') izero
         write (lfname(i+3:i+3),'(i1)') ign
        else
         write (lfname(i+2:i+3),'(i2)') ign
        end if
        lfname(i+4:i+4) = '.'
        lfname(i+5:i+5) = ' '
        ierr = 0
  11    continue

        end


      subroutine set_scales 

c       In model units for SUPERBOX we have by default
c       Gs = 2 ; Ms = 1 = Rs   (s = superbox)  
c       I write the physical quantities G, M, R and t such that 
c       Gs = a x G  ; Ms = b x M ; Rs = c x R and ts = d x t 
c
c       The physical unit of mass is the solar mass, time is in years
c       and length has been defined above. Then by construction 
c         
c       [ts]^-2 = GsMs/Rs^3 = ab/c^3 [t]^-2 => [ts] = sqrt(c^3/(ab)) [t]
c       [Vs]^2  = GsMs/Rs   = ab/c GM/R   => [Vs] = sqrt(ab/c) [V] 
c                                                 = cv [V] 
c       Also 
c       [Ms] = ts^-2Rs^3/Gs = d^-2c^3/a  t^-2R/G => [Ms] = d^-2c^3/a [M]
c       and of course d = sqrt(c^3/(ab)).
c 
c      fh(54) = gp 
c      fh(55) = a 
c      fh(56) = b 
c      fh(57) = c
c      fh(58) = d 
c      fh(59) = cv = c/d
c      fh(60) = kmspcy 

      integer   galnum
      parameter (galnum=2)  ! max.number ofgalaxies

      integer      ih(60),gnum, unitl,cnum, model  
      real         fh(60,galnum),sl,ml,tl,vl
      real         a,b,c,cv,d,gp,kmspcy 

      character    rpc*4, munit*5, tunit*5

      common /scales/ sl,ml,tl,vl, gp,kmspcy,
     &                model,unitl,tunit,munit,rpc  

      tunit = 'Myrs '
      munit = 'solar'
      rpc = 'pc' 
      if( unitl .gt. 0.5 ) rpc = 'Kpc' 
      if( model .gt. 0.5 ) rpc = 'm.u.'
      if( rpc .eq. 'm.u.' ) tunit = rpc 
      if( rpc .eq. 'm.u.' ) munit = rpc 
 
      if( rpc .eq. 'm.u.' ) then 
        sl = 1. 
        ml = 1. 
        tl = 1.
        vl = 1. 
      endif 

      return 

      end 

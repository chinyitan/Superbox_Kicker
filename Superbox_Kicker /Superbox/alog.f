      character*100    fname,lfname,outname
      logical          ex
      real     atot,afirst,alog,tlog,time,apro
      integer  ipos,ihlen,ignum,ind
      real     h(60),fpos,fhlen,fgnum
      do i=1,100
         lfname(i:i) = ' '
         outname(i:i) = ' '
      enddo
      write (*,FMT='($,a)') 'input filename without extension : '
      read (*,'(a)') fname
      call makeext (fname,lfname,'.HEAD')
      call makeext (fname,outname,'.alog')
      inquire (file=lfname,exist=ex)
      if (ex .eqv. .false.) then
         write (*,'(a,a,a)') 'file ',lfname,' does not exist...'
         stop
      endif
      open (1,FILE=lfname,ACCESS='DIRECT',RECL=4,IOSTAT=iso,
     &     STATUS = 'OLD')
      read (1,rec = 1) fpos
      read (1,rec = 2) fhlen
      read (1,rec = 3) fgnum
      ipos  = int(fpos)
      ihlen = int(fhlen)
      ignum = int(fgnum)
      close (1)
      open (1,FILE=lfname,ACCESS='DIRECT',RECL=ihlen*4,IOSTAT=iso,
     &         STATUS = 'OLD')
      open (2,file=outname)
      do i = 1,(ipos/ignum)-1
         ind = (i-1)*ignum+2
         read (1,rec = ind) h
         if (i .eq. 1) then
            afirst = sqrt(h(30)**2 + h(31)**2 + h(32)**2)
            goto 10
         endif
         time = h(1)
         tlog = log10(time)
         atot = sqrt(h(30)**2 + h(31)**2 + h(32)**2)
         apro = (afirst - atot) / afirst
         alog = log10(abs(apro))
         write (2,'(2f12.8)') tlog,alog
 10      continue
      enddo
      close (2)
      close (1)
      end

      SUBROUTINE makeext (f1,f2,ext)

      character*(*) f1,f2,ext

      do i = 1,100
         f2(i:i) = ' '
      enddo
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


c
c******************************************************************
c
      subroutine out_time (step,cstep,time,unit)
c
c     note: this subroutine only calculates how much CPU-time 
c           the simulation has needed ; it's not essential for
c           the simulation, so if it doesn't work on your machine
c           simply skip it
c
      integer step,cstep, unit 
      real    time
c
      fdiv  = real(step) - real(cstep)
      xtime = time
      ih1   = int(time / 3600.0)
      time  = time - ih1*3600.0
      im    = int(time / 60)
      time  = time - im*60
      is    = int(time)
      ihs   = int((time - is) * 100.0)
c
      write (unit,'(a)') ' '
      write (unit,'(a)') '**************************
     $     *************************'
      write (unit,'(a,i3,a,i2,a,i2,a,i2,a)') ' CPU-time      :  ',
     &  IH1,' h  ',IM,' m  ',IS,' s  '
      write (unit,'(a,f7.3)') ' sec per step  :  ',xtime/fdiv
      write (unit,'(a)') '**************************
     $     *************************'
      write (unit,'(a)') ' '
c
      return
c
      end

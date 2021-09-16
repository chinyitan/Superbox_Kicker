c
c******************************************************************
c
      subroutine out_error (num)
c
c     this routine simply prints the error-messages if something
c     is wrong with the input data
c
      include 'super.cb'
c
      integer  num
c
      if (num .eq. 100) then
         write (*,'(a)') 'Missing command line argument'
         write (*,'(a)') 'Program terminated !!'
      endif
c
      if (num.eq.101) then
         write (unit,'(a)') 'The parameter MSTEP is smaller or equal'
         write (unit,'(a)') 'than the current integration step !!'
         write (unit,'(a)') 'please change MSTEP '
         write (unit,'(a)') 
     $        'Use <define> to set the parameters right !'
         write (unit,'(a)') 'Program terminated !!'
      endif
c
      if (num.eq.102) then
         write (unit,'(a)') 'Timestep is less or equal 0 !'
         write (unit,'(a)')
     $        'Use <define> to set the parameters right !'
         write (unit,'(a)') 'Program terminated !!'
      endif
      close (unit)
      stop
c
      end

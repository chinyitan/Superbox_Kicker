c      
c******************************************************************
c
c     If a file named 'stop' exist in your HOME-directory,
c     stop execuction
c
      subroutine dostop (flag)
c
      integer flag
      logical ex
c      
      inquire (file='stop',exist=ex)
c      
      if (ex .eqv. .true.) then ! user want to stop execuction
        flag = 1
      else                      ! goon
         flag = 0
      endif
c
      return
c
      end

c**********************************************************************
c
c     I O F U N C . F
c
c**********************************************************************
c
c     Read data for the simulation from xxx.CONT file :
c     =================================================
c
      subroutine rdata
c      
      include 'super.cb'
c      
      real*4    hfh(60),hstar(6), r, izz, count  
      integer   i,j,iso,ioff,i_com,iih(60)
      character lfname*100
      logical   ex
c
c     fname  = name of model without extension (fname <-- super.cb)
c     lfname = name of model with extension '.CONT'
c
      do i = 1,100
         lfname(i:i) = ' '
      enddo
c
      call makeext (fname,lfname,'.CONT')
c     
      inquire (file=lfname,exist=ex)
c      
      if (ex .eqv. .false.) then
         write (99,'(a,a,a)') 'file ',lfname,' does not exist...'
         write (99,'(a)') 'program terminated !!'
         stop
      endif
c      
      open (1,file=lfname,access='direct',recl=60*lori,iostat=iso)
c      
c     First, read  Integer-Header :
c     =============================
c
      read (1,rec = 1) ih       ! ih defined in super.cb
 
            
c+++
c      print*,' iofunc.f: ih finished'

c     Next, read REAL-Headers of all galaxies :
c     =========================================
c     
      do i = 1,ih(1)            ! ih(1) = number of galaxies 
         read (1,rec = 1+i) hfh
 
         do j = 1,60
            fh (j,i) = hfh(j)   ! fh defined in super.cb            
         enddo
      enddo
      
c+++
c      print*,' iofunc.f: fh finished'
      
      close (1)
     
c
c     Reopen file with recl = 6*real :
c     ================================
c
c     x,y,z,u,v,w of a star
c
      open (1,file=lfname,access='direct',recl=6*lori,iostat=iso)
     
c
c     Read data of stars :
c     ====================
c
      ioff = 10 + 10*ih(1)        ! skip header data
     
c     
      do i = 1,ih(2)              ! ih(2) = TOTAL number of stars
         read (1,rec = ioff+i) hstar         
         do j = 1,6
            star(j,i) = hstar(j)
         enddo
         if (mod(i,100000).eq.0) print*,i,' iofunc.f: stars finished'
      enddo
c
      close (1)
c
c++     Check input file : debugging incompatible binaries .. (cmb.26.03.99)
c
c     open( unit = 32, file = 'supercheck.d', status = 'unknown' )

      count = 0. 
      izz   = 0. 
      do i = 1,ih(2)              ! ih(2) = TOTAL number of stars

      r = star(1,i)**2 + star(2,i)**2  + star(3,i)**2  
      r = sqrt(r) 

      count = count + 1 
      izz = izz + r*r-star(3,i)*star(3,i)
c      write( 32,'(i6,1x,1e12.7,6(1x,f10.6), 1x,1e12.7)') 
c     &    i, count,( star(j,i), j=1,6), r
     
      end do 

      write( 6,* ) 'Total number of input stars,Izz = ', 
     &   ih(2),count,izz/count  
c++ 
      return
c
      end

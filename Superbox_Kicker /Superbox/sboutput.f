      implicit real (a-h,o-z)
c      real x(80000),y(80000),z(80000)
c      real vx(80000),vy(80000),vz(80000)
      real fh(60),x,y,z,u,v,w
      integer ih(60)
      integer nbodies,ngalaxies

      ngalaxies=1
      nbodies=10000


      open (2,file='build40000')
      open (1,file='build4.STARTUP',access='direct',status='old',
     $     recl=240)
      
      read (1,rec=1) ih
      do i=1,ngalaxies
         read(1,rec=1+i) fh
      enddo
      
      close(1)
      
      open (1,file='build4.STARTUP',access='direct',status='old',
     $     recl=24)
      
      do i=1,nbodies
         read(1,rec=10+ngalaxies*10+i) x,y,z,u,v,w
         write(2,'(6f20.10)')  x,y,z,u,v,w
      enddo 

      close(1)
      close(2)

      end

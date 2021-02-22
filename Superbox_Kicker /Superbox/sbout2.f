      implicit real (a-h,o-z)
      real x,y,z,u,v,w
      integer nbodies,ngalaxies

      ngalaxies=1
      nbodies=10000
      
      open (2,file='build44500')
      open (1,file='build4-g01.04500',access='direct',status='old',
     &     recl=24)
      
      do i=1,nbodies
         read(1,rec=i) x,y,z,u,v,w
         write(2,'(6f20.10)')  x,y,z,u,v,w
      enddo
      
      close(1)
      close(2)
      
      end

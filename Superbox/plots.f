      implicit real (a-h,o-z)
      real r(10,5)
      integer instep,in

      instep=5
      
      open(10,file='lagr.g01',status='old')
      do i=1,instep
         read(10,*) in,r(1,i),r(2,i),r(3,i),r(4,i),r(5,i),
     &        r(6,i),r(7,i),r(8,i),r(9,i),r(10,i)
         write(6,*) in,r(1,i),r(2,i),r(3,i),r(4,i),r(5,i),
     &        r(6,i),r(7,i),r(8,i),r(9,i),r(10,i)
      enddo
      close(10)

      open(11,file='plotsm.dat')
      do j=2,10
         write(11,*) 0.1*(j-1),r(j,1),r(j,2),r(j,3),r(j,4),r(j,5)
      enddo
      close(11)
      end


c     nota:vease que solo vamos del 10% al 90%

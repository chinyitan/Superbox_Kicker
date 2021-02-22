      implicit real(a-h,o-z)
      real x(80000),y(80000),z(80000)
      real vx(80000),vy(80000),vz(80000)
      real mlx,mly,mlz,mpx,mpy,mpz
      
      nbodies=80000

      open (10,file='build5out.dat',status='old')
      do i=1,nbodies
         read(10,*) x(i),y(i),z(i),vx(i),vy(i),vz(i)
      enddo
      close (10)

       mlx=0.
       mly=0.
       mlz=0.
       mpx=0.
       mpy=0.
       mpz=0.
       
       do i=1,nbodies
           mpx=mpx+vx(i)
           mpy=mpy+vy(i)
           mpz=mpz+vz(i)
           mlx=mlx+y(i)*vz(i)-vy(i)*z(i)
           mly=mly+z(i)*vx(i)-x(i)*vz(i)
           mlz=mlz+x(i)*vy(i)-y(i)*vx(i)
            
       enddo

       write(6,*) ' '
       write(6,*) 'px',mpx/nbodies,'  py',mpy/nbodies,'  pz',mpz/nbodies
       write(6,*) 'lx',mlx/nbodies,'  ly',mly/nbodies,'  lz',mlz/nbodies 
       write(6,*) 'total nbodies in buildjorge.dat',nbulge*8
       write(6,*) ' '

       end

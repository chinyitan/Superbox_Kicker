      
      implicit real (a-h,o-z)
      real x(80000),y(80000),z(80000)
      real xxim2(10),yyim2(10),zzim2(10)
      real ba2(10),ca2(10),cb2(10),tau2(10)
      integer num,nbodies

c     calculus of the axis ratio of our system.
c     calculus of center of mass

      nbodies=10000

      open(10,file='build40000',status='old')

      do i=1,nbodies
         read(10,*,end=10) x(i),y(i),z(i) 
      enddo   
      
 10   continue

      close(10)

      open(11,file='inmom.dat')

      xxprod=0.
      yyprod=0.
      zzprod=0.
      xsum=0.
      ysum=0.
      zsum=0.

      do i=1,nbodies
c     inertia momments
         xxprod=x(i)*x(i)
         yyprod=y(i)*y(i)
         zzprod=z(i)*z(i)
         
         xxim=xxim+xxprod
         yyim=yyim+yyprod
         zzim=zzim+zzprod

c     center of mass
         xsum=xsum+x(i)
         ysum=ysum+y(i)
         zsum=zsum+z(i)
      
      enddo

      ba=(yyim/xxim)**0.5
      ca=(zzim/xxim)**0.5
      cb=(zzim/yyim)**0.5
      tau=(ba-ca)/(1.-ca)

      write(11,*) 'b/a total',ba
      write(11,*) 'c/a total',ca
      write(11,*) 'c/b total',cb
      write(11,*) 'tau total',tau
      write(11,*) ' '

      call sort (nbodies,x)
      call sort (nbodies,y)
      call sort (nbodies,z)

      do j=1,10
         num=j*0.1*nbodies

         xxprod=0.
         yyprod=0.
         zzprod=0.

         do i=1,num
            xxprod=x(i)*x(i)
            yyprod=y(i)*y(i)
            zzprod=z(i)*z(i)

            xxim2(j)=xxim2(j)+xxprod
            yyim2(j)=yyim2(j)+yyprod
            zzim2(j)=zzim2(j)+zzprod
         enddo

         ba2(j)=(yyim2(j)/xxim2(j))**0.5
         ca2(j)=(zzim2(j)/xxim2(j))**0.5
         cb2(j)=(zzim2(j)/yyim2(j))**0.5
         tau2(j)=(ba2(j)-ca2(j))/(1.-ca2(j))


         write(11,*) 'masa',j*10,'%'
         write(11,*) 'b/a',ba2(j)
         write(11,*) 'c/a',ca2(j)
         write(11,*) 'c/b',cb2(j)
         write(11,*) 'tau',tau2(j)
         write(11,*)'---------------'
      
      enddo

      write(11,*) ''
      write(11,*) 'CENTER OF MASS'
      write(11,*) 'Xcm',xsum/nbodies
      write(11,*) 'Ycm',ysum/nbodies
      write(11,*) 'Zcm',zsum/nbodies
      
      close(11)

      end






      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL*4 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*4 a,temp
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then
         do 12 j=l+1,ir
            a=arr(j)
           do 11 i=j-1,1,-1
              if(arr(i).le.a)goto 2
             arr(i+1)=arr(i)
 11       continue
         i=0
 2         arr(i+1)=a
 12     continue
       if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
           temp=arr(l+1)
           arr(l+1)=arr(ir)
           arr(ir)=temp
       endif
       if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
       endif
       i=l+1
       j=ir
       a=arr(l)
 3     continue
       i=i+1
       if(arr(i).lt.a)goto 3
 4     continue
       j=j-1
       if(arr(j).gt.a)goto 4
       if(j.lt.i)goto 5
       temp=arr(i)
       arr(i)=arr(j)
       arr(j)=temp
       goto 3
 5     arr(l)=arr(j)
       arr(j)=a
       jstack=jstack+2
       if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
      endif
      goto 1
      END

c      
c******************************************************************
c      
c     quicksort-algorithm from
c     Numerical Recipies 2nd edition
c      
      subroutine sort(n,arr)
c      
      integer n,m,nstack
      real    arr(n)
c      
      parameter (m = 7 , nstack = 50)
c      
      integer i,ir,j,jstack,k,l,istack(nstack)
      real    a,temp
c      
      jstack = 0
      l      = 1
      ir     = n
c      
 1    if (ir-l .lt. m) then
c         
         do 12 j = l+1,ir
c            
            a        = arr(j)
c            
            do 11 i = j-1,1,-1
c               
               if (arr(i) .le. a) goto 2
c               
               arr(i+1) = arr(i)
c               
 11         continue
c            
            i        = 0
 2          arr(i+1) = a
c            
 12      continue
c         
         if (jstack .eq. 0) return
c         
         ir       = istack(jstack)
         l        = istack(jstack-1)
         jstack   = jstack - 2
c         
      else
c         
         k        = (l+ir)/2
         temp     = arr(k)
         arr(k)   = arr(l+1)
         arr(l+1) = temp
c         
         if (arr(l+1) .gt. arr(ir)) then
            temp     = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir)  = temp
         endif
c         
         if (arr(l) .gt. arr(ir)) then
            temp     = arr(l)
            arr(l)   = arr(ir)
            arr(ir)  = temp
         endif
c         
         if (arr(l+1) .gt. arr(l)) then
            temp     = arr(l+1)
            arr(l+1) = arr(l)
            arr(l)   = temp
         endif
c         
         i = l+1
         j = ir
         a = arr(l)
c         
 3       continue
c        
         i = i+1
c         
         if (arr(i) .lt. a) goto 3
c         
 4       continue
c         
         j = j-1
c         
         if (arr(j) .gt. a) goto 4
         if (j .lt. i)      goto 5
c         
         temp   = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
c         
         goto 3
c         
 5       arr(l) = arr(j)
         arr(j) = a
         jstack = jstack+2
c         
         if (jstack .gt. NSTACK) pause 'NSTACK too small in sort'
c         
         if (ir-i+1 .ge. j-l) then
            istack(jstack)   = ir
            istack(jstack-1) = i
            ir               = j-1
         else
            istack(jstack)   = j-1
            istack(jstack-1) = l
            l                = i
         endif
c         
      endif
c
      goto 1
c      
      return
c      
      end

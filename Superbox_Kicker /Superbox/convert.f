
c
c******************************************************************
c
      subroutine convert (mode)
c
c     this subroutine converts all data from phsical units into
c     model units or vice versa
c
      include 'super.cb'
c
      real    scm1,scl1,scv1,sct1
      integer mode
      integer i,j

c
      if (mode.eq.0) then ! physical --> model
         scm1 = scm
         scl1 = scl
         sct1 = sct
         scv1 = scv
      else                ! model --> physical
         scm1 = 1.0 / scm
         scl1 = 1.0 / scl
         sct1 = 1.0 / sct
         scv1 = 1.0 / scv
      endif
c
      do i = 1,istno      ! transform data of all stars
         if (star(1,i).lt.1.e10) then
            star(1,i) = star(1,i) * scl1
            star(2,i) = star(2,i) * scl1
            star(3,i) = star(3,i) * scl1
            star(4,i) = star(4,i) * scv1
            star(5,i) = star(5,i) * scv1
            star(6,i) = star(6,i) * scv1
         endif
      enddo
c     
      do i = 1,gnum       ! start loop for all galaxies
c
         mtot(i)    = mtot(i)    * scm1    ! total mass of galaxy i
         mass(i)    = mass(i)    * scm1    ! mass of one star
c
         rcore(i)   = rcore(i)   * scl1    ! grid - lengths
         rout(i)    = rout(i)    * scl1
         rsystem(i) = rsystem(i) * scl1
c
         dt(i)      = dt(i)      * sct1    ! timestep of galaxy i
         ctime(i)   = ctime(i)   * sct1    ! current time
c
         do j = 1,3                        
            gcms(j,i)  = gcms(j,i)  * scl1 ! center of mass
            dgcms(j,i) = dgcms(j,i) * scl1 ! center of density
         enddo
         do j = 4,6
            gcms(j,i)  = gcms(j,i)  * scv1
            dgcms(j,i) = dgcms(j,i) * scv1
         enddo
c
      enddo
c      
      return
c
      end

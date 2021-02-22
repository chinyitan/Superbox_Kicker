c      
c******************************************************************
c
c     Calculate velocity dispersion tensor of all stars within
c     radius rmax
c     ========================================================
c
c     ign = which galaxy, rmax = maximale distance for a star from
c     the center of density
c
      subroutine disptens (step,ign,rmax,time,iflag)
c
      include 'super.cb'
c      
      integer i,j,nrot,ierr, indx(3), shape_iter 
      integer ign,start,istop,step,ihlen,ipos,iflag
      real    tt(3,3),d(3),v(3,3),v_in(3,3),rmax,count,time
      real    hoh(hlend),fpos,fhlen, mstar, vxm, vym, vzm 
      real    a,b,c,dd,vprod, a1,b1,cc1, delta, vol
      double precision dtt(3,3),vcmsx,vcmsy,vcmsz,scl1,scv1
      double precision cmsx,cmsy,cmsz,vx,vy,vz,r,dx,dy,dz,rmax2
      character*100    lfname,lfname1
      logical          ex
c      
      count = 0.d0
c
      do i = 1,3
         do j = 1,3
            dtt(i,j) = 0.d0
         enddo
      enddo
c      
      if (ign .eq. 1) then
         start  = 1
         istop  = gstno(1)
      else
         start  = 0
         do i = 1,ign-1
            start = start + gstno(i)
         enddo
         istop  = start + gstno(ign)
         start  = start + 1
      endif
c
      cmsx  = gcms(1,ign) ! fh(14,ign))  
      cmsy  = gcms(2,ign) ! center of density or of mass
      cmsz  = gcms(3,ign)   
      vcmsx = gcms(4,ign) 
      vcmsy = gcms(5,ign) 
      vcmsz = gcms(6,ign) ! fh(19,ign))

c++     Add small value to maximum radius to fight off roundoff errors.
      rmax2 = dble(rmax)*dble(rmax) + 2.*epsilon*rmax 
c     
c++     First compute the stream velocity, ie mean values in each component. 

      count = 0. 

      vxm = 0.
      vym = 0. 
      vzm = 0. 

      do i = start,istop
        
         if (abs(star(1,i)) .lt. 1.e10 .and. abs(star(2,i)) .lt. 1.e10 
     &       .and. abs(star(3,i)) .lt. 1.e10 ) then

            dx = star(1,i) - cmsx
            dy = star(2,i) - cmsy
            dz = star(3,i) - cmsz
            r  = dx*dx + dy*dy +dz*dz
c            
            if (r .lt. rmax2) then
               count = count + 1.d0

               vx = star(4,i) - vcmsx
               vy = star(5,i) - vcmsy
               vz = star(6,i) - vcmsz
c               
               vxm = vxm + vx
               vym = vym + vy
               vzm = vzm + vz

               dtt(1,1) = dtt(1,1) + vx*vx
               dtt(1,2) = dtt(1,2) + vx*vy
               dtt(1,3) = dtt(1,3) + vx*vz
               dtt(2,2) = dtt(2,2) + vy*vy
               dtt(2,3) = dtt(2,3) + vy*vz
               dtt(3,3) = dtt(3,3) + vz*vz
            endif
         endif
      enddo
c      
c++     Subtract streaming motion to total tensor: 

         vxm = vxm / count 
         vym = vym / count 
         vzm = vzm / count 

         dtt(1,1) = dtt(1,1)/count - vxm*vxm 
         dtt(2,2) = dtt(2,2)/count - vym*vym 
         dtt(3,3) = dtt(3,3)/count - vzm*vzm 

         dtt(1,2) = dtt(1,2)/count - vxm*vym 
         dtt(1,3) = dtt(1,3)/count - vxm*vzm 

         dtt(2,3) = dtt(2,3)/count - vym*vzm 
      
         dtt(2,1) = dtt(1,2)
         dtt(3,1) = dtt(1,3)
         dtt(3,2) = dtt(2,3)
c      
      if (count .eq. 0.d0) count = 1.d0
c
      do i = 1,3
         do j = 1,3
            tt(i,j) = dtt(i,j)  
         enddo
      enddo
c
c     Compute eigenvalues and eigenvectors of tt :
c     ============================================
c
c     using Num.Recipes routine:
c
      call jacobi (tt,3,3,d,v,nrot)
c
c     eigenvalues are stored in d(1),d(2),d(3)
c     the corresponding eigenvectors are in 
c     d(1) --> v(1,1),v(2,1),v(3,1)
c     d(2) --> v(1,2),v(2,2),v(3,2)
c     d(3) --> v(1,3),v(2,3),v(3,3)
c
c     Write data in file xxx-gyy.DTENS :
c     ==================================
c

! Convert back to units in case it is required - vectors v taken to be velocities.

       do i = 1,3
          d(i) = d(i) / scv/scv
         do j = 1,3
            v(i,j) = v(i,j) / scv 
         enddo
      enddo

      call createna (fname,lfname,ign,ierr)
      call makeext  (lfname,lfname1,'DTENS')
c      
      inquire (file=lfname1,exist=ex)
c      
      if (ex .eqv. .true.) then ! file already exists
         open (1,file=lfname1,access='direct',recl=lori,status='old')
         read (1,rec = 1) fpos 
         ipos  = int(fpos)
         ihlen = hlend 
         close (1)
         open (1,file=lfname1,access='direct',recl=ihlen*lori,
     $        status='old')
      else                      ! we have to create a new file
         ihlen = hlend
*
*       Initialiase arary position; store conversion factors : 
         ipos   = 1
         hoh(1) = 1. 
         hoh(2) = float( hlend ) 
         hoh(3) = float( ih(5) ) 
         hoh(4) = fh(54,1) ! = gp 

         hoh(5) = fh(55,1) ! = g/gp 
         hoh(6) = fh(56,1) ! = scale of mass 
         hoh(7) = fh(57,1) ! = scale of length 
         hoh(8) = fh(58,1) ! = scale of time 
         hoh(9) = fh(59,1) ! = scale of velocity 
         hoh(10)= fh(60,1) ! = conversion factor km/sec -> pc or Kpc /Myrs 
         hoh(11)= float(ih(13)) ! choice of kpc or pc as unit of length. 

         do j = 12,ihlen 
            hoh(j) = 0. 
         end do 

         open (1,file=lfname1,access='direct',recl=ihlen*lori,
     $        status='new')
         write (1,rec = 1) hoh
      endif
c      
      hoh(1) = time
      hoh(2) = real(step)
      hoh(3) = rmax
      hoh(4) = count
      hoh(5) = real(iflag)
      do j = 1,3
         hoh(5+j) = d(j)
      enddo
      do i = 1,3
         do j = 1,3
            hoh(8+(i-1)*3+j) = v(i,j)
         enddo
      enddo
c     
      ipos = ipos + 1
      write (1,rec = ipos) hoh
c      
      close(1) 
c      
* Reopen to store only first the updated number of  entries : 

      open (1,file=lfname1,access='direct',recl=lori,
     $        status='old') 
      write (1,rec = 1) real(ipos) 
      close (1)
c      
      return
c
      end

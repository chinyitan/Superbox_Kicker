c**********************************************************************
c
c     E V A L . F
c
c**********************************************************************
c
c This subroutine calculates moment of inertia wrt centre of mass inside rmax :
c =============================================================================
c
      subroutine ritens (step,ign,rmax,time,iflag)
c
      include 'super.cb'
c
      integer i,j,ierr, indx(3), shape_iter 
      integer step,ign,start,istop,nrot,ipos,ihlen,iflag
      real    tt(3,3),d(3),v(3,3),rmax,count,hoh(hleni+6),time
      real    scl1,fhlen,fpos, v_in(3,3), a,b,c,dd, vprod,a1,b1,cc1  
      real    nstar(6,numstars), sum1, sum2, sum3, vol, delta 
      real    vrot(3,3), v1(3), v2(3) , v3(3), phi, psi, teta, d1,d2,d3  
      real    subcm(3)  
      double precision  dtt(3,3),mstar
      double precision  dx,dy,dz,r,rmax2
c
      logical ex
      character*100 lfname,lfname1
      external ludcmp, lubksb 
c
      do i = 1,3
         do j = 1,3 

            dtt(i,j) = 0.d0
            v(i,j) = 0. 
            v_in(i,j) = 0. 
            vrot(i,j) = 0.
         enddo
            vrot(i,i) = 1.
            subcm(i)  = 0. 
      enddo
c ign = galaxy which is considered here, e.g. galaxy#3, ign=3
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

c++     Transfer particle data to new matrix - safety procedure cmb 26.07.1999
c       fh(..) = center of mass -> gcsm(1-3,ign) 
      do i = start, istop 
         do j = 1,3 
         nstar(j,i) =  star(j,i) - gcms(j,ign)
         end do 
      end do 

c++     Add small value to maximum radius to fight off roundoff errors.

      rmax2 = dble(rmax)*dble(rmax) + 2.*epsilon*rmax 
    
      shape_iter = 0 

c++     Initialise volume to be scanned : spherical density profile. 

      a = rmax 
      b = a 
      c = a 

90210 continue
  
c++     Compute the centre of mass of the stars selected ONLY 

      count = 0. 
     
      do i = 1,3 
         subcm(i) = 0.
      end do 

      do i = start,istop

         if (abs(nstar(1,i)) .lt. 1.e10 .and. abs(nstar(2,i)) .lt. 1.e10 
     &       .and. abs(nstar(3,i)) .lt. 1.e10 ) then
             dx = nstar(1,i)  
             dy = nstar(2,i) 
             dz = nstar(3,i) 
             r  = dx*dx/a/a + dy*dy/b/b + dz*dz/c/c
            if (r .lt. 1.0) then
             count = count + 1 
             subcm(1) = subcm(1) + dx  
             subcm(2) = subcm(2) + dy 
             subcm(3) = subcm(3) + dz 
            end if 
         end if 
      end do 
      
      do i = 1,3 
         subcm(i) = subcm(i) / count 
      end do 

c++     Reset particle positions according to centre of mass of selected particles.

      do i = start, istop 
         do j = 1,3 
         nstar(j,i) =  star(j,i) - subcm(j)
         end do 
      end do 

      count = 0 

      do i = start,istop

         if (abs(nstar(1,i)) .lt. 1.e10 .and. abs(nstar(2,i)) .lt. 1.e10 
     &       .and. abs(nstar(3,i)) .lt. 1.e10 ) then
             dx = nstar(1,i)  
             dy = nstar(2,i)  
             dz = nstar(3,i) 
             r  = dx*dx/a/a + dy*dy/b/b + dz*dz/c/c
            if (r .lt. 1.0) then

             count    = count + 1.d0
C See Landau&Lifschitz, 11th ed. (1984), p.121 : F.Scheck, "Mechanik", Springer Verlag, p.112
             r  = dx*dx + dy*dy + dz*dz
             dtt(1,1) = dtt(1,1) + r - dx*dx
             dtt(2,2) = dtt(2,2) + r - dy*dy
             dtt(3,3) = dtt(3,3) + r - dz*dz 
             dtt(1,2) = dtt(1,2) - dx*dy
             dtt(1,3) = dtt(1,3) - dx*dz
             dtt(2,3) = dtt(2,3) - dy*dz
            endif
         endif
      enddo
c Use symmetry argument to complete the matrix . 
      dtt(2,1) = dtt(1,2)
      dtt(3,1) = dtt(1,3)
      dtt(3,2) = dtt(2,3)
c Mass of a single star:  [gstno(ign) = total number of stars in galaxy ign]
      mstar = mtot(ign) / real( gstno(ign) )
      do i = 1,3
         do j = 1,3
            tt(i,j) = dtt(i,j)*mstar
            dtt(i,j) = 0. 
         enddo
      enddo

c     Compute eigenvalues and eigenvectors of dtt :
c     =============================================
c
c     using Num.Recipes routine:
     
      call jacobi (tt,3,3,d,v,nrot)

c     eigenvalues are stored in d(1),d(2),d(3)
c     the corresponding eigenvectors are in 
c     d(1) --> v(1,1),v(2,1),v(3,1)
c     d(2) --> v(1,2),v(2,2),v(3,2)
c     d(3) --> v(1,3),v(2,3),v(3,3)

!    Transform coordinates to new frame, ie frame of eigenvectors. 

c++  Firstly, normalise eigenvectors to unity: check only, as it should be
c    done in Jacobi. (v1 x v2) . v3 = 1 (works, 27.07.1999) 

       vprod = ( v(1,1)*v(2,2)-v(2,1)*v(1,2) ) * v(3,3) - 
     c         ( v(1,1)*v(3,2)-v(3,1)*v(1,2) ) * v(2,3) - 
     c         ( v(3,1)*v(2,2)-v(2,1)*v(3,2) ) * v(1,3)

c++       write( 6,* ) ' In itens - vprod = ', vprod 

c    Secondly, compute the inverse matrix 

       do i = 1,3 
          do j = 1,3 

             v_in(i,j) = 0. 
          end do 

             v_in(i,i) = 1. 
       end do 

       call ludcmp(v,3,3,indx,dd) 
c
       do j = 1,3 
          call lubksb( v,3,3,indx,v_in(1,j) ) 
       end do 

c++  Firstly, normalise eigenvectors to unity: check only, as it should be
c    done in Jacobi. (v1 x v2) . v3 = 1 (works, 27.07.1999)

       vprod = ( v_in(1,1)*v_in(2,2)-v_in(2,1)*v_in(1,2) ) * v_in(3,3) - 
     c         ( v_in(1,1)*v_in(3,2)-v_in(3,1)*v_in(1,2) ) * v_in(2,3) - 
     c         ( v_in(3,1)*v_in(2,2)-v_in(2,1)*v_in(3,2) ) * v_in(1,3)

c++       write( 6,* ) ' In itens - vprod = ', vprod 

c++     3rdly, Transform data to new vector coordinates. note the same would 
c       apply to the velocity field. 

       do i = start, istop 

          sum1 = v_in(1,1)*nstar(1,i)+v_in(1,2)*nstar(2,i)+ 
     &                v_in(1,3)*nstar(3,i) 
          sum2 = v_in(2,1)*nstar(1,i)+v_in(2,2)*nstar(2,i)+ 
     &                v_in(2,3)*nstar(3,i)      
          sum3 = v_in(3,1)*nstar(1,i)+v_in(3,2)*nstar(2,i)+ 
     &                v_in(3,3)*nstar(3,i)   

          nstar(1,i) = sum1 
          nstar(2,i) = sum2
          nstar(3,i) = sum3 
       end do 

! Build rotation matrix : 

       do j = 1,3
 
       sum1 =v_in(1,1)*vrot(1,j)+v_in(1,2)*vrot(2,j)+v_in(1,3)*vrot(3,j)
       sum2 =v_in(2,1)*vrot(1,j)+v_in(2,2)*vrot(2,j)+v_in(2,3)*vrot(3,j)
       sum3 =v_in(3,1)*vrot(1,j)+v_in(3,2)*vrot(2,j)+v_in(3,3)*vrot(3,j)

       vrot(1,j) = sum1 
       vrot(2,j) = sum2 
       vrot(3,j) = sum3 

       end do 

! Iterative scheme: from the eigenvalues, work out quadratic-weighted 
! axes a,b,c: a>b>c. 

       a1 = a 
       b1 = b 
       cc1 = c

       a = sqrt( d(3) + d(2) - d(1) ) / 1.4142 /sqrt(mstar*count)
       b = sqrt( d(3) + d(1) - d(2) ) / 1.4142 /sqrt(mstar*count)
       c = sqrt( d(2) + d(1) - d(3) ) / 1.4142 /sqrt(mstar*count)

c++     Compare volumes spanned: make sure volume spanned is unchanged

       vol = rmax**3 / a/b/c 
       a = a * vol**0.3333
       b = b * vol**0.3333
       c = c * vol**0.3333 

c++     Error estimate on the changes to a,b,c : 

       delta = sqrt( (a-a1)**2+(b-b1)**2+(c-cc1)**2 )/(a+b+c)

c++     Define volume for triaxial structure, iterate : 

       shape_iter = shape_iter + 1 

      if( shape_iter .lt. 5 .and. 
     &    delta .gt. 1./sqrt(count) ) go to 90210 

c++    Re-organise vectors & eigenvalues in increasing order of the principal axes:

       a1 = max( a, max(b,c) ) 
       cc1 = min( a, min(b,c) ) 
       if( b.ne.a1 .and. c.ne. a1 ) b1 = min( a, max(b,c) )
       if( b.eq. a1 .or. c.eq. a1 ) b1 = max( a, min(b,c) ) 

       do j = 1,3 

       if( a1 .eq. a ) then 

          v1(j) = vrot(j,1) 
          d1 = d(1) 
          if( b1 .eq. b ) then 
             v2(j) = vrot(j,2)
             v3(j) = vrot(j,3)
             d2 = d(2) 
             d3 = d(3) 
          else 
             v2(j) = vrot(j,3) 
             v3(j) = vrot(j,2) 
             d2 = d(3) 
             d3 = d(2) 
          endif 
       end if 

       if( a1 .eq. b ) then 
          v1(j) = vrot(j,2) 
          d1 = d(2) 
          if( b1 .eq. a ) then 
             v2(j) = vrot(j,1) 
             v3(j) = vrot(j,3) 
             d2 = d(1) 
             d3 = d(3) 
          else 
             v2(j) = vrot(j,3)
             v3(j) = vrot(j,1)
             d2 = d(3) 
             d3 = d(1) 
          endif 
       endif 

       if( a1 .eq. c ) then 
          v1(j) = vrot(j,3)
          d1 = d(3) 
          if( b1 .eq. a ) then 
             v2(j) = vrot(j,1)
             v3(j) = vrot(j,2) 
             d2 = d(1) 
             d3 = d(2) 
          else 
             v2(j) = vrot(j,2)
             v3(j) = vrot(j,1) 
             d2 = d(2) 
             d3 = d(1) 
          endif 
       endif 

       end do 

c++     Re-shuffle things ... 

       d(1) = d1 
       d(2) = d2 
       d(3) = d3 

c++     Compute the three Euler angles of the final transformation, cf Goldstein p.146
c++     Actually, it is sufficient to compute the angles of the major axis alone (the 
c++     other two lie in a plane orthogonal to it. Angles: theta, phi. 

       teta  = acos( v1(3) ) 
       phi   = acos( v1(1)/sqrt( v1(1)**2 + v1(2)**2 ) ) 
       psi   = 0. 

c     Save data in file modelname-gnr.ITENS :
c     ======================================

! Convert back to units in case it is required -vectors v taken to be lengths. 

       do i = 1,3
          d(i) = d(i) /scm / scl/scl 

            v(i,1) = v1(i) / scl
            v(i,2) = v2(i) / scl 
            v(i,3) = v3(i) / scl 
      enddo
c     
      call createna (fname,lfname,ign,ierr)
      call makeext  (lfname,lfname1,'ITENS')
c
      inquire (file=lfname1,exist=ex)
c
      if (ex .eqv. .true.) then ! file already exists
         open (1,file=lfname1,access='direct',recl=lori,status='old')
         read (1,rec = 1) fpos
         ipos  = int(fpos)
         ihlen = hleni + 6 
         close (1)
         open (1,file=lfname1,access='direct',recl=ihlen*lori,
     $        status='old')
      else                    ! we have to create a new file
         ihlen = hleni + 6 
         open (1,file=lfname1,access='direct',recl=ihlen*lori,
     $        status='new')
*
*       Initialiase array position; store conversion factors : 
         ipos   = 1
         hoh(1) = 1. 
         hoh(2) = float( hleni ) 
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
         write (1,rec = 1) hoh
      endif
c      
      hoh(1) = time/sct
      hoh(2) = real(step)
      hoh(3) = rmax/scl
      hoh(4) = count
      hoh(5) = real(iflag)
      do j = 1,3
         hoh(5+j) = d(j)
      enddo

         hoh(9) = a1
         hoh(10)= b1 
         hoh(11)= cc1 

         hoh(12)= teta
         hoh(13)= phi
         hoh(14)= psi 

      do i = 1,3
         do j = 1,3
C hoh(9) =v(1,1), hoh(10)=v(1,2), hoh(11)=v(1,3)
C hoh(12)=v(2,1), hoh(13)=v(2,2), hoh(14)=v(2,3)
C hoh(15)=v(3,1), hoh(16)=v(3,2), hoh(17)=v(3,3)
            hoh(14+(i-1)*3+j) = v(i,j)
         enddo
      enddo
c      
      ipos = ipos + 1
      write (1,rec = ipos) hoh         

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


c      
c******************************************************************
c      
c     updates the positions of all stars and gets the output data
c
      subroutine pusher2 (step)
c      
      include 'super.cb'
c      
      real    rmax,myr,sf4,fnh,x14,y14,z14
      real    rad(numstars),x,y,z,r
      integer i,j,k,ipr,ip
      integer step,rest,start,istop,ix,iy,iz
c      
      double precision  dbx,dby,dbz,dbvx,dbvy,dbvz,dbgcms(6)
c
      tmass = 0.0
      do j = 1,6
         tcms(j) = 0.0
      enddo
c
c     Start the loop over all galaxies :
c     ==================================
c
      do k = 1,gnum             ! starting galaxy loop
c
         if (k .eq. 1) then     ! getting the number of the first
            start = 1           ! and the last star of galaxy k
            istop = gstno(1)
         else
            start = 0
            do i = 1,k-1
               start = start + gstno(i)
            enddo
            istop = start + gstno(k)
            start = start + 1
         endif
c         
         rest   = 0
c
         do i = 1,6
            gcms(i,k)  = 0.0
            dbgcms(i)  = 0.d0
         enddo
c
         sf4    = enh4(k)
         fnh    = real(n/2)
c         
         do  i = start,istop    ! starting particle loop
c
c     Skip stars out of the simulation :
c     ==================================
c
            if (star(1,i) .gt. 1.e+10) goto 21
c
c     Change to double precission :
c     =============================
c
            dbx  = dble(star(1,i))
            dby  = dble(star(2,i))
            dbz  = dble(star(3,i))
            dbvx = dble(star(4,i))
            dbvy = dble(star(5,i))
            dbvz = dble(star(6,i))
c
c     Leap-Frog integration :
c     =======================
c
            dbx = dbx + dbvx * dble(dt(k))
            dby = dby + dbvy * dble(dt(k))
            dbz = dbz + dbvz * dble(dt(k))
c
c     Change back to single precision :
c     =================================
c
            star(1,i) = real(dbx)
            star(2,i) = real(dby)
            star(3,i) = real(dbz)
c            
c     Looking for stars that fall out of simulation :
c     ===============================================
c
            x14 = sf4 * star(1,i) + fnh
            y14 = sf4 * star(2,i) + fnh
            z14 = sf4 * star(3,i) + fnh
            ix = nint(x14)
            iy = nint(y14)
            iz = nint(z14)
c
            if ((ix .lt. 2) .or. (ix .ge. (n-2))) then
               star(1,i) = 1.e+30
               star(2,i) = 1.e+30
               star(3,i) = 1.e+30
               goto 21
            endif
            if ((iy .lt. 2) .or. (iy .ge. (n-2))) then
               star(1,i) = 1.e+30
               star(2,i) = 1.e+30
               star(3,i) = 1.e+30
               goto 21
            endif
            if ((iz .lt. 2) .or. (iz .ge. (n-2))) then
               star(1,i) = 1.e+30
               star(2,i) = 1.e+30
               star(3,i) = 1.e+30
               goto 21
            endif
c
            rest   = rest+1
c
c     Calculating new center of mass for all stars :
c     ==============================================
c
            do j = 1,6
               dbgcms(j) = dbgcms(j) + dble(star(j,i))
            enddo
c
 21         continue
c
         enddo                  ! end particle loop
c
         grest(k) = rest        ! particles left in simulation
c
         if (grest(k) .gt. 0) then
            do j = 1,6
               gcms(j,k)  = real( dbgcms(j) / dble(grest(k)) )
            enddo
         endif
         
c     Calculating center of density :
c     ===============================
c
c     with newcms if focus is on center of density, else center of mass
c     is taken as focus
c
         if (origin.eq.0) then
            do j = 1,6
               dgcms(j,k) = gcms(j,k)
            enddo
         else
            call newcms (step,k,start,istop)
         endif
c         
         rest   = 0
c
c     Now the radii of all stars of galaxy k are calculated :
c     =======================================================
c
         if (gstno(k).gt.0) then
c     
            do i = start,istop  ! starting 2nd particle loop
c     
               if (star(1,i).gt.1.e+10) goto 22
c
               x = star(1,i) - dgcms(1,k)
               y = star(2,i) - dgcms(2,k)
               z = star(3,i) - dgcms(3,k)
               rest = rest+1
               r = x*x + y*y + z*z
               rad(rest) = r
c               
 22            continue
c
            enddo               ! end 2nd particle loop
c
            call sort (rest,rad) 
c            
            ipr   = rest / 10
c
c     Calculating 10%...90% Lagrange-radii :
c     ======================================
c
            do j = 1,9
               outhead(32+j,k) = sqrt(rad(j*ipr)) / scl
            enddo

c
         endif
c     
c     Calculate new center of simulation :
c     ====================================
c
         tmass   = tmass + mass(k) * grest(k)
         do j = 1,6
            tcms(j) = tcms(j) + gcms(j,k) * mass(k) * grest(k)
         enddo
c
      enddo                     ! end galaxy loop
c
      do j = 1,6
         tcms(j) = tcms(j) / tmass
      enddo
c      
c     Now all the output data is collected :
c     ======================================
c
      ih(3) = step
c      
      do k = 1,gnum             ! start 2nd galaxy loop
c
         do i = 1,3
            fh(7+i,k)  =  gcms(i,k)  / scl
            fh(13+i,k) =  dgcms(i,k) / scl
         enddo
         do i = 4,6
            fh(7+i,k)  =  gcms(i,k)  / scv
            fh(13+i,k) =  dgcms(i,k) / scv
         enddo
c
         ctime(k) = ctime(k) + dt(k)
         fh(20,k) = ctime(k) / sct
c
         fh(7,k)  = grest(k)
c         
c     Filling up output array outhead:
c     ================================
c
c     length of array is 60 (hlen) * real per galaxy
c
c     outhead (01) : stepnumber
c     outhead (02) : length of timestep
c     outhead (03) : time of simulation
c     outhead (04-09) : center of mass and velocity (all)
c     outhead (10-15) : center of density
c     outhead (16-21) : center of the simulation
c     outhead (22) : stars left in simulation
c     outhead (23) : kinetic energy of galaxy
c     outhead (24) : potential energy of galaxy
c     outhead (25) : total kinetic energie of simulation
c     outhead (26) : total potential energy of simulation
c     outhead (27-29) : internal angular momentum of galaxy 
c     outhead (30-32) : total angular momentum
c     outhead (33-41) : radii of mass-shells (10-90 %) 
c     outhead (42-44) : hcms
c     outhead (45-60) : free for some more data (e.g. sticky particle)
c
         outhead(1,k) = real(step)
         outhead(2,k) = dt(k) /sct 
         outhead(3,k) = ctime(k) /sct
c         
         do i = 1,3
            outhead(3+i,k)  = gcms(i,k)  / scl
            outhead(9+i,k)  = dgcms(i,k) / scl
            outhead(15+i,k) = tcms(i) / scl
         enddo
         do i = 4,6
            outhead(3+i,k)  = gcms(i,k)  / scv
            outhead(9+i,k)  = dgcms(i,k) / scv
            outhead(15+i,k) = tcms(i) / scv
            outhead(39+i,k) = hcms(i,k) / scv
         enddo
c         
         outhead(22,k) = grest(k)
c         
         outhead(23,k) = gekin(k)
         outhead(24,k) = gepot(k)
         outhead(25,k) = ekin
         outhead(26,k) = epot
c
         outhead(27,k) = lgx(k)
         outhead(28,k) = lgy(k)
         outhead(29,k) = lgz(k)
         outhead(30,k) = lx
         outhead(31,k) = ly
         outhead(32,k) = lz
c         
c     If tifreq(k) > 0 calculate moment of inertia :
c     ==============================================
c
         if (tifreq(k) .gt. 0) then
            if (mod(step,tifreq(k)) .eq. 0) then
               if (frmaxi(k) .lt. 2) then
                  rmax = rmaxi(k)
               else
                  ip = int((rmaxi(k)/100.0)*gstno(k))
                  rmax = sqrt(rad(ip)) 
               endif 
               call ritens (step,k,rmax,ctime(k),frmaxi(k))
            endif 
         endif
c         
c     If tdfreq(k) > 0 calculate velocity dispersion :
c     ================================================
c
         if (tdfreq(k) .gt. 0) then
            if (mod (step,tdfreq(k)) .eq. 0) then
               if (frmaxd(k) .lt. 2) then
                  rmax = rmaxd(k)
               else
                  ip = int((rmaxd(k)/100.0)*gstno(k))
                  rmax = sqrt(rad(ip)) 
               endif 
               call disptens (step,k,rmax,ctime(k),frmaxd(k))
            endif 
         endif
c
      enddo                     ! end 2nd galaxy loop
c      
      call outdata (gnum)       ! writing outhead in xxx.HEAD
      return
c
      end

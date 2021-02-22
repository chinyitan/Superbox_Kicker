c**********************************************************************
c
c     G E T P X X . REC. F
c     Version komplett umgestellt auf neue fft - same file as getpxx.f
c     but now using the new numerical recipes fast fourier transform 
c     call to REALFT. must be compiled with GETH.RECIPE.F 
c
c**********************************************************************
c
      subroutine getpxx (grnr)
c
      include 'super.cb'
c
      real    rh(0:n-1,0:2*n-1,0:2*n-1),daten(2*n),z(0:2*n-1)
      real    odnm2
      integer grnr,i,j,k, nm2 
c
      odnm2 = 1.0 / real(n)
      nm2   = 2*n 
c
      do i = 0,n-1              ! for all x lines
         do j = 0,n-1           ! for all y lines
c
            do k = 0,n-1        ! load z-line into daten
               daten(k+1)   = rho(i,j,k,grnr)
               daten(n+k+1) = 0.0 ! and fill inactive part with 0
            enddo
c
            call realft (daten,nm2,1) ! 1-D-FFT of array daten
c
c     write array rh with Fourier-transform
c
            rh(i,j,0) = daten(1) ! extract the 2 real valued 
            rh(i,j,n) = daten(2) ! frequencies
            do k = 1,n-1
               rh(i,j,  k) = daten(2*k+1) ! place real-values
                                ! from 0 to n-1
               rh(i,j,n+k) = daten(2*k+2) ! place imag.-values 
                                ! from n+1 to nm2-1
            enddo
c
         enddo
      enddo
c
c     Fourier-Transformation in y-Direction:
c     --------------------------------------
c
      do i = 0,n-1              ! for all x lines
         do k = 0,nm2-1         ! for all z lines (now 2*n)
c
            do j = 0,n-1          ! load y-line into daten
               daten(j+1)   = rh(i,j,k)
               daten(n+j+1) = 0.0 ! and fill inactive part with 0
            enddo
c
            call realft (daten,nm2,1) ! 1-D-FFT of array daten
c
c     write array rh with Fourier-transform
c
            rh(i,0,k) = daten(1) ! see above
            rh(i,n,k) = daten(2)
            do j = 1,n-1
               rh(i,  j,k) = daten(2*j+1)
               rh(i,n+j,k) = daten(2*j+2)
            enddo
c
         enddo
      enddo
c
c     Fourier-Transformation in x-Direction:
c     --------------------------------------
c
      do j = 0,nm2-1              ! for all y lines (now 2*n)
         do k = 0,nm2-1           ! for all z lines (now 2*n)
c
            do i = 0,n-1          ! load x-line into daten
               daten(i+1)   = rh(i,j,k)
               daten(n+i+1) = 0.0 ! fill up with 0
            enddo
c
            call realft (daten,nm2,1) ! 1-D-FFT of array daten
c
c     write array Fourier-transform in array z
c
            z(0) = daten(1) ! see above
            z(n) = daten(2)
            do i = 1,n-1
               z(  i) = daten(2*i+1)
               z(n+i) = daten(2*i+2)
            enddo
c
c     Multiplication with Green's Function h:
c     ----------------------------------------
c
c     this is done immediately within the outer loops (j,k)
c
            if ((j .gt. n) .and. (k .le. n)) then 
c
               z(0) = z(0) * h(0,j-n,k)
               z(n) = z(n) * h(n,j-n,k)

               do i = 1,n-1
                  z(  i) = z(  i) * h(i,j-n,k)
                  z(n+i) = z(n+i) * h(i,j-n,k)

               enddo
c
            else if ((j .le. n) .and. (k .gt. n)) then
c
               z(0) = z(0) * h(0,j,k-n)
               z(n) = z(n) * h(n,j,k-n)

               do i = 1,n-1
                  z(  i) = z(  i) * h(i,j,k-n)
                  z(n+i) = z(n+i) * h(i,j,k-n)

               enddo
c
            else if ((j .gt. n) .and. (k .gt. n)) then
c
               z(0) = z(0) * h(0,j-n,k-n)
               z(n) = z(n) * h(n,j-n,k-n)

               do i = 1,n-1
                  z(  i) = z(  i) * h(i,j-n,k-n)
                  z(n+i) = z(n+i) * h(i,j-n,k-n)

               enddo
c
            else                ! j and k are <= n
c
               z(0) = z(0) * h(0,j,k)
               z(n) = z(n) * h(n,j,k)
               do i = 1,n-1
                  z(  i) = z(  i) * h(  i,j,k)
                  z(n+i) = z(n+i) * h(i,j,k)

               enddo
c
            endif
c
c     Back-Transformation in x-Direction:
c     -----------------------------------
c
            daten(1) = z(0)
            daten(2) = z(n)
            do i = 1,n-1      ! load daten with
               daten(2*i+1) = z(  i) ! real parts
               daten(2*i+2) = z(n+i) ! imag. parts
            enddo               ! of rho
c
            call realft (daten,nm2,-1) ! Back-1-D-FFT of array daten
c
            do i = 0,n-1        ! Multiply with Fourier-Factor 2/n
               rh(i,j,k) = odnm2 * daten(i+1)
            enddo
c
         enddo
      enddo
c
c     Back-Transformation in y-Direction:
c     -----------------------------------
c
      do i = 0,n-1              ! for all x lines
         do k = 0,nm2-1         ! for all z lines
c
            daten(1) = rh(i,0,k)
            daten(2) = rh(i,n,k)
            do j = 1,n-1      ! load daten with
               daten(2*j+1) = rh(i,  j,k) ! real parts
               daten(2*j+2) = rh(i,n+j,k) ! imag. parts
            enddo               ! of rho
c
            call realft (daten,nm2,-1) ! Back-1-D-FFT of array daten
c
            do j = 0,n-1        ! Multiply with Fourier-Factor 2/n
               rh(i,j,k) = odnm2 * daten(j+1)
            enddo
c
         enddo
      enddo
c
c     Back-Transformation in z-Direction:
c     -----------------------------------
c
      do i = 0,n-1              ! for all x lines
         do j = 0,n-1           ! for all y lines
c
            daten(1) = rh(i,j,0)
            daten(2) = rh(i,j,n)
            do k = 1,n-1      ! load daten with
               daten(2*k+1) = rh(i,j,  k) ! real parts
               daten(2*k+2) = rh(i,j,n+k) ! imag. parts
            enddo               ! of rho
c
            call realft (daten,nm2,-1) ! Back-1-D-FFT of array daten
c
            do k = 0,n-1        ! Multiply with Fourier-Factor 2/n
               rho(i,j,k,grnr) = odnm2 * daten(k+1)
            enddo
c
         enddo
      enddo
c
      return
      end
c

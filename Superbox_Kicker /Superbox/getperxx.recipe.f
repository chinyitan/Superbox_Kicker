
c**********************************************************************
c
c     G E T PER  X X . REC. F
c     Version komplett umgestellt auf neue fft
c
c**********************************************************************
c
c	implements periodic cosmological boundary conditions 
c**********************************************************************
c
      subroutine getperxx
c
      include 'super.cb'
c
      integer i,j,k
      real    daten(n),tdn
c
      tdn = 2.0 / real(n)
c
c     Fourier-Transformation in z-Direction:
c     --------------------------------------
c
      do i = 0,n-1              ! for all x lines
         do j = 0,n-1           ! for all y lines
c
            do k = 0,n-1        ! load z-line into daten
               daten(k+1) = rho(i,j,k,6)
            enddo
c
            call realft (daten,n,1) ! 1-D-FFT of array daten
c
c     rewrite array rho with Fourier-transform
c
            rho(i,j,0  ,6) = daten(1) ! extract the 2 real valued 
            rho(i,j,nd2,6) = daten(2) ! frequencies
            do k = 1,nd2-1
               rho(i,j,    k,6) = daten(2*k+1) ! place real-values
                                ! from 0 to nd2
               rho(i,j,nd2+k,6) = daten(2*k+2) ! place imag.-values 
                                ! from nd2+1 to n-1
            enddo
c
         enddo
      enddo
c
c     Fourier-Transformation in y-Direction:
c     --------------------------------------
c
      do i = 0,n-1              ! for all x lines
         do k = 0,n-1           ! for all z lines
c
            do j = 0,n-1        ! load y-line into daten
               daten(j+1) = rho(i,j,k,6)
            enddo
c
            call realft (daten,n,1) ! 1-D-FFT of array daten
c
c     rewrite array rho with Fourier-transform
c
            rho(i,0  ,k,6) = daten(1) ! see above
            rho(i,nd2,k,6) = daten(2)
            do j = 1,nd2-1
               rho(i,    j,k,6) = daten(2*j+1)
               rho(i,nd2+j,k,6) = daten(2*j+2)
            enddo
c
         enddo
      enddo
c
c     Fourier-Transformation in x-Direction:
c     --------------------------------------
c
      do j = 0,n-1              ! for all y lines
         do k = 0,n-1           ! for all z lines
c
            do i = 0,n-1        ! load x-line into daten
               daten(i+1) = rho(i,j,k,6)
            enddo
c
            call realft (daten,n,1) ! 1-D-FFT of array daten
c
c     rewrite array rho with Fourier-transform
c
            rho(0  ,j,k,6) = daten(1) ! see above
            rho(nd2,j,k,6) = daten(2)
            do i = 1,nd2-1
               rho(    i,j,k,6) = daten(2*i+1)
               rho(nd2+i,j,k,6) = daten(2*i+2)
            enddo
c
c     Multiplication with Green's Function hp:
c     ----------------------------------------
c
c     this is done immediately within the outer loops (j,k)
c
c     only values of hp from 0 to nd2 are used
c
            if ((j .gt. nd2) .and. (k .le. nd2)) then 
                                ! we have to subtract nd2 from j
c
               rho(0  ,j,k,6) = rho(0  ,j,k,6) * hp(0  ,j-nd2,k)
               rho(nd2,j,k,6) = rho(nd2,j,k,6) * hp(nd2,j-nd2,k)
               do i = 1,nd2-1
                  rho(    i,j,k,6) = rho(    i,j,k,6) * hp(i,j-nd2,k)
                  rho(nd2+i,j,k,6) = rho(nd2+i,j,k,6) * hp(i,j-nd2,k)
               enddo
c
            else if ((j .le. nd2) .and. (k .gt. nd2)) then
                                ! we have to subtract nd2 from k
c
               rho(0  ,j,k,6) = rho(0  ,j,k,6) * hp(0  ,j,k-nd2)
               rho(nd2,j,k,6) = rho(nd2,j,k,6) * hp(nd2,j,k-nd2)
               do i = 1,nd2-1
                  rho(    i,j,k,6) = rho(    i,j,k,6) * hp(i,j,k-nd2)
                  rho(nd2+i,j,k,6) = rho(nd2+i,j,k,6) * hp(i,j,k-nd2)
               enddo
c
            else if ((j .gt. nd2) .and. (k .gt. nd2)) then
                                ! we have to subtract nd2 from j and k
c
               rho(0  ,j,k,6) = rho(0  ,j,k,6) * hp(0  ,j-nd2,k-nd2)
               rho(nd2,j,k,6) = rho(nd2,j,k,6) * hp(nd2,j-nd2,k-nd2)
               do i = 1,nd2-1
                  rho(    i,j,k,6) = rho(    i,j,k,6) 
     $                 * hp(i,j-nd2,k-nd2)
                  rho(nd2+i,j,k,6) = rho(nd2+i,j,k,6) 
     $                 * hp(i,j-nd2,k-nd2)
               enddo
c
            else                ! j and k are <= nd2
c
               rho(0  ,j,k,6) = rho(0  ,j,k,6) * hp(0  ,j,k)
               rho(nd2,j,k,6) = rho(nd2,j,k,6) * hp(nd2,j,k)
               do i = 1,nd2-1
                  rho(    i,j,k,6) = rho(    i,j,k,6) * hp(i,j,k)
                  rho(nd2+i,j,k,6) = rho(nd2+i,j,k,6) * hp(i,j,k)
               enddo
c
            endif
c
c     Back-Transformation in x-Direction:
c     -----------------------------------
c
            daten(1) = rho(0  ,j,k,6)
            daten(2) = rho(nd2,j,k,6)
            do i = 1,nd2-1      ! load daten with
               daten(2*i+1) = rho(    i,j,k,6) ! real parts
               daten(2*i+2) = rho(nd2+i,j,k,6) ! imag. parts
            enddo               ! of rho
c
            call realft (daten,n,-1) ! Back-1-D-FFT of array daten
c
            do i = 0,n-1        ! Multiply with Fourier-Factor 2/n
               rho(i,j,k,6) = tdn * daten(i+1) 
            enddo
c
         enddo
      enddo
c
c     Back-Transformation in y-Direction:
c     -----------------------------------
c
      do i = 0,n-1              ! for all x lines
         do k = 0,n-1           ! for all z lines
c
            daten(1) = rho(i,0  ,k,6)
            daten(2) = rho(i,nd2,k,6)
            do j = 1,nd2-1      ! load daten with
               daten(2*j+1) = rho(i,    j,k,6) ! real parts
               daten(2*j+2) = rho(i,nd2+j,k,6) ! imag. parts
            enddo               ! of rho
c
            call realft (daten,n,-1) ! Back-1-D-FFT of array daten
c
            do j = 0,n-1        ! Multiply with Fourier-Factor 2/n
               rho(i,j,k,6) = tdn * daten(j+1)
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
            daten(1) = rho(i,j,0  ,6)
            daten(2) = rho(i,j,nd2,6)
            do k = 1,nd2-1      ! load daten with
               daten(2*k+1) = rho(i,j,    k,6) ! real parts
               daten(2*k+2) = rho(i,j,nd2+k,6) ! imag. parts
            enddo               ! of rho
c
            call realft (daten,n,-1) ! Back-1-D-FFT of array daten
c
            do k = 0,n-1        ! Multiply with Fourier-Factor 2/n
               rho(i,j,k,6) = tdn * daten(k+1)
            enddo
c
         enddo
      enddo
c
      return
      end
c
c******************************************************************
c     end of file
c******************************************************************

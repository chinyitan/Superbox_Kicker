
c******************************************************************
c      
      subroutine set_scales( cnum,ih,fh ) 

c       In model units for SUPERBOX we have by default
c       Gs = 2 ; Ms = 1 = Rs   (s = superbox)  
c       I write the physical quantities G, M, R and t such that 
c       Gs = a x G  ; Ms = b x M ; Rs = c x R and ts = d x t 
c
c       The physical unit of mass is the solar mass, time is in years
c       and length has been defined above. Then by construction 
c         
c       [ts]^-2 = GsMs/Rs^3 = ab/c^3 [t]^-2 => [ts] = sqrt(c^3/(ab)) [t]
c       [Vs]^2  = GsMs/Rs   = ab/c GM/R   => [Vs] = sqrt(ab/c) [V] 
c                                                 = cv [V] 
c       Also 
c       [Ms] = ts^-2Rs^3/Gs = d^-2c^3/a  t^-2R/G => [Ms] = d^-2c^3/a [M]
c       and of course d = sqrt(c^3/(ab)).
c 
c      fh(54) = gp 
c      fh(55) = a 
c      fh(56) = b 
c      fh(57) = c
c      fh(58) = d 
c      fh(59) = cv = c/d
c      fh(60) = kmspcy 

      integer   galnum
      parameter (galnum=2)  ! max.number ofgalaxies

      integer      ih(60),gnum, unitl,cnum, model  
      real         fh(60,galnum),sl,ml,tl,vl
      real         a,b,c,cv,d,gp,kmspcy 

      character    rpc*4, munit*5, tunit*5

      common /scales/ sl,ml,tl,vl, gp,kmspcy,
     &                model,unitl,tunit,munit,rpc  

      unitl = ih(13) 
      tunit = 'Myrs '
      munit = 'solar'
      rpc = 'pc' 
      if( unitl .gt. 0.5 ) rpc = 'Kpc' 
      if( model .gt. 0.5 ) rpc = 'm.u.'
      if( rpc .eq. 'm.u.' ) tunit = rpc 
      if( rpc .eq. 'm.u.' ) munit = rpc 

      a = fh(55,cnum) 
      b = fh(56,cnum) 
      c = fh(57,cnum)
      d = fh(58,cnum)
      cv =fh(59,cnum) 
      gp =fh(54,cnum) 
      kmspcy = fh(60,cnum)

      sl = 1. 
      ml = 1. 
      tl = 1.
      vl = 1. 

      if( rpc .ne. 'm.u.' ) sl = c
      if( rpc .ne. 'm.u.' ) ml = b
      if( rpc .ne. 'm.u.' ) tl = d
      if( rpc .ne. 'm.u.' ) vl = cv   

      return 

      end 
